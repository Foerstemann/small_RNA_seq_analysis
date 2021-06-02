#!/bin/bash

set -u
set -o pipefail

# This script loops through all .fastq files in the current directory, 
# removes the adapter and size-selects the reads as given on the command line.

# It then performs the standard mapping and analysis tasks.

# It is recommended to create a sub-folder ../scripts/ in your data directory
# copy this script into it and run it from there
# calling syntax is then: bash scripts/small_RNA_seq_[date].sh [minlength] [maxlength]






# create logfile
echo "start">sRNA_processing.log
echo "start"
date>>sRNA_processing.log
minlength=$1
maxlength=$2

#start trimming
mkdir adapter_trimmed



echo "trimming adapter CTGTAGGCA and readlengths of $minlength nt to $maxlength nt. ">>sRNA_processing.log
echo "trimming adapter CTGTAGGCA and readlengths of $minlength nt to $maxlength nt. "


for xx in *.fastq; do
    echo >>sRNA_processing.log
    echo $xx>>sRNA_processing.log
    perl /home/exchange/sRNA_pipeline/process_3prime_rN_170406.plx -file:$xx -fq -adapter:CTGTAGGCA -path:adapter_trimmed >>sRNA_processing.log 2>>sRNA_processing.log

done
for xy in adapter_trimmed/*_trim3.fastq; do
    echo >>sRNA_processing.log
    echo $xy>>sRNA_processing.log

    perl /home/exchange/sRNA_pipeline/length_select_170406.plx -file:$xy -min:"$minlength" -max:"$maxlength" -fq >>sRNA_processing.log 2>>sRNA_processing.log

done

rm adapter_trimmed/*trim3.fastq

echo "trimming finished">>sRNA_processing.log 
echo "trimming finished"
date>>sRNA_processing.log 



# map reads with bowtie onto Drosophila genome
echo "start mapping reads with bowtie">>sRNA_processing.log
echo "start mapping reads with bowtie"
mkdir bowtie
mkdir SAM
mkdir garbage
for y in adapter_trimmed/*.fastq; do
    echo $y>>sRNA_processing.log
    echo $y
   # get basename
    filename=$(basename "$y")
    # get rid of extension  
    sample="${filename%.*}"
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster r5.33">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log -t -p4 --sam -v0 --un "garbage/${sample}.dump"  dme_chromosomes_533 "$y" "SAM/${sample}.sam"
    echo >>sRNA_processing.log 
#    echo "mapping to all cDNAs r5.33">>sRNA_processing.log 
#    bowtie 2>>sRNA_processing.log -t -p4 -v0 --un "garbage/${sample}.dump"  dmel-all-CDS-r5.33 "$y" "bowtie/${sample}_dmel-all-CDS-r5.33.map"
#    echo >>sRNA_processing.log  
#    echo "mapping to all transcripts r5.33">>sRNA_processing.log 
#    bowtie 2>>sRNA_processing.log -t -p4 -v0 --un "garbage/${sample}.dump"  dmel-all-transcript-r5.33 "$y" "bowtie/${sample}_dmel-all-transcript-r5.33.map"
#    echo >>sRNA_processing.log 
    echo "mapping to Drosophila melanogaster mature miRNAs, bowtie format">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log   -t -p4 -v0  maturemiRNAs_mirbase181022 "$y" "bowtie/${sample}_mature_miRNAs.map"
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster transposon consensus, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -v0  dme_transposon_cons "$y" "bowtie/${sample}_dme_transposons.map" 
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster mapping_loci_190627, bowtie format, all possible matches reported">>sRNA_processing.log     
    bowtie 2>>sRNA_processing.log  -t -p4 -v0 -a  mapping_loci_190627 "$y" "bowtie/${sample}_mapping_loci_190627.map" 
    echo >>sRNA_processing.log
    echo "mapping to blanks_dependent loci">>sRNA_processing.log
    bowtie 2>>sRNA_processing.log  -t -p4 -v0 blanks_dep_VN181018 "$y" "bowtie/${sample}_blanks_dep_VN181018.map"
    echo >>sRNA_processing.log
    echo "mapping to extended gene regions">>sRNA_processing.log
    bowtie 2>>sRNA_processing.log  -t -p4 -v0 dmel_extendedgenes_602 "$y" "bowtie/${sample}_extendedgenes_602.map"
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster tRNA collection, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -v0  dmel-all-tRNA-r5.33 "$y" "bowtie/${sample}_dmel-all-tRNA-r5.33.map" 
    echo >>sRNA_processing.log
    echo "mapping to E. coli genome, bowtie format">>sRNA_processing.log     
    bowtie 2>>sRNA_processing.log  -t -p4 --sam -v0  e_coli "$y" "SAM/${sample}_e_coli.sam_no_tdf" 
done


echo "mapping finished">>sRNA_processing.log
echo "mapping finished"
date>>sRNA_processing.log


# convert SAM files into BAM and BAI files
echo "start converting SAM files">>sRNA_processing.log
echo "start converting SAM files"
mkdir BAM
for z in SAM/*.sam; do
    echo $z>>sRNA_processing.log
     # get basename
    filename=$(basename "$z")
    # get rid of extension  
    sample="${filename%.*}"

	samtools view -bSF4 $z > "BAM/${sample}.bam"
	samtools sort "BAM/${sample}.bam" "BAM/${sample}"
	samtools index "BAM/${sample}.bam" "BAM/${sample}.bai"
done


# convert Bam Files to Bed files and to tdf files
mkdir Bed
mkdir tdf
for za in BAM/*.bam; do
    echo $za>>sRNA_processing.log
    echo $za
    # get basename
    filename=$(basename "$za")
    # get rid of extension  
    sample="${filename%.*}"
    # Bed-file
    bedtools bamtobed -i $za | sort -k1,1 -k2,2n > Bed/${sample}.sorted.bed
    bedtools merge -s -n -i Bed/${sample}.sorted.bed > Bed/${sample}.unique.bed
    # count genome matching reads for normalization
    wc -l  Bed/${sample}.sorted.bed >bowtie/${sample}.genome_matching	
	
    # tdf-file 
    /home/exchange/igvtools/igvtools count $za tdf/${sample}.tdf /home/exchange/igvtools/dmel_r5.33.chrom.sizes
done

# read in the number of genome matching reads
echo >>sRNA_processing.log
echo "----------------------" >>sRNA_processing.log
echo >> sRNA_processing.log
echo "reading in normalization values for ppm conversion">>sRNA_processing.log

if ls bowtie/*.genome_matching 1>/dev/null 2>&1
	then
	declare -A normalization_values;
	for a in bowtie/*.genome_matching; do
	basename=$(echo "$a" | cut -f 1 -d '.')
	normalization_values[$basename]=$(cat $a | awk '{print $1}')	
	echo "normalization read: $a"
	echo "$basename" >>sRNA_processing.log	
	echo "${normalization_values["$basename"]}">>sRNA_processing.log
	done
	echo >>sRNA_processing.log
	echo "----------------------">>sRNA_processing.log
	echo >>sRNA_processing.log
fi

# calculate reads per fasta line and size distribution
echo "treating mapping_loci_190627.fa">>sRNA_processing.log
echo "treating mapping_loci_190627.fa"
for b in bowtie/*_mapping_loci_190627.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_mapping_loci_190627.map//')
    echo $normname>>sRNA_processing.log
    echo $normname
 	normvalue=${normalization_values["$normname"]}
    echo $normvalue>>sRNA_processing.log
    echo $normvalue
    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_190627.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_190627.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_190627.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/loci_counts_ppm
mv bowtie/*.ppmcounts bowtie/loci_counts_ppm/

mkdir bowtie/binned_loci
mv bowtie/*.binned bowtie/binned_loci/

mkdir bowtie/binned_loci_ppm
mv bowtie/*.binppm bowtie/binned_loci_ppm/

mkdir bowtie/plotted_loci_ppm
mv bowtie/*_ppm.svg bowtie/plotted_loci_ppm/
mv bowtie/*_ppm.png bowtie/plotted_loci_ppm/

mkdir bowtie/plotted_loci
mv bowtie/*.svg bowtie/plotted_loci/
#mv bowtie/*.png bowtie/plotted_loci/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
#echo "treating dmel-all_CDS_r5.33.fa">>sRNA_processing.log
#echo "treating dmel-all_CDS_r5.33.fa"
#for b in bowtie/*_dmel-all-CDS-r5.33.map; do
#    echo >>sRNA_processing.log
#    echo $b>>sRNA_processing.log
#    echo $b#
#	normname=$(echo "$b" | sed 's/_dmel-all-CDS-r5.33.map//')
# 	normvalue=${normalization_values["$normname"]}#

#    echo "converting reads per fasta line">>sRNA_processing.log
#    echo "converting reads per fasta line"
#        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-CDS-r5.33.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

#    echo "mapping reads on bases per fasta line">>sRNA_processing.log
#    echo "mapping reads on bases per fasta line"
#    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/dmel-all-CDS-r5.33.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

#    echo "computing size distribution per fasta line">>sRNA_processing.log
#    echo "computing size distribution per fasta line"
#    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-CDS-r5.33.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
#done

# file cleanup
#mkdir bowtie/all_CDS_counts_ppm
#mv bowtie/*.ppmcounts bowtie/all_CDS_counts_ppm/

#mkdir bowtie/binned_all_CDS
#mv bowtie/*.binned bowtie/binned_all_CDS/

#mkdir bowtie/binned_all_CDS_ppm
#mv bowtie/*.binppm bowtie/binned_CDS_ppm/

#mkdir bowtie/plotted_all_CDS_ppm
#mv bowtie/*_ppm.svg bowtie/plotted_all_CDS_ppm/
#mv bowtie/*_ppm.png bowtie/plotted_all_CDS_ppm/

#mkdir bowtie/plotted_all_CDS
#mv bowtie/*.svg bowtie/plotted_all_CDS/
#mv bowtie/*.png bowtie/plotted_all_CDS/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
#echo "treating dmel-all-transcript-r5.33.fa">>sRNA_processing.log
#echo "treating dmel-all-transcript-r5.33.fa"
#for b in bowtie/*_dmel-all-transcript-r5.33.map; do
#    echo >>sRNA_processing.log
#    echo $b>>sRNA_processing.log
#    echo $b
#	normname=$(echo "$b" | sed 's/_dmel-all-transcript-r5.33.map//')
# 	normvalue=${normalization_values["$normname"]}
#
#    echo "converting reads per fasta line">>sRNA_processing.log
#    echo "converting reads per fasta line"
#        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-transcript-r5.33.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

#    echo "mapping reads on bases per fasta line">>sRNA_processing.log
#    echo "mapping reads on bases per fasta line"
#    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/dmel-all-transcript-r5.33.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

#    echo "computing size distribution per fasta line">>sRNA_processing.log
#    echo "computing size distribution per fasta line"
#    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-transcript-r5.33.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
#done

# file cleanup
#mkdir bowtie/all_transcript_counts_ppm
#mv bowtie/*.ppmcounts bowtie/all_transcript_counts_ppm/

#mkdir bowtie/binned_all_transcript
#mv bowtie/*.binned bowtie/binned_all_transcript/

#mkdir bowtie/binned_all_transcript_ppm
#mv bowtie/*.binppm bowtie/binned_all_transcript_ppm/

#mkdir bowtie/plotted_all_transcript_ppm
#mv bowtie/*_ppm.svg bowtie/plotted_all_transcript_ppm/
#mv bowtie/*_ppm.png bowtie/plotted_all_transcript_ppm/

#mkdir bowtie/plotted_all_transcript
#mv bowtie/*.svg bowtie/plotted_all_transcript/
#mv bowtie/*.png bowtie/plotted_all_trascript/

# calculate reads per fasta line and size distribution
echo "treating dmel-all_tRNA_r5.33.fa">>sRNA_processing.log
echo "treating dmel-all_tRNA_r5.33.fa"
for b in bowtie/*_dmel-all-tRNA-r5.33.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-tRNA-r5.33.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-tRNA-r5.33.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

#    echo "mapping reads on bases per fasta line">>sRNA_processing.log
#    echo "mapping reads on bases per fasta line"
#    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/dmel-all-CDS-r5.33.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-tRNA-r5.33.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_tRNA_counts_ppm
mv bowtie/*.ppmcounts bowtie/all_tRNA_counts_ppm/

#mkdir bowtie/binned_all_CDS
#mv bowtie/*.binned bowtie/binned_all_CDS/

#mkdir bowtie/binned_all_CDS_ppm
#mv bowtie/*.binppm bowtie/binned_CDS_ppm/

#mkdir bowtie/plotted_all_CDS_ppm
#mv bowtie/*_ppm.svg bowtie/plotted_all_CDS_ppm/
#mv bowtie/*_ppm.png bowtie/plotted_all_CDS_ppm/

#mkdir bowtie/plotted_all_CDS
#mv bowtie/*.svg bowtie/plotted_all_CDS/
#mv bowtie/*.png bowtie/plotted_all_CDS/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log
echo "treating D_mel_transposon_sequence_cons_set_v941.fasta">>sRNA_processing.log
echo "treating D_mel_transposon_sequence_cons_set_v941.fasta"
for c in bowtie/*_dme_transposons.map; do
    echo >>sRNA_processing.log
    echo $c>>sRNA_processing.log
    echo $c
	normname=$(echo "$c" | sed 's/_dme_transposons.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
    perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" -bin:1 -norm:"$normvalue" >>sRNA_processing.log  2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" >>sRNA_processing.log 2>>sRNA_processing.log

done

# file cleanup

mkdir bowtie/TE_counts_ppm
mv bowtie/*.ppmcounts bowtie/TE_counts_ppm/

mkdir bowtie/binned_TE
mv bowtie/*.binned bowtie/binned_TE/

mkdir bowtie/TE_ppm
mv bowtie/*.binppm bowtie/TE_ppm/

mkdir bowtie/plotted_TE_ppm
mv bowtie/*_ppm.svg bowtie/plotted_TE_ppm/
mv bowtie/*_ppm.png bowtie/plotted_TE_ppm/

mkdir bowtie/plotted_TE
mv bowtie/*.svg bowtie/plotted_TE/
#mv bowtie/*.png bowtie/plotted_TE/


echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log



echo "treating mature_miRNAs_mirbase181022.fasta">>sRNA_processing.log
echo "treating mature_miRNAs_mirbase181022.fasta"
for d in bowtie/*_mature_miRNAs.map; do
    echo >>sRNA_processing.log
    echo $d>>sRNA_processing.log
    echo $d
	# compute correct normalization parameter
	normname=$(echo "$d" | sed 's/_mature_miRNAs.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
    perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/mature_miRNAs_mirbase181022.fasta -map:"$d" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/mature_miRNAs_mirbase181022.fasta -map:"$d"  >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup

mkdir bowtie/miRNA_ppm_counts
mv bowtie/*.ppmcounts bowtie/miRNA_ppm_counts/


echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

echo "treating blanks_dep_VN181018.fa">>sRNA_processing.log
echo "treating blanks_dep_VN181018.fa"
for e in bowtie/*_blanks_dep_VN181018.map; do
    echo >>sRNA_processing.log
    echo $e>>sRNA_processing.log
    echo $e
	normname=$(echo "$e" | sed 's/_blanks_dep_VN181018.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
    perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/blanks_dep_VN181018.fa -map:"$e" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/blanks_dep_VN181018.fa -map:"$e" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log


    echo "computing size distribution per fasta line">>sRNA_processing.log 
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/blanks_dep_VN181018.fa -map:"$e" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup

mkdir bowtie/blanks_loci_counts_ppm
mv bowtie/*.ppmcounts bowtie/blanks_loci_counts_ppm/

mkdir bowtie/binned_blanks_loci
mv bowtie/*.binned bowtie/binned_blanks_loci/

mkdir bowtie/blanks_loci_ppm
mv bowtie/*.binppm bowtie/blanks_loci_ppm/

mkdir bowtie/plotted_blanks_loci_ppm
mv bowtie/*_ppm.svg bowtie/plotted_blanks_loci_ppm/
mv bowtie/*_ppm.png bowtie/plotted_blanks_loci_ppm/

mkdir bowtie/plotted_blanks_loci
mv bowtie/*.svg bowtie/plotted_blanks_loci/
#mv bowtie/*.png bowtie/plotted_blanks_loci/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

echo "mapfile analysis finished">>sRNA_processing.log
echo "mapfile analysis finished"
date>>sRNA_processing.log



mkdir bowtie/counts
mv bowtie/*.counts bowtie/counts

mkdir bowtie/sizes
mv bowtie/*.sizes bowtie/sizes




echo "file conversion finished">>sRNA_processing.log
echo "file conversion finished"
date>>sRNA_processing.log

mkdir logfiles
mv -f sRNA_processing.log logfiles/
