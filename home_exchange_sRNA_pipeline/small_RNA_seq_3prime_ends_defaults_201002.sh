#!/bin/bash

set -u
set -o pipefail

# This script loops through all .fastq files in the current directory, 
# removes the adapter and size-selects the reads as given on the command line.

# It then performs the standard mapping and analysis tasks.

# It is recommended to create a sub-folder ../scripts/ in your data directory
# copy this script into it and run it from there
# calling syntax is then: bash scripts/small_RNA_seq_[date].sh [minlength] [maxlength] [keep_files]
# all options must be given as integer numbers, only keep_files=0 will delete the less informative files


# some options can be set here, they should only be changed for good reason - therefore no command line query!
genome_version="dmel_chromosomes_602"
bowtie_mismatch=0
bowtie_trim5=0
bowtie_trim3=0
mapping_loci="mapping_loci_200920"

# create logfile
echo "start">sRNA_processing.log
echo "start"
date>>sRNA_processing.log

# read in command line options
minlength=$1
maxlength=$2
keep_files=$3

#start trimming
mkdir adapter_trimmed



echo "trimming adapter CTGTAGGCA and readlengths of $minlength nt to $maxlength nt. ">>sRNA_processing.log
echo "trimming adapter CTGTAGGCA and readlengths of $minlength nt to $maxlength nt. "


for xx in *.fastq; do
    echo >>sRNA_processing.log
    echo $xx>>sRNA_processing.log
    perl /home/exchange/sRNA_pipeline/process_3prime_rN_170406.plx -file:$xx -fq -adapter:CTGTAGGCA -path:adapter_trimmed >>sRNA_processing.log 2>>sRNA_processing.log &
done

wait

for xy in adapter_trimmed/*_trim3.fastq; do
    echo >>sRNA_processing.log
    echo $xy>>sRNA_processing.log

    perl /home/exchange/sRNA_pipeline/length_select_170406.plx -file:$xy -min:"$minlength" -max:"$maxlength" -fq >>sRNA_processing.log 2>>sRNA_processing.log &

done

wait

rm adapter_trimmed/*trim3.fastq

echo "trimming finished">>sRNA_processing.log 
echo "trimming finished"
date>>sRNA_processing.log 



# map reads with bowtie onto Drosophila genome
echo "start mapping reads with bowtie with $bowtie_mismatch mismatche(s), $bowtie_trim5 nt 5'-trim and $bowtie_trim3 nt 3'-trim of read.">>sRNA_processing.log
echo "start mapping reads with bowtie with $bowtie_mismatch mismatche(s), $bowtie_trim5 nt 5'-trim and $bowtie_trim3 nt 3'-trim of read."
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
     echo "mapping to $genome_version">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log -t -p4 --sam -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best --un "garbage/${sample}.dump"  "$genome_version" "$y" "SAM/${sample}.sam"
    echo >>sRNA_processing.log 
    echo "mapping to all cDNAs r6.19">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best --un "garbage/${sample}.dump"  dmel-all-CDS-r6.19 "$y" "bowtie/${sample}_dmel-all-CDS-r6.19.map"
    echo >>sRNA_processing.log  
    echo "mapping to all transcripts r6.19">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best --un "garbage/${sample}.dump"  dmel-all-transcript-r6.19 "$y" "bowtie/${sample}_dmel-all-transcript-r6.19.map"
    echo >>sRNA_processing.log 
    echo "mapping to Drosophila melanogaster mature miRNAs, bowtie format">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log   -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best  maturemiRNAs_mirbase181022 "$y" "bowtie/${sample}_mature_miRNAs.map"
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster transposon consensus, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best --un "adapter_trimmed/${sample}_no_TE.fq" dme_transposon_cons "$y" "bowtie/${sample}_dme_transposons.map" 
    echo >>sRNA_processing.log
    echo "mapping without TE matching reads to $genome_version">>sRNA_processing.log 
    bowtie 2>>sRNA_processing.log -t -p4 --sam -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best "$genome_version" "adapter_trimmed/${sample}_no_TE.fq" "SAM/${sample}_no_TE.sam"
    echo >>sRNA_processing.log 
    echo "mapping to Drosophila melanogaster "$mapping_loci", bowtie format, all possible matches reported">>sRNA_processing.log     
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best -a  "$mapping_loci" "$y" "bowtie/${sample}_$mapping_loci.map" 
    echo >>sRNA_processing.log
    echo "mapping to extended gene regions">>sRNA_processing.log
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best dmel_extendedgenes_602 "$y" "bowtie/${sample}_extendedgenes_602.map"
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster tRNA collection, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch"  dmel-all-tRNA-r6.19 "$y" "bowtie/${sample}_dmel-all-tRNA-r6.19.map" 
        echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster ncRNA collection r6.19, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch"  dmel-all-ncRNA-r6.19 "$y" "bowtie/${sample}_dmel-all-ncRNA-r6.19.map" 
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster intron collection r6.19, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch"  dmel-all-intron-r6.19 "$y" "bowtie/${sample}_dmel-all-intron-r6.19.map" 
    echo >>sRNA_processing.log
    echo "mapping to Drosophila melanogaster intergenic collection r6.19, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch"  dmel-all-intergenic-r6.19 "$y" "bowtie/${sample}_dmel-all-intergenic-r6.19.map" 
    echo >>sRNA_processing.log   
    echo "mapping to Drosophila melanogaster miscRNA collection r6.19, bowtie format">>sRNA_processing.log  
    bowtie 2>>sRNA_processing.log  -t -p4 -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch"  dmel-all-miscRNA-r6.19 "$y" "bowtie/${sample}_dmel-all-miscRNA-r6.19.map" 
    echo >>sRNA_processing.log

    echo >>sRNA_processing.log
    echo "mapping to E. coli genome, bowtie format">>sRNA_processing.log     
    bowtie 2>>sRNA_processing.log  -t -p4 --sam -5 "$bowtie_trim5" -3 "$bowtie_trim3" -v "$bowtie_mismatch" --best  e_coli "$y" "SAM/${sample}_e_coli.sam" 
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
	samtools sort "BAM/${sample}.bam" -o "BAM/${sample}.sorted.bam"
	mv "BAM/${sample}.sorted.bam" "BAM/${sample}.bam"
	samtools index "BAM/${sample}.bam" "BAM/${sample}.bam.bai" # this naming is needed for IGV
	# remove SAM files to free up space
	rm "SAM/${sample}.sam"
done


# convert Bam Files to Bed files (tdf no longer possible with free floating genome version)
mkdir Bed
# mkdir tdf
for za in BAM/*.bam; do
    echo $za>>sRNA_processing.log
    echo $za
    # get basename
    filename=$(basename "$za")
    # get rid of extension  
    sample="${filename%.*}"
    # Bed-file
    bedtools bamtobed -i $za | sort -k1,1 -k2,2n > Bed/${sample}.sorted.bed
    bedtools merge -s -c 1 -o count -i Bed/${sample}.sorted.bed > Bed/${sample}.unique.bed
    # count genome matching reads for normalization
    wc -l  Bed/${sample}.sorted.bed >bowtie/${sample}.genome_matching	

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
echo "treating $mapping_loci.fa">>sRNA_processing.log
echo "treating $mapping_loci.fa"
for b in bowtie/*_"$mapping_loci".map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_'$mapping_loci'.map//')
    echo $normname>>sRNA_processing.log
    echo $normname
 	normvalue=${normalization_values["$normname"]}
    echo $normvalue>>sRNA_processing.log
    echo $normvalue
    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/"$mapping_loci".fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/"$mapping_loci".fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/"$mapping_loci".fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup 
# we use xargs because files may be too numerous for direct use of mv
# as a precaution, this is done even in cases where mv currently would not be overloaded

mkdir bowtie/loci_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/loci_counts_ppm/

mkdir bowtie/binned_loci
echo bowtie/*.binned | xargs mv -t bowtie/binned_loci/

mkdir bowtie/binned_loci_ppm
echo bowtie/*.binppm | xargs mv -t bowtie/binned_loci_ppm/

mkdir bowtie/plotted_loci_ppm
echo bowtie/*_ppm.svg | xargs mv -t bowtie/plotted_loci_ppm/
echo bowtie/*_ppm.png | xargs mv -t bowtie/plotted_loci_ppm/

mkdir bowtie/plotted_loci
echo bowtie/*.svg | xargs mv -t bowtie/plotted_loci/
echo bowtie/*.png | xargs mv -t bowtie/plotted_loci/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
echo "treating dmel-all_CDS_r6.19.fa">>sRNA_processing.log
echo "treating dmel-all_CDS_r6.19.fa"
for b in bowtie/*_dmel-all-CDS-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-CDS-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-CDS-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-CDS-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_CDS_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_CDS_counts_ppm/


echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
echo "treating dmel-all-transcript-r6.19.fa">>sRNA_processing.log
echo "treating dmel-all-transcript-r6.19.fa"
for b in bowtie/*_dmel-all-transcript-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-transcript-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-transcript-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-transcript-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_transcript_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_transcript_counts_ppm/


# calculate reads per fasta line and size distribution
echo "treating dmel-all_tRNA_r6.19.fa">>sRNA_processing.log
echo "treating dmel-all_tRNA_r6.19.fa"
for b in bowtie/*_dmel-all-tRNA-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-tRNA-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-tRNA-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/dmel-all-tRNA-r6.19.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-tRNA-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_tRNA_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_tRNA_counts_ppm/

mkdir bowtie/binned_all_tRNA
echo bowtie/*.binned | xargs mv -t bowtie/binned_all_tRNA/
 
mkdir bowtie/binned_all_tRNA_ppm
echo bowtie/*.binppm | xargs mv -t bowtie/binned_all_tRNA_ppm/

mkdir bowtie/plotted_all_tRNA_ppm
echo bowtie/*_ppm.svg | xargs mv -t bowtie/plotted_all_tRNA_ppm/
echo bowtie/*_ppm.png | xargs mv -t bowtie/plotted_all_tRNA_ppm/

mkdir bowtie/plotted_all_tRNA
echo bowtie/*.svg | xargs mv -t bowtie/plotted_all_tRNA/
echo bowtie/*.png | xargs mv -t bowtie/plotted_all_tRNA/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
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
echo bowtie/*.ppmcounts | xargs mv -t bowtie/TE_counts_ppm/

mkdir bowtie/binned_TE
echo bowtie/*.binned | xargs mv -t bowtie/binned_TE/

mkdir bowtie/TE_ppm
echo bowtie/*.binppm | xargs mv -t bowtie/TE_ppm/

mkdir bowtie/plotted_TE_ppm
echo bowtie/*_ppm.svg | xargs mv -t bowtie/plotted_TE_ppm/
echo bowtie/*_ppm.png | xargs mv -t bowtie/plotted_TE_ppm/

mkdir bowtie/plotted_TE
echo bowtie/*.svg | xargs mv -t bowtie/plotted_TE/
echo bowtie/*.png | xargs mv -t bowtie/plotted_TE/


# calculate reads per fasta line and size distribution
echo "treating dmel-all-ncRNA-r6.19.fa">>sRNA_processing.log
echo "treating dmel-all-ncRNA-r6.19.fa"
for b in bowtie/*_dmel-all-ncRNA-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-ncRNA-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-ncRNA-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/dmel-all-ncRNA-r6.19.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-ncRNA-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_ncRNA_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_ncRNA_counts_ppm/

mkdir bowtie/binned_ncRNA
echo bowtie/*.binned | xargs mv -t bowtie/binned_ncRNA/

mkdir bowtie/all_ncRNA_ppm
echo bowtie/*.binppm | xargs mv -t bowtie/all_ncRNA_ppm/

mkdir bowtie/plotted_ncRNA_ppm
echo bowtie/*_ppm.svg | xargs mv -t bowtie/plotted_ncRNA_ppm/
echo bowtie/*_ppm.png | xargs mv -t bowtie/plotted_ncRNA_ppm/

mkdir bowtie/plotted_ncRNA
echo bowtie/*.svg | xargs mv -t bowtie/plotted_ncRNA/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
echo "treating dmel-all-intron-r6.19.fa">>sRNA_processing.log
echo "treating dmel-all-intron-r6.19.fa"
for b in bowtie/*_dmel-all-intron-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-intron-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-intron-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-intron-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_intron_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_intron_counts_ppm/


echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log

# calculate reads per fasta line and size distribution
echo "treating dmel-all-intergenic-r6.19.fa">>sRNA_processing.log
echo "treating dmel-all-intergenic-r6.19.fa"
for b in bowtie/*_dmel-all-intergenic-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-intergenic-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-intergenic-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log


    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-intergenic-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_intergenic_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_intergenic_counts_ppm/

echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log


# calculate reads per fasta line and size distribution
echo "treating dmel-all-miscRNA-r6.19.fa">>sRNA_processing.log
echo "treating dmel-all-miscRNA-r6.19.fa"
for b in bowtie/*_dmel-all-miscRNA-r6.19.map; do
    echo >>sRNA_processing.log
    echo $b>>sRNA_processing.log
    echo $b
	normname=$(echo "$b" | sed 's/_dmel-all-miscRNA-r6.19.map//')
 	normvalue=${normalization_values["$normname"]}

    echo "converting reads per fasta line">>sRNA_processing.log
    echo "converting reads per fasta line"
        perl /home/exchange/sRNA_pipeline/count_hits_to_fasta_lines_190524.plx -ref:/home/exchange/bowtie/genomes/dmel-all-miscRNA-r6.19.fa -map:"$b" -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "mapping reads on bases per fasta line">>sRNA_processing.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/sRNA_pipeline/bin_3prime_individual_fasta_lines_plot_svg_190529.plx -ref:/home/exchange/bowtie/genomes/dmel-all-misc-r6.19.fa -map:"$b" -bin:1 -norm:"$normvalue" >>sRNA_processing.log 2>>sRNA_processing.log

    echo "computing size distribution per fasta line">>sRNA_processing.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/sRNA_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dmel-all-miscRNA-r6.19.fa -map:"$b" >>sRNA_processing.log 2>>sRNA_processing.log
done

# file cleanup
mkdir bowtie/all_miscRNA_counts_ppm
echo bowtie/*.ppmcounts | xargs mv -t bowtie/all_miscRNA_counts_ppm/

mkdir bowtie/binned_miscRNA
echo bowtie/*.binned | xargs mv -t bowtie/binned_miscRNA/

mkdir bowtie/all_miscRNA_ppm
echo bowtie/*.binppm | xargs mv -t bowtie/all_miscRNA_ppm/

mkdir bowtie/plotted_miscRNA_ppm
echo bowtie/*_ppm.svg | xargs mv -t bowtie/plotted_miscRNA_ppm/
echo bowtie/*_ppm.png | xargs mv -t bowtie/plotted_miscRNA_ppm/

mkdir bowtie/plotted_miscRNA
echo bowtie/*.svg | xargs mv -t bowtie/plotted_miscRNA/

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
echo bowtie/*.ppmcounts | xargs mv -t bowtie/miRNA_ppm_counts/


echo >>sRNA_processing.log
echo "----------------------">>sRNA_processing.log
echo >>sRNA_processing.log


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

if [ "$keep_files" = "0" ] 
	then
	rm -r garbage/
	rm bowtie/*.map
	echo "removed garbage and  mapfiles no longer needed" >>sRNA_processing.log
fi


mkdir logfiles
mv -f sRNA_processing.log logfiles/
