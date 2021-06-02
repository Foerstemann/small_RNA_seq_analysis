#!/bin/bash
# This scripts maps the reads contained in a subdirectory adapter_trimmed/
# onto several loci. For the whole genome, a SAM-format output is enforced, 
# while for the others (transposons, miRNAs, model loci) bowtie-format is used.
# The latter files are then summarized by our PERL scripts and size distribution etc. is calculated.

set -u
set -o pipefail
mkdir logfiles
echo "start">sRNA_mapping.log
echo "start"
date >>sRNA_mapping.log
echo >>sRNA_mapping.log
echo "----------------------">>sRNA_mapping.log
echo >>sRNA_mapping.log 
# map reads with bowtie onto Drosophila genome
echo "start mapping reads with bowtie">>sRNA_mapping.log
echo "start mapping reads with bowtie"
mkdir Dme_533_SAM
mkdir bowtie
mkdir garbage

for y in adapter_trimmed/*.fastq; do
    echo >>sRNA_mapping.log
    echo $y>>sRNA_mapping.log
    echo $y
   # get basename
    filename=$(basename "$y")
    # get rid of extension  
    sample="${filename%.*}"
    echo >>sRNA_mapping.log
    echo "mapping to Drosophila melanogaster r5.33, SAM format">>sRNA_mapping.log 
    bowtie 2>>sRNA_mapping.log   -t -p4 --sam -v0 --un "garbage/${sample}.dump"  dme_chromosomes_533 "$y" "Dme_533_SAM/${sample}_dme533.sam"
    echo >>sRNA_mapping.log
    echo "mapping to Drosophila melanogaster mature miRNAs, bowtie format">>sRNA_mapping.log 
    bowtie 2>>sRNA_mapping.log   -t -p4 -v0  dme_miR_mature "$y" "bowtie/${sample}_dme_miR_mature.map"
    echo >>sRNA_mapping.log
    echo "mapping to Drosophila melanogaster transposon consensus, bowtie format">>sRNA_mapping.log  
    bowtie 2>>sRNA_mapping.log  -t -p4 -v0  dme_transposon_cons "$y" "bowtie/${sample}_dme_transposons.map" 
    echo >>sRNA_mapping.log
    echo "mapping to Drosophila melanogaster mapping_loci_171124, bowtie format">>sRNA_mapping.log     
bowtie 2>>sRNA_mapping.log  -t -p4 -v0  mapping_loci_171124 "$y" "bowtie/${sample}_mapping_loci_170406.map" 

done


echo "mapping finished">>sRNA_mapping.log
echo "mapping finished"
date>>sRNA_mapping.log
echo >>sRNA_mapping.log
echo "----------------------">>sRNA_mapping.log
echo >>sRNA_mapping.log

# calculate reads per fasta line and size distribution
echo "treating mapping_loci_171124.fa">>sRNA_mapping.log
echo "treating mapping_loci_171124.fa"
for b in bowtie/*_mapping_loci_170406.map; do
    echo >>sRNA_mapping.log
    echo $b>>sRNA_mapping.log
    echo $b

    echo "converting reads per fasta line">>sRNA_mapping.log
    echo "converting reads per fasta line"
        perl /home/exchange/scripts_pipeline/count_hits_to_fasta_lines_170419.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_171124.fa -map:"$b" >>sRNA_mapping.log 2>>sRNA_mapping.log
    echo "mapping reads on bases per fasta line">>sRNA_mapping.log
    echo "mapping reads on bases per fasta line"
        perl /home/exchange/scripts_pipeline/map_in_interval_170419.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_171124.fa -map:"$b" -bin:1 >>sRNA_mapping.log 2>>sRNA_mapping.log
    echo "computing size distribution per fasta line">>sRNA_mapping.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/scripts_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_171124.fa -map:"$b" >>sRNA_mapping.log 2>>sRNA_mapping.log
done

echo >>sRNA_mapping.log
echo "----------------------">>sRNA_mapping.log
echo >>sRNA_mapping.log

echo "treating D_mel_transposon_sequence_cons_set_v941.fasta">>sRNA_mapping.log
echo "treating D_mel_transposon_sequence_cons_set_v941.fasta"
for c in bowtie/*_dme_transposons.map; do
    echo >>sRNA_mapping.log
    echo $c>>sRNA_mapping.log
    echo $c

    echo "converting reads per fasta line">>sRNA_mapping.log
    echo "converting reads per fasta line"
    perl /home/exchange/scripts_pipeline/count_hits_to_fasta_lines_170419.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" >>sRNA_mapping.log 2>>sRNA_mapping.log
    echo "mapping reads on bases per fasta line">>sRNA_mapping.log
    echo "mapping reads on bases per fasta line"
    perl /home/exchange/scripts_pipeline/map_in_interval_170419.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" -bin:1 >>sRNA_mapping.log 2>>sRNA_mapping.log
    echo "computing size distribution per fasta line">>sRNA_mapping.log 
    echo "computing size distribution per fasta line"
    perl /home/exchange/scripts_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" >>sRNA_mapping.log 2>>sRNA_mapping.log
done

echo >>sRNA_mapping.log
echo "----------------------">>sRNA_mapping.log
echo >>sRNA_mapping.log

echo "treating dme_miR_mature_miRBase_170119_U_to_T_170123.fasta">>sRNA_mapping.log
echo "treating dme_miR_mature_miRBase_170119_U_to_T_170123.fasta"
for d in bowtie/*_dme_miR_mature.map; do
    echo >>sRNA_mapping.log
    echo $d>>sRNA_mapping.log
    echo $d

    echo "converting reads per fasta line">>sRNA_mapping.log
    echo "converting reads per fasta line"
    perl /home/exchange/scripts_pipeline/count_hits_to_fasta_lines_170419.plx -ref:/home/exchange/bowtie/genomes/dme_miR_mature_miRBase_170119_U_to_T_170123.fasta -map:"$d" >>sRNA_mapping.log 2>>sRNA_mapping.log
    echo "computing size distribution per fasta line">>sRNA_mapping.log
    echo "computing size distribution per fasta line"
    perl /home/exchange/scripts_pipeline/quantify_size_distribution_sense_as_170419.plx -ref:/home/exchange/bowtie/genomes/dme_miR_mature_miRBase_170119_U_to_T_170123.fasta -map:"$d"  >>sRNA_mapping.log 2>>sRNA_mapping.log
done

echo >>sRNA_mapping.log
echo "----------------------">>sRNA_mapping.log
echo >>sRNA_mapping.log

echo "mapfile analysis finished">>sRNA_mapping.log
echo "mapfile analysis finished"
date>>sRNA_mapping.log
mv -f sRNA_mapping.log logfiles/
rm bowtie/*.log

mkdir bowtie/counts
mv bowtie/*.counts bowtie/counts

mkdir bowtie/sizes
mv bowtie/*.sizes bowtie/sizes

mkdir bowtie/binned
mv bowtie/*.binned bowtie/binned


