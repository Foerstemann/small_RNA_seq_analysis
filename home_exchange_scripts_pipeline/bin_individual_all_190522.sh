#!/bin/bash

set -u
set -o pipefail
if ls *_mapping_loci_190517.map 1>/dev/null 2>&1
then

	for b in *_mapping_loci_190517.map; do
	echo >>sRNA_binning.log
    	echo $b>>sRNA_binning.log
    	echo $b

 
    	echo "mapping reads on bases per fasta line">>sRNA_binning.log
    	echo "mapping reads on bases per fasta line"
        perl /home/exchange/sRNA_pipeline/bin_individual_fasta_lines_190522.plx -ref:/home/exchange/bowtie/genomes/mapping_loci_190517.fa -map:"$b" -bin:1 >>sRNA_binning.log 2>>sRNA_binning.log

	done	

	mkdir binned_loci
	mv *.binned binned_loci/
fi

if ls *_dme_transposons.map  1>/dev/null 2>&1
then

for c in *_dme_transposons.map; do
    echo >>sRNA_binning.log
    echo $c>>sRNA_binning.log
    echo $c

 
    echo "mapping reads on bases per fasta line">>sRNA_binning.log
    echo "mapping reads on bases per fasta line"
        perl /home/exchange/sRNA_pipeline/bin_individual_fasta_lines_190522.plx -ref:/home/exchange/bowtie/genomes/D_mel_transposon_sequence_cons_set_v941.fasta -map:"$c" -bin:1 >>sRNA_binning.log 2>>sRNA_binning.log

done

mkdir binned_TE
mv *.binned binned_TE/
fi
