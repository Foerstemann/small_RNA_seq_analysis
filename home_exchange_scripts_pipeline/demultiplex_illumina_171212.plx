#!/usr/bin/perl
# demultiplex_illumina.plx
#
# command line arguments: -index: [index fastq file], -reads: [reads fastq file] -l: [integer defining the length of the index to be analyzed]
# The index sequences for demultiplexing are hard-coded at the beginning of the script. 
# The script populates a hash with the read names, then  assigns a value according to the index read.
# The reads-containing fastq-file is parsed in parallel and split according to the indexes using the pre-populated hash as a reference.
# If the two files are "in phase", then there will be rather little memory use and the script runs fast.
# It can deal, however, with out-of-phase index and read files as well as corrupted fastq files containing additional lines or missing lines. 
# In these cases, however, one or two reads surrounding the problematic lines will be dropped. 
# Furthermore, the script may require a lot of memory because the hash is largely populated before the reads are then assigned in a second pass. 
# Author: Klaus Foerstemann, Date: December 12, 2017


use warnings;
# use strict;

# initialize variables
my @index_name;
my @index;
my $z = 0;
my %read_to_index = ();
my $index_file;
my $reads_file;
my @outfile;
my $index_read_counter = 0;
my $sequence_read_counter = 0;
my @demultiplexed_reads_counter;
my $firstline_reads;
my $sequence_line_reads;
my $thirdline_reads;
my $fourthline_reads;
my $fh_reads;
my $second_pass_reads = 0;
my $raw_counts=0;

# define index read sequences and corresponding names
my $index_number = 9;
my $index_length = 6;
$index[0] = "NNNNNNNN";
$index_name[0]= "unassigned";
$index[1] = "CGTGATAT";
$index_name[1] = "i1";
$index[2] = "ACATCGAT";
$index_name[2] = "i2";
$index[3] = "GCCTAAAT";
$index_name[3] = "i3";
$index[4] = "TGGTCAAT";
$index_name[4] = "i4";
$index[5] = "ATTCTGAT";
$index_name[5] = "i5";
$index[6] = "TACAGCAT";
$index_name[6] = "i6";
$index[7] = "GTGCGAAT";
$index_name[7] = "i7";
$index[8] = "CCACTGAT";
$index_name[8] = "i8";
$index[9] = "GATACGAT";
$index_name[9] = "i9";

for ($z=0; $z<=$index_number; $z++){
	$demultiplexed_reads_counter[$z]=0;
	}

my $screen_copy_flag = 0;

# print newline to separate output better from command promt
print "\n\n";

# process command-line arguments
foreach (@ARGV) {
		if ($_ =~ /-index:/) {
		my @split_entry = split (/:/, $_); 
		$index_file = $split_entry[1];
		if (!-e $index_file) {
			print "\n Sorry, index reads file not found.\n";
			undef $index_file;
			}
		}
	if ($_ =~ /-reads:/) {
		my @split_entry = split (/:/, $_);
		$reads_file = $split_entry[1];
		if (!-e $reads_file) {
			print "Sorry, sequence reads file not found.\n";
			undef $reads_file;
			}
		}
	if ($_ =~ /-l:/i) {
		my @split_entry = split (/:/, $_);
		$index_length = $split_entry[1];
		}
	} 



# get filename and adapter sequence if not given on command line
while (!defined $index_file) {
	print "\nPlease enter the filename for the index file:";
	$index_file = <STDIN>;
	chomp $index_file;
	if (!-e $index_file) {
		print "\n Sorry, FASTA reference file not found.\n";
		undef $index_file;
		next;
		}
	}
	
while (!defined $reads_file) {
	print "\nPlease enter the filename for the reads file:";
	$reads_file = <STDIN>;
	chomp $reads_file;
	if (!-e $reads_file) {
		print "\n Sorry, reads file not found.\n";
		undef $reads_file;
		next;
		}
	}
	
	
	
	
my $start_time = localtime (time);

# open filehandles
open INDEXFILE, $index_file or die "$index_file not found!\n", $!, "\n";
open READSFILE, $reads_file or die "$reads_file not found!\n", $!, "\n";


# open one outfile for each index and outfile[0] for unassigned reads 
my @split_readname = split /\./, $reads_file;
for ($z=0; $z<=$index_number; $z++) {
	open ($index_name[$z], ">", $split_readname[0]."_".$index_name[$z].".fastq");
	}

#open temporary file for reads that are not "in phase"

open SECOND_PASS, ">temporary_read_storage.fastq"; 	 


# populate the hash with sequence names and assign value according to index sequence, parse reads file in parallel

while (<INDEXFILE>) 
{
	my $match_flag = 0;
	my $firstline =$_;
	if ($firstline=~/^@/) {
		my @splitline= split ' ', $firstline;
		$firstline=$splitline[0];
		my $sequence_line = <INDEXFILE>;
		my $thirdline = <INDEXFILE>;
		my $fourthline = <INDEXFILE>;
		for ($z=1; $z<=$index_number; $z++) {
		
			if (substr($sequence_line,0,$index_length) =~ substr($index[$z],0,$index_length)) 
				{$read_to_index{$firstline}=$z;
				$match_flag = 1;
				last;
				}
			}
			if ($match_flag == 0) {
				$read_to_index{$firstline}=0;
				}
		$index_read_counter ++;
		# check if the readsfile still has some reads, if not go next in while loop for index
		if (eof(READSFILE)) {next};

		$firstline_reads = <READSFILE>;
	        while ($firstline_reads !~/^@/) {
			$firstline_reads = <READSFILE>;
		}				
			
		my @splitline_reads= split ' ', $firstline_reads;
		$firstline_reads=$splitline_reads[0];
		$sequence_line_reads = <READSFILE>;
		$thirdline_reads = <READSFILE>;
		$fourthline_reads = <READSFILE>;
		if (exists($read_to_index{$firstline_reads})) {		
			$fh_reads=$index_name[$read_to_index{$firstline_reads}];
			print $fh_reads $firstline_reads,"\n", $sequence_line_reads, $thirdline_reads, $fourthline_reads;
			
			$demultiplexed_reads_counter[$read_to_index{$firstline_reads}] ++;
			delete $read_to_index{$firstline_reads}; # reomve hash value to recover memory
			$raw_counts++;
			}
		else	{
		print SECOND_PASS $firstline_reads,"\n", $sequence_line_reads, $thirdline_reads, $fourthline_reads;
		}
	}
}
# now check in readsfile was completely parsed, copy rest into SECOND_PASS if not
if (!eof(READSFILE)) {
		$firstline_reads = <READSFILE>;
	        while ($firstline_reads !~/^@/) {
			$firstline_reads = <READSFILE>;
		}				
			
		my @splitline_reads= split ' ', $firstline_reads;
		$firstline_reads=$splitline_reads[0];
		$sequence_line_reads = <READSFILE>;
		$thirdline_reads = <READSFILE>;
		$fourthline_reads = <READSFILE>;
		if (exists($read_to_index{$firstline_reads})) {		
			$fh_reads=$index_name[$read_to_index{$firstline_reads}];
			print $fh_reads $firstline_reads,"\n", $sequence_line_reads, $thirdline_reads, $fourthline_reads;
			
			$demultiplexed_reads_counter[$read_to_index{$firstline_reads}] ++;
			delete $read_to_index{$firstline_reads}; # reomve hash value to recover memory
			$raw_counts++;
			}
		else	{
		print SECOND_PASS $firstline_reads,"\n", $sequence_line_reads, $thirdline_reads, $fourthline_reads;
		}
	}
close INDEXFILE;
close READSFILE;
close SECOND_PASS;
my $hash_time=localtime;
print $hash_time, "\n";	
print "index file completely processed. \n\n";
print "reads after first pass:\t", $raw_counts, "\n";
	
# go through temporary seuqence file one by one and assign to index according to the read name
open SECOND_PASS, "temporary_read_storage.fastq"; 
print unassigned "second_pass\n";
while (<SECOND_PASS>)
	{
	my $firstline_second = $_;
		 
		while ($firstline_second !~/^@/){
		if (!eof (SECOND_PASS)) {$firstline_second = <SECOND_PASS>};
		}
	my @splitline= split ' ', $firstline_second;
	$firstline_second=$splitline[0];
	my $sequence_line_second = <SECOND_PASS>;
	my $thirdline_second = <SECOND_PASS>;
	my $fourthline_second = <SECOND_PASS>;
	if (exists($read_to_index{$firstline_second})) {		
		my $fh=$index_name[$read_to_index{$firstline_second}];
		print $fh $firstline_second,"\n", $sequence_line_second, $thirdline_second, $fourthline_second;
		$demultiplexed_reads_counter[$read_to_index{$firstline_second}] ++;
		$raw_counts++;
		$reads_second_pass ++;
		}
	
	else {
		my $fh = $index_name[0];
		print $fh $firstline_second, "\n", $sequence_line_second, $thirdline_second, $fourthline_second;
		$demultiplexed_reads_counter[0] ++;
		$raw_counts++;
		$reads_second_pass ++;
		}
	}

close SECOND_PASS;	
	

# close demultiplexed outfiles
for ($z=0; $z<=$index_number; $z++) {
	close ($index_name[$z]);
	}

# sum up all demultiplexed reads
$sequence_read_counter += $_ foreach(@demultiplexed_reads_counter);

my $end_time = localtime (time);

# screen copy of information
print "\nstart time: \t", $start_time, "\n";
print "end time: \t", $end_time, "\n";
print "index length analzyed for sorting:\t", $index_length, "\n";
print "index file: \t\t", $index_file, "\n";
print "number of indexes: \t", $index_read_counter, "\n";
print "finish time hash construction: \t", $hash_time, "\n";
print "reads file: \t\t", $reads_file, "\n";
print "number of total reads: \t", $raw_counts, "\n";
print "sorted reads: \t", $sequence_read_counter, "\n";
print "reads sorted in second pass: \t", $reads_second_pass, "\n";

for (my $z=0; $z<=$index_number; $z++) {
	print "index name:\t", $index_name[$z], "\tindex sequence:\t",$index[$z],"\t sorted reads:\t", $demultiplexed_reads_counter[$z], "\n";
	}
	


# logfile generation
open LOGFILE, ">", $split_readname[0].".log";
print LOGFILE "\nstart time: \t", $start_time, "\n";
print LOGFILE "end time: \t", $end_time, "\n";
print LOGFILE "index length analzyed for sorting:\t", $index_length, "\n";
print LOGFILE "index file: \t\t", $index_file, "\n";
print LOGFILE "number of indexes: \t", $index_read_counter, "\n";
print LOGFILE "finish time hash construction: \t", $hash_time, "\n";
print LOGFILE "reads file: \t\t", $reads_file, "\n";
print LOGFILE "number of total reads: \t", $raw_counts, "\n";
print LOGFILE "sorted reads: \t", $sequence_read_counter, "\n";
print LOGFILE "reads sorted in second pass: \t", $reads_second_pass, "\n";


for ($z=0; $z<=$index_number; $z++) {
	print LOGFILE "index name:\t", $index_name[$z], "\tindex sequence:\t",$index[$z],"\t sorted reads:\t", $demultiplexed_reads_counter[$z], "\n";
	}
close LOGFILE;	

# delete temporary file
unlink "temporary_read_storage.fastq";
