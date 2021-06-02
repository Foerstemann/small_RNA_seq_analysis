#!/usr/bin/perl
# process_3prime.plx
# trims the 3'-adapter sequences from deep sequencing reads of cloned small RNAs
# this version must receive parameters on command line, cannot be run interactively
# defaults: fastq format, adapter ctgtaggca (miRNA linker 1 from IDT)
# author: Klaus Foerstemann, April 6, 2017
use warnings;
use strict;

# initialize variables
my $path;
my $end_match = 0;
my $reads_retained =0;
my $no_match = 0;
my $fully_removed = 0;
my $counter = 0;
my $multimatch = 0;
my $matchnumber = 0;

$counter = 1;
my $fasta_flag = 0;

my $fastq_flag = 0;
my $error_flag = 0;
my $in_file_original ;
my $adapter_seq = "ctgtaggca";

# print newline to separate output better from command promt
print "\n\n";
# process command-line arguments
foreach (@ARGV) {
	if ($_ eq "-fa") {
		$fasta_flag = 1;
		$fastq_flag = 0;
		}	
	if ($_ eq "-fq") {
		$fasta_flag = 0;
		$fastq_flag = 1;
		}
	if ($_ =~ /-adapter/) {
		my @split_entry = split (/:/, $_); 
		$adapter_seq = $split_entry[1];
		if ($adapter_seq =~ /[^ACGTacgt]/) {
			print "Sorry, adapter sequence contains bad character (allowed are: ACGTacgt).\n";
			undef $adapter_seq;
			}
		}
	if ($_ =~ /-file/) {
		my @split_entry = split (/:/, $_);
		$in_file_original = $split_entry[1];
		if (!-e $in_file_original) {
			print "Sorry, input file does not exist. \n";
			undef $in_file_original;
			}
		}
	if ($_ =~ /-path/) {
		my @split_entry = split (/:/, $_);
		$path = $split_entry[1];
		
		}
	} 



my @in_file_split = split (/\./, $in_file_original);
my $in_file = $in_file_split[0];


open INFILE, $in_file_original or die "$in_file_original not found!\n", $!, "\n";
while (-e $in_file . "_trim3_" . $adapter_seq . ".txt") {
	++ $counter;
	if ($counter >2) {
		my @in_file = split (/\_/, $in_file);
		pop @in_file;
		$in_file = join "_", @in_file;
		}
	$in_file = $in_file . "_" . $counter;
	next
	}
$counter = 1;
open OUTFILE, ">./$path/"."$in_file\_trim3.fastq";

my $start_time = localtime (time);

if ($fasta_flag == 1) {
	trim_fasta();
	}
elsif ($fastq_flag == 1) {
	trim_fastq();
	}

else {
	print "Couldn't find appropriate routine!\n";
	$error_flag = 1;
	}

	


# process INFILE fastq
sub trim_fastq {

	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		chomp $sequence_line;
		my $thirdline = <INFILE>;
		my $fourthline = <INFILE>;
		if ($firstline eq "\n") {last};	
		
		# replace match to adapter only upon the last occurrence
		if ($sequence_line =~ /$adapter_seq/i)
			{
			my $line_buffer = $sequence_line;
			my $line_buffer_2 = $fourthline;	
			$matchnumber = 0;
			my @split_line = split (/$adapter_seq/i, $sequence_line);
			$matchnumber = @split_line;
			if ($matchnumber >> 2) {$multimatch = $multimatch + 1};
			pop @split_line;
			$sequence_line = join ($adapter_seq, @split_line);
			$fourthline = substr ($fourthline,0,length($sequence_line));
			$end_match = $end_match + 1;
			if ($sequence_line =~/N/i) 
				{
				
				++ $fully_removed;
				}
			elsif (length($sequence_line) ==0) 
				{
				
				++ $fully_removed;
				}
			else {
				print OUTFILE $firstline;
				print OUTFILE $sequence_line , "\n";
				print OUTFILE $thirdline;
				print OUTFILE $fourthline, "\n";
				++ $reads_retained;
				}	
			}	
		else	
		{
					
			$no_match = $no_match + 1;
		}
	}
}
# process INFILE fasta
sub trim_fasta {
	
	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		chomp $sequence_line;
		if ($firstline eq "\n") {last};
		
		# replace match to adapter only upon the last occurrence
		if ($sequence_line =~ /$adapter_seq/i)
			{
			my $line_buffer = $sequence_line;
			$matchnumber = 0;
			my @split_line = split (/$adapter_seq/i, $sequence_line);
			$matchnumber = @split_line;
			if ($matchnumber >> 2) {$multimatch = $multimatch + 1};
			pop @split_line;
			$sequence_line = join ($adapter_seq, @split_line);
			
			$end_match = $end_match + 1;
			if ($sequence_line =~/N/i) 
				{
				
				++ $fully_removed;
				}
			elsif (length($sequence_line) ==0 ) 
				{
				
				++ $fully_removed;
				}
			else {
				print OUTFILE $firstline;
				print OUTFILE $sequence_line , "\n";
				++ $reads_retained;
				}	
			}	
		else	
		{
					
			$no_match = $no_match + 1;
		}
	}
}



my $end_time = localtime (time);

	# screen copy;
	print "processed file: $in_file\n", "trimmed sequence tag (3'): $adapter_seq\n\n", "sequence type processed:\n", "fasta ", $fasta_flag, "\n", "fastq ", $fastq_flag, "\n\n";
	print "lines with matches: \t", $end_match, "\n";
	print "lines with multiple matches (only last one removed): \t", $multimatch, "\n";
 	print "lines fully removed or containing N: \t", $fully_removed, "\n";
	print "lines without matches: \t", $no_match, "\n";
	print "reads retained after clipping: \t", $reads_retained, "\n\n";
	print "start time: \t", $start_time, "\n";
	print "end time: \t", $end_time, "\n\n";
	print "\t\t\tbefore\tafter\n";
	
	print "\n"; #separate output better from next command prompt




 







