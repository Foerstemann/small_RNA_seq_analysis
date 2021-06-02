#!/usr/bin/perl
# length_select.plx
# selects deep sequencing reads accoring to length criteria
# this version must receive parameters on command line, cannot be run interactively
# defaults: fastq format, length 19-29 nt
# author: Klaus Foerstemann, April 6, 2017

use warnings;
use strict;

# initialize variables
my @length_before = 0;
my @length_after = 0;
my $counter = 1;
my $lines_before = 0;
my $lines_after = 0;
my $fasta_flag = 0;
my $fastq_flag = 1;
my $clip_length = 19;
my $clipped_end = 5;
my $error_flag = 0;
my $in_file_original;
my $in_file;
while ($counter <=80) {
	$length_before[$counter] = 0;
	$length_after[$counter] = 0;
	$counter++;
	}
$counter =1;
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

	if ($_ =~ /-length/) {
		my @split_entry = split (/:/, $_); 
		$clip_length = $split_entry[1];
		}

	if ($_ =~ /-end/) {
		my @split_entry = split (/:/, $_); 
		$clipped_end = $split_entry[1];
		}

		
	if ($_ =~ /-file/) {
		my @split_entry = split (/:/, $_);
		$in_file_original = $split_entry[1];
		if (!-e $in_file_original) {
			print "Sorry, input file does not exist. \n";
			undef $in_file_original;
			}
		}
	} 



	



open INFILE, $in_file_original or die "$in_file_original not found!\n", $!, "\n";
my @split_file = split (/\./, $in_file_original);
$in_file = $split_file[0];
while (-e $in_file . "_" .$clipped_end."'_" . $clip_length . "clip" . ".fastq") {
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


open OUTFILE, "> $in_file\_$clip_length\_clip\.fastq";


my $start_time = localtime (time);

if ($fasta_flag == 1) {
	if ($clipped_end == 3) {
		clip_fasta_3();
		}
	else clip_fasta_5();
	
	}
elsif ($fastq_flag == 1) {
	if ($clipped_end == 3) {
		clip_fastq_3();
		}
	else clip_fastq_5();	
	}

else {
	print "Couldn't find appropriate routine!\n";
	$error_flag = 1;
	}


# process INFILE
sub clip_fasta_5 {

	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		if ($firstline eq "\n") {last};
		chomp $sequence_line;
		++ $lines_before;
		++ $length_before[length($sequence_line)]; 

				print OUTFILE $firstline;
				print OUTFILE substr($sequence_line, 0, $clip_length), "\n";
			++ $length_after[length(substr($sequence_line, 0, $clip_length))];
			++ $lines_after;
		}
			
	}

sub clip_fasta_3 {

	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		if ($firstline eq "\n") {last};
		chomp $sequence_line;
		++ $lines_before;
		++ $length_before[length($sequence_line)]; 

				print OUTFILE $firstline;
				print OUTFILE substr($sequence_line, -$clip_length, $clip_length), "\n";
			++ $length_after[length(substr($sequence_line, -$clip_length, $clip_length))];
			++ $lines_after;
		}
			
	}	


sub clip_fastq_5 {

	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		chomp $sequence_line;
		my $thirdline = <INFILE>;
		my $fourthline = <INFILE>;
		if ($firstline eq "\n") {last};
		++ $lines_before;
		++ $length_before[length($sequence_line)]; 
				print OUTFILE $firstline;
				print OUTFILE substr($sequence_line, 0, $clip_length), "\n";
				print OUTFILE $thirdline;
				print OUTFILE substr($fourthline, 0, $clip_length), "\n";
			++ $length_after[length(substr($sequence_line, 0, $clip_length))];
			++ $lines_after;
			
			
	}
}

sub clip_fastq_3 {

	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		chomp $sequence_line;
		my $thirdline = <INFILE>;
		my $fourthline = <INFILE>;
		if ($firstline eq "\n") {last};
		++ $lines_before;
		++ $length_before[length($sequence_line)]; 
				print OUTFILE $firstline;
				print OUTFILE substr($sequence_line, -$clip_length, $clip_length), "\n";
				print OUTFILE $thirdline;
				print OUTFILE substr($fourthline, -$clip_length, $clip_length), "\n";
			++ $length_after[length(substr($sequence_line, -$clip_length, $clip_length))];
			++ $lines_after;
			
			
	}
}	

$counter = 1;

my $end_time = localtime;
print "processed file: ", $in_file_original, "\n";
print "start time: ", $start_time, "\n";
print "end time: ", $end_time, "\n"; 
print "clip length: \t", $clip_length, " from the ", $clipped_end, "'end"\n";
print "reads before clipping: \t", $lines_before, "\n";
print "reads after clipping: \t", $lines_after, "\n";

print "\t\t\tbefore\tafter\n";

while ($counter <= 80) 
{
	print "length $counter nt\t\t", $length_before[$counter], "\t", $length_after[$counter], "\n";

	$counter++;
	}



 







