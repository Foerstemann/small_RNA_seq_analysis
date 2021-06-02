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
my $min_length = 19;
my $max_length = 29;
my $error_flag = 0;
my $in_file_original;
my $in_file;
while ($counter <=40) {
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

	if ($_ =~ /-min/) {
		my @split_entry = split (/:/, $_); 
		$min_length = $split_entry[1];
		}
	if ($_ =~ /-max/) {
		my @split_entry = split (/:/, $_); 
		$max_length = $split_entry[1];
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
while (-e $in_file . "_" . $min_length . "_" . $max_length . ".fastq") {
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


open OUTFILE, "> $in_file\_$min_length\_$max_length\.fastq";


my $start_time = localtime (time);

if ($fasta_flag == 1) {
	clip_fasta();
	}
elsif ($fastq_flag == 1) {
	clip_fastq();
	}

else {
	print "Couldn't find appropriate routine!\n";
	$error_flag = 1;
	}


# process INFILE
sub clip_fasta {

	while (<INFILE>) 
	{
		my $firstline = $_;	
		my $sequence_line = <INFILE>;
		if ($firstline eq "\n") {last};
		chomp $sequence_line;
		++ $lines_before;
		++ $length_before[length($sequence_line)]; 
		if (length($sequence_line) >= $min_length) 
			{
			if (length($sequence_line) <= $max_length)
			{
				print OUTFILE $firstline;
				print OUTFILE $sequence_line, "\n";
			++ $length_after[length($sequence_line)];
			++ $lines_after;
			}
			}
	}
}	


sub clip_fastq {

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
		if (length($sequence_line) >= $min_length) 
			{
			if (length($sequence_line) <= $max_length)
			{
				print OUTFILE $firstline;
				print OUTFILE $sequence_line, "\n";
				print OUTFILE $thirdline;
				print OUTFILE $fourthline;
			++ $length_after[length($sequence_line)];
			++ $lines_after;
			}
			}
	}
}	

$counter = 1;

my $end_time = localtime;
print "processed file: ", $in_file_original, "\n";
print "start time: ", $start_time, "\n";
print "end time: ", $end_time, "\n"; 
print "minimum length: \t", $min_length, "\n";
print "maximum length: \t", $max_length, "\n";
print "reads before clipping: \t", $lines_before, "\n";
print "reads after clipping: \t", $lines_after, "\n";

print "\t\t\tbefore\tafter\n";

while ($counter <= 40) 
{
	print "length $counter nt\t\t", $length_before[$counter], "\t", $length_after[$counter], "\n";

	$counter++;
	}



 







