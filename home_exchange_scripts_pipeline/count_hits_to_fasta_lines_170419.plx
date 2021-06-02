#!/usr/bin/perl
# count_htis_to_fasta_lines.plx
# The script counts read matches to each of the sequences in the FASTA reference file. 
# Script needs the FASTA reference file and a BOWTIE map-file from a run using the same FASTA reference for mapping.
# Author: Klaus Foerstemann, Date: Jul. 17, 2009


use warnings;
use strict;

# initialize variables

my $total_count = 0;
my $sequence_counter = 0;
my @miR_name;
my $z = 0;
my %quant_hash_pos = ();
my %quant_hash_neg = ();
my $lib_file;
my $out_file;
my $bowtie_file;
my $screen_copy_flag = 0;
my $counter = 1;


# print newline to separate output better from command promt
print "\n\n";

# process command-line arguments
foreach (@ARGV) {
		if ($_ =~ /-ref:/) {
		my @split_entry = split (/:/, $_); 
		$lib_file = $split_entry[1];
		if (!-e $lib_file) {
			print "\n Sorry, FASTA reference file not found.\n";
			undef $lib_file;
			}
		}
	if ($_ =~ /-map:/) {
		my @split_entry = split (/:/, $_);
		$bowtie_file = $split_entry[1];
		if (!-e $bowtie_file) {
			print "Sorry, BOWTIE map file does not exist. \n";
			undef $bowtie_file;
			}
		}
	if ($_ =~ /-s/i) {
		$screen_copy_flag = 1;
		}
	} 



# get filename and adapter sequence if not given on command line
while (!defined $lib_file) {
	print "\nPlease enter the filename for the FASTA reference used for mapping:";
	$lib_file = <STDIN>;
	chomp $lib_file;
	if (!-e $lib_file) {
		print "\n Sorry, FASTA reference file not found.\n";
		undef $lib_file;
		next;
		}
	}
	
while (!defined $bowtie_file) {
	print "\nPlease enter the filename for the BOWTIE map file:";
	$bowtie_file = <STDIN>;
	chomp $bowtie_file;
	if (!-e $lib_file) {
		print "\n Sorry, BOWTIE map file not found.\n";
		undef $bowtie_file;
		next;
		}
	}
	

# open filehandles
open LIBFILE, $lib_file;
my @split_bowtie = split /\./, $bowtie_file;
$out_file = $split_bowtie[0];
chomp $out_file;
while (-e $out_file . ".counts") {
	++ $counter;
	if ($counter >2) {
		my @out_file_split = split (/\_/, $out_file);
		pop @out_file_split;
		$out_file = join "_", @out_file_split;
		print $out_file;
		}
	$out_file = $out_file . "_" . $counter;
	next;
	}

open OUTFILE, "> $out_file\.counts";
open LOGFILE, "> $out_file\_counts\.log";

my $start_time = localtime (time);

# get all miRNAs, make array for names to function as keys in hash

while (<LIBFILE>) 
{
	if ($_ =~ />/) 
	{
	my $nameline = $_;	
	my $sequence_line = <LIBFILE>;
	chomp $nameline;
	my @splitname = split />/, $nameline;
	my @splitname2 = split / /, $splitname[1];
	chomp $splitname2[0];
	$miR_name[$sequence_counter] = $splitname2[0];
	$quant_hash_pos{$splitname2[0]} = 0;	
	$quant_hash_neg{$splitname2[0]} = 0;
	++ $sequence_counter;	
	}	
}
close LIBFILE;

# go through bowtie mapfile, get names for hits and increment corresponding hash value (separate for + strand and - strand hits)

open BOWFILE, $bowtie_file or die "$bowtie_file not found!\n", $!, "\n";
while (<BOWFILE>)
	{
	my $mapline = $_;
	my @splitline = split /\t/, $mapline;
	chomp $splitline[2];
	if ($splitline[1] =~ /\+/)
		{
		++ $quant_hash_pos{$splitline[2]};	       
		}
	if ($splitline[1] =~ /\-/)
		{
		++ $quant_hash_neg{$splitline[2]};
		}
	}
close BOWFILE;

my $end_time = localtime (time);

# generate output (screen copy and file)


print OUTFILE "seq. name \t", "pos. strand hits \t", "neg. strand hits\t", "total hits\n";
$z = 0;
while  ($z < $sequence_counter)
{ 
	my $combined_hits = $quant_hash_pos{$miR_name[$z]} + $quant_hash_neg{$miR_name[$z]};
	print OUTFILE $miR_name[$z], "\t", $quant_hash_pos{$miR_name[$z]}, "\t", $quant_hash_neg{$miR_name[$z]}, "\t", $combined_hits, "\n";
	$total_count = $total_count + $quant_hash_pos{$miR_name[$z]} + $quant_hash_neg{$miR_name[$z]};
	++ $z;
}

# screen copy of logfile

print "BOWTIE map file: \t", $bowtie_file, "\n";
print "FASTA reference file: \t", $lib_file, "\n";
print "File with hit counts: \t", $out_file . ".counts\n";
print "number of sequences in FASTA reference file: \t", $sequence_counter, "\n";
print "total numer of matching reads in BOWTIE map file: \t", $total_count, "\n";
print "\nstart time: \t", $start_time, "\n";
print "end time: \t", $end_time, "\n";

# logfile generation

print LOGFILE "BOWTIE map file: \t", $bowtie_file, "\n";
print LOGFILE "FASTA reference file: \t", $lib_file, "\n";
print LOGFILE "File with hit counts: \t", $out_file . ".counts\n";
print LOGFILE "number of sequences in FASTA reference file: \t", $sequence_counter, "\n";
print LOGFILE "total numer of matching reads in BOWTIE map file: \t", $total_count, "\n";
print LOGFILE "\nstart time: \t", $start_time, "\n";
print LOGFILE "end time: \t", $end_time, "\n";

# screen copy of counts (only if selected wiht -s option)
if ($screen_copy_flag ==1) {
	print "seq. name \t", "pos. strand hits \t", "neg. strand hits\t", "total hits\n";
	$z = 0;
	$total_count = 0;
	while  ($z < $sequence_counter)
	{ 
		my $combined_hits = $quant_hash_pos{$miR_name[$z]} + $quant_hash_neg{$miR_name[$z]};
		print $miR_name[$z], "\t", $quant_hash_pos{$miR_name[$z]}, "\t", $quant_hash_neg{$miR_name[$z]}, "\t", $combined_hits, "\n";
		$total_count = $total_count + $quant_hash_pos{$miR_name[$z]} + $quant_hash_neg{$miR_name[$z]};
		++ $z;
	}
	
}

close BOWFILE;
close OUTFILE;
close LOGFILE;

