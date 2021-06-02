#!/usr/bin/perl
# quantify_size_distribution.plx
# from an sequence file (FASTA) and the Bowtie Output, matching read lengths are calculated
# author: Klaus Foerstemann, Jul. 21, 2009


use warnings;
use strict;

# initialize variables
my $lib_file;
my $bowtie_file;
my $out_file;
my $counter;
my $total_count = 0;
my @global_length_count_s = 0;
my @global_length_count_as = 0;
my $sequence_counter = 0;
my @target_name;
my @target_length;
my @splitline;
my $z = 0;
my $i;
$i = 14;
while ($i <=32)	
	{
	$global_length_count_s[$i] = 0;
	$global_length_count_as[$i] = 0;
	++ $i;
	} 
my %nt14_hash_s = ();
my %nt15_hash_s = ();
my %nt16_hash_s =();	
my %nt17_hash_s = ();
my %nt18_hash_s = ();
my %nt19_hash_s = ();
my %nt20_hash_s = ();
my %nt21_hash_s = ();
my %nt22_hash_s = ();
my %nt23_hash_s = ();
my %nt24_hash_s = ();
my %nt25_hash_s = ();
my %nt26_hash_s = ();
my %nt27_hash_s = ();
my %nt28_hash_s = ();
my %nt29_hash_s = ();
my %nt30_hash_s = ();
my %nt31_hash_s = ();
my %nt32_hash_s = ();


my %nt14_hash_as = ();
my %nt15_hash_as = ();
my %nt16_hash_as =();	
my %nt17_hash_as = ();
my %nt18_hash_as = ();
my %nt19_hash_as = ();
my %nt20_hash_as = ();
my %nt21_hash_as = ();
my %nt22_hash_as = ();
my %nt23_hash_as = ();
my %nt24_hash_as = ();
my %nt25_hash_as = ();
my %nt26_hash_as = ();
my %nt27_hash_as = ();
my %nt28_hash_as = ();
my %nt29_hash_as = ();
my %nt30_hash_as = ();
my %nt31_hash_as = ();
my %nt32_hash_as = ();
my @sequence_line;
my $screen_copy_flag = 0;

# process command line arguments
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
open LIBFILE, $lib_file or die $!;
my @split_bowtie = split /\./, $bowtie_file;
$out_file = $split_bowtie[0];
chomp $out_file;
while (-e $out_file . ".sizes") {
	++ $counter;
	if ($counter >2) {
		my @out_file_split = split (/\_/, $out_file);
		pop @out_file_split;
		$out_file = join "_", @out_file_split;
		}
	$out_file = $out_file . "_" . $counter;
	next;
	}

open OUTFILE, "> $out_file\.sizes";
open BOWFILE, $bowtie_file or die $!;
	


my $start_time = localtime (time);

# get all miRNAs, make array for names to function as keys in hash

while (<LIBFILE>) 
{	
	if ($_ =~ />/) 
	{
	my $nameline = $_;	
	$sequence_line[$sequence_counter] = <LIBFILE>;
	chomp $nameline;
	chomp $sequence_line[$sequence_counter];
	my @splitname = split />/, $nameline;
	my @splitname2 = split / /, $splitname[1];
	chomp $splitname2[0];
	$target_name[$sequence_counter] = $splitname2[0];
	$target_length[$sequence_counter] = length ($sequence_line[$sequence_counter]);

	my $name = $target_name[$sequence_counter];
	$nt14_hash_s{$name} = 0;
	$nt15_hash_s{$name} = 0;
	$nt16_hash_s{$name} = 0;
	$nt17_hash_s{$name} = 0;
	$nt18_hash_s{$name} = 0;
	$nt19_hash_s{$name} = 0;
	$nt20_hash_s{$name} = 0;
	$nt21_hash_s{$name} = 0;
	$nt22_hash_s{$name} = 0;
	$nt23_hash_s{$name} = 0;
	$nt24_hash_s{$name} = 0;
	$nt25_hash_s{$name} = 0;
	$nt26_hash_s{$name} = 0;
	$nt27_hash_s{$name} = 0;
	$nt28_hash_s{$name} = 0;	
	$nt29_hash_s{$name} = 0;
	$nt30_hash_s{$name} = 0;
	$nt31_hash_s{$name} = 0;
	$nt32_hash_s{$name} = 0;

	$nt14_hash_as{$name} = 0;
	$nt15_hash_as{$name} = 0;
	$nt16_hash_as{$name} = 0;
	$nt17_hash_as{$name} = 0;
	$nt18_hash_as{$name} = 0;
	$nt19_hash_as{$name} = 0;
	$nt20_hash_as{$name} = 0;
	$nt21_hash_as{$name} = 0;
	$nt22_hash_as{$name} = 0;
	$nt23_hash_as{$name} = 0;
	$nt24_hash_as{$name} = 0;
	$nt25_hash_as{$name} = 0;
	$nt26_hash_as{$name} = 0;
	$nt27_hash_as{$name} = 0;
	$nt28_hash_as{$name} = 0;	
	$nt29_hash_as{$name} = 0;
	$nt30_hash_as{$name} = 0;
	$nt31_hash_as{$name} = 0;
	$nt32_hash_as{$name} = 0;	
	++ $sequence_counter;	
	}	
}
close LIBFILE;

# go through bowtie mapfile, get names for hits and increment corresponding hash value


print "going through bowtie mapfile...";
while (<BOWFILE>)
	{
	++ $total_count;
	my $mapline = $_;
	@splitline = split /\t/, $mapline;
	chomp $splitline[1];
if ($splitline[1] eq "+")
	{	
	chomp $splitline[2];
	chomp $splitline[4];	
	my $readlength = length($splitline[4]);
	++ $global_length_count_s[$readlength];
	my $name = $splitline[2];
	if ($readlength == 14)
			{
			$nt14_hash_s{$name} ++;				
			}
	elsif ($readlength == 15)
			{
			$nt15_hash_s{$name} ++;				
			}	
	elsif ($readlength == 16)
			{
			$nt16_hash_s{$name} ++;				
			}	
	elsif ($readlength == 17)
			{
			$nt17_hash_s{$name} ++;				
			}	
	elsif ($readlength == 18)
			{
			$nt18_hash_s{$name} ++;				
			}			
	elsif ($readlength == 19)
			{
			++ $nt19_hash_s{$name};				
			}			
	elsif ($readlength == 20)
			{
			++ $nt20_hash_s{$name};				
			}			
	elsif ($readlength == 21)
			{
			++ $nt21_hash_s{$name};				
			}			
	elsif ($readlength == 22)
			{
			++ $nt22_hash_s{$name};				
			}		
	elsif ($readlength == 23)
			{
			++ $nt23_hash_s{$name};				
			}			
	elsif ($readlength == 24)
			{
			++ $nt24_hash_s{$name};				
			}
	elsif ($readlength == 25)
			{
			++ $nt25_hash_s{$name};				
			}
	elsif ($readlength == 26)
			{
			++ $nt26_hash_s{$name};				
			}
	elsif ($readlength == 27)
			{
			++ $nt27_hash_s{$name};				
			}
	elsif ($readlength == 28)
			{
			++ $nt28_hash_s{$name};				
			}
	elsif ($readlength == 29)
			{
			++ $nt29_hash_s{$name};				
			}
	elsif ($readlength == 30)
			{
			++ $nt30_hash_s{$name};				
			}
	elsif ($readlength == 31)
			{
			++ $nt31_hash_s{$name};				
			}
	elsif ($readlength == 32)
			{
			++ $nt32_hash_s{$name};				
			}
		}
if ($splitline[1] eq "-")
	{	
	chomp $splitline[2];
	chomp $splitline[4];	
	my $readlength = length($splitline[4]);
	++ $global_length_count_as[$readlength];
	my $name = $splitline[2];
	if ($readlength == 14)
			{
			$nt14_hash_as{$name} ++;				
			}
	elsif ($readlength == 15)
			{
			$nt15_hash_as{$name} ++;				
			}	
	elsif ($readlength == 16)
			{
			$nt16_hash_as{$name} ++;				
			}	
	elsif ($readlength == 17)
			{
			$nt17_hash_as{$name} ++;				
			}	
	elsif ($readlength == 18)
			{
			$nt18_hash_as{$name} ++;				
			}			
	elsif ($readlength == 19)
			{
			++ $nt19_hash_as{$name};				
			}			
	elsif ($readlength == 20)
			{
			++ $nt20_hash_as{$name};				
			}			
	elsif ($readlength == 21)
			{
			++ $nt21_hash_as{$name};				
			}			
	elsif ($readlength == 22)
			{
			++ $nt22_hash_as{$name};				
			}		
	elsif ($readlength == 23)
			{
			++ $nt23_hash_as{$name};				
			}			
	elsif ($readlength == 24)
			{
			++ $nt24_hash_as{$name};				
			}
	elsif ($readlength == 25)
			{
			++ $nt25_hash_as{$name};				
			}
	elsif ($readlength == 26)
			{
			++ $nt26_hash_as{$name};				
			}
	elsif ($readlength == 27)
			{
			++ $nt27_hash_as{$name};				
			}
	elsif ($readlength == 28)
			{
			++ $nt28_hash_as{$name};				
			}
	elsif ($readlength == 29)
			{
			++ $nt29_hash_as{$name};				
			}
	elsif ($readlength == 30)
			{
			++ $nt30_hash_as{$name};				
			}
	elsif ($readlength == 31)
			{
			++ $nt31_hash_as{$name};				
			}
	elsif ($readlength == 32)
			{
			++ $nt32_hash_as{$name};				
			}
		}
	}
close BOWFILE;

my $end_time = localtime (time);

# generate output 


$z = 0;
while  ($z < $sequence_counter)
{ 
			

	print OUTFILE "\n", $target_name[$z], "\n";
	
	print OUTFILE "length \tsense\tantis.\ttotal\n";
	
	print OUTFILE "14 \t", $nt14_hash_s{$target_name[$z]},"\t", $nt14_hash_as{$target_name[$z]}, "\t", ($nt14_hash_s{$target_name[$z]}+$nt14_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "15 \t", $nt15_hash_s{$target_name[$z]},"\t", $nt15_hash_as{$target_name[$z]}, "\t", ($nt15_hash_s{$target_name[$z]}+$nt15_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "16 \t", $nt16_hash_s{$target_name[$z]},"\t", $nt16_hash_as{$target_name[$z]}, "\t", ($nt16_hash_s{$target_name[$z]}+$nt16_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "17 \t", $nt17_hash_s{$target_name[$z]},"\t", $nt17_hash_as{$target_name[$z]}, "\t", ($nt17_hash_s{$target_name[$z]}+$nt17_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "18 \t", $nt18_hash_s{$target_name[$z]},"\t", $nt18_hash_as{$target_name[$z]}, "\t", ($nt18_hash_s{$target_name[$z]}+$nt18_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "19 \t", $nt19_hash_s{$target_name[$z]},"\t", $nt19_hash_as{$target_name[$z]}, "\t", ($nt19_hash_s{$target_name[$z]}+$nt19_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "20 \t", $nt20_hash_s{$target_name[$z]},"\t", $nt20_hash_as{$target_name[$z]}, "\t", ($nt20_hash_s{$target_name[$z]}+$nt20_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "21 \t", $nt21_hash_s{$target_name[$z]},"\t", $nt21_hash_as{$target_name[$z]}, "\t", ($nt21_hash_s{$target_name[$z]}+$nt21_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "22 \t", $nt22_hash_s{$target_name[$z]},"\t", $nt22_hash_as{$target_name[$z]}, "\t", ($nt22_hash_s{$target_name[$z]}+$nt22_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "23 \t", $nt23_hash_s{$target_name[$z]},"\t", $nt23_hash_as{$target_name[$z]}, "\t", ($nt23_hash_s{$target_name[$z]}+$nt23_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "24 \t", $nt24_hash_s{$target_name[$z]},"\t", $nt24_hash_as{$target_name[$z]}, "\t", ($nt24_hash_s{$target_name[$z]}+$nt24_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "25 \t", $nt25_hash_s{$target_name[$z]},"\t", $nt25_hash_as{$target_name[$z]}, "\t", ($nt25_hash_s{$target_name[$z]}+$nt25_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "26 \t", $nt26_hash_s{$target_name[$z]},"\t", $nt26_hash_as{$target_name[$z]}, "\t", ($nt26_hash_s{$target_name[$z]}+$nt26_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "27 \t", $nt27_hash_s{$target_name[$z]},"\t", $nt27_hash_as{$target_name[$z]}, "\t", ($nt27_hash_s{$target_name[$z]}+$nt27_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "28 \t", $nt28_hash_s{$target_name[$z]},"\t", $nt28_hash_as{$target_name[$z]}, "\t", ($nt28_hash_s{$target_name[$z]}+$nt28_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "29 \t", $nt29_hash_s{$target_name[$z]},"\t", $nt29_hash_as{$target_name[$z]}, "\t", ($nt29_hash_s{$target_name[$z]}+$nt29_hash_as{$target_name[$z]}), "\n";	
	print OUTFILE "30 \t", $nt30_hash_s{$target_name[$z]},"\t", $nt30_hash_as{$target_name[$z]}, "\t", ($nt30_hash_s{$target_name[$z]}+$nt30_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "31 \t", $nt31_hash_s{$target_name[$z]},"\t", $nt31_hash_as{$target_name[$z]}, "\t", ($nt31_hash_s{$target_name[$z]}+$nt31_hash_as{$target_name[$z]}), "\n";
	print OUTFILE "32 \t", $nt32_hash_s{$target_name[$z]},"\t", $nt32_hash_as{$target_name[$z]}, "\t", ($nt32_hash_s{$target_name[$z]}+$nt32_hash_as{$target_name[$z]}), "\n";

	++ $z;
}

if ($screen_copy_flag ==1) {	
$z = 0;
while  ($z < $sequence_counter)
{ 
			
	print  "\n", $target_name[$z], "\n";
	
	print  "length \tsense\tantis.\ttotal\n";
	
	print  "14 \t", $nt14_hash_s{$target_name[$z]},"\t", $nt14_hash_as{$target_name[$z]}, "\t", ($nt14_hash_s{$target_name[$z]}+$nt14_hash_as{$target_name[$z]}), "\n";
	print  "15 \t", $nt15_hash_s{$target_name[$z]},"\t", $nt15_hash_as{$target_name[$z]}, "\t", ($nt15_hash_s{$target_name[$z]}+$nt15_hash_as{$target_name[$z]}), "\n";
	print  "16 \t", $nt16_hash_s{$target_name[$z]},"\t", $nt16_hash_as{$target_name[$z]}, "\t", ($nt16_hash_s{$target_name[$z]}+$nt16_hash_as{$target_name[$z]}), "\n";
	print  "17 \t", $nt17_hash_s{$target_name[$z]},"\t", $nt17_hash_as{$target_name[$z]}, "\t", ($nt17_hash_s{$target_name[$z]}+$nt17_hash_as{$target_name[$z]}), "\n";
	print  "18 \t", $nt18_hash_s{$target_name[$z]},"\t", $nt18_hash_as{$target_name[$z]}, "\t", ($nt18_hash_s{$target_name[$z]}+$nt18_hash_as{$target_name[$z]}), "\n";
	print  "19 \t", $nt19_hash_s{$target_name[$z]},"\t", $nt19_hash_as{$target_name[$z]}, "\t", ($nt19_hash_s{$target_name[$z]}+$nt19_hash_as{$target_name[$z]}), "\n";
	print  "20 \t", $nt20_hash_s{$target_name[$z]},"\t", $nt20_hash_as{$target_name[$z]}, "\t", ($nt20_hash_s{$target_name[$z]}+$nt20_hash_as{$target_name[$z]}), "\n";
	print  "21 \t", $nt21_hash_s{$target_name[$z]},"\t", $nt21_hash_as{$target_name[$z]}, "\t", ($nt21_hash_s{$target_name[$z]}+$nt21_hash_as{$target_name[$z]}), "\n";
	print  "22 \t", $nt22_hash_s{$target_name[$z]},"\t", $nt22_hash_as{$target_name[$z]}, "\t", ($nt22_hash_s{$target_name[$z]}+$nt22_hash_as{$target_name[$z]}), "\n";
	print  "23 \t", $nt23_hash_s{$target_name[$z]},"\t", $nt23_hash_as{$target_name[$z]}, "\t", ($nt23_hash_s{$target_name[$z]}+$nt23_hash_as{$target_name[$z]}), "\n";
	print  "24 \t", $nt24_hash_s{$target_name[$z]},"\t", $nt24_hash_as{$target_name[$z]}, "\t", ($nt24_hash_s{$target_name[$z]}+$nt24_hash_as{$target_name[$z]}), "\n";
	print  "25 \t", $nt25_hash_s{$target_name[$z]},"\t", $nt25_hash_as{$target_name[$z]}, "\t", ($nt25_hash_s{$target_name[$z]}+$nt25_hash_as{$target_name[$z]}), "\n";
	print  "26 \t", $nt26_hash_s{$target_name[$z]},"\t", $nt26_hash_as{$target_name[$z]}, "\t", ($nt26_hash_s{$target_name[$z]}+$nt26_hash_as{$target_name[$z]}), "\n";
	print  "27 \t", $nt27_hash_s{$target_name[$z]},"\t", $nt27_hash_as{$target_name[$z]}, "\t", ($nt27_hash_s{$target_name[$z]}+$nt27_hash_as{$target_name[$z]}), "\n";
	print  "28 \t", $nt28_hash_s{$target_name[$z]},"\t", $nt28_hash_as{$target_name[$z]}, "\t", ($nt28_hash_s{$target_name[$z]}+$nt28_hash_as{$target_name[$z]}), "\n";
	print  "29 \t", $nt29_hash_s{$target_name[$z]},"\t", $nt29_hash_as{$target_name[$z]}, "\t", ($nt29_hash_s{$target_name[$z]}+$nt29_hash_as{$target_name[$z]}), "\n";	
	print  "30 \t", $nt30_hash_s{$target_name[$z]},"\t", $nt30_hash_as{$target_name[$z]}, "\t", ($nt30_hash_s{$target_name[$z]}+$nt30_hash_as{$target_name[$z]}), "\n";
	print  "31 \t", $nt31_hash_s{$target_name[$z]},"\t", $nt31_hash_as{$target_name[$z]}, "\t", ($nt31_hash_s{$target_name[$z]}+$nt31_hash_as{$target_name[$z]}), "\n";
	print  "32 \t", $nt32_hash_s{$target_name[$z]},"\t", $nt32_hash_as{$target_name[$z]}, "\t", ($nt32_hash_s{$target_name[$z]}+$nt32_hash_as{$target_name[$z]}), "\n";

	++ $z;
}
}
	



# screen copy of log

print "\nBOWTIE map file: \t", $bowtie_file, "\n";
print "FASTA file of sequences: \t", $lib_file, "\n";
print "number of sequences in FASTA input file: \t", $sequence_counter, "\n";
print "total numer of matching reads in BOWTIE map file: \t", $total_count, "\n";
print "\nglobal length distribution:\n";
print "length \tsense\tantis.\ttotal\n";
$i = 14;
while ($i <=32)	
	{
	print $i, "\t", $global_length_count_s[$i], "\t", $global_length_count_as[$i], "\t", ($global_length_count_s[$i]+ $global_length_count_as[$i]),"\n";
	++ $i;
	}  
print "\nstart time: \t", $start_time, "\n";
print "end time: \t", $end_time, "\n";

# log appended to OUTFILE

print OUTFILE "\nBOWTIE map file: \t", $bowtie_file, "\n";
print OUTFILE "FASTA file of sequences: \t", $lib_file, "\n";
print OUTFILE "number of sequences in FASTA input file: \t", $sequence_counter, "\n";
print OUTFILE "total numer of matching reads in BOWTIE map file: \t", $total_count, "\n";
print OUTFILE "\nglobal length distribution:\n";
print OUTFILE "length\tsense\tantis.\ttotal\n";
$i = 14;
while ($i <=32)	
	{
	print OUTFILE $i, "\t", $global_length_count_s[$i], "\t", $global_length_count_as[$i], "\t", ($global_length_count_s[$i]+ $global_length_count_as[$i]),"\n";
	++ $i;
	} 
print OUTFILE "\nstart time: \t", $start_time, "\n";
print OUTFILE "end time: \t", $end_time, "\n";


close BOWFILE;
close OUTFILE;
close OUTFILE;

