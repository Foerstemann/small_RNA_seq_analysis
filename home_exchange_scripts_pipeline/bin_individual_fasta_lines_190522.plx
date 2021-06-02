#!/usr/bin/perl
# quantify_in_bins.plx


use warnings;
# use strict;

# initialize variables

my $lib_file;
my $bowtie_file;
my $out_file;
my $counter;
my $i;
my $sequence_counter = 0;
my $bin_size;
my @items;
my $bin_number = 0;
my @pos_strand;
my @neg_strand;
my @total_bins;
my @seq_name;
my @seq_length;
my $temp_hashname_pos;
my $temp_hashname_neg;
my @pos_count;
my @neg_count;
my %name_hash;
my $in_line;
my $next_name;
my $current_name;
my @splitname;
my @splitname2;


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
	if ($_ =~ /-bin:/) {
		my @split_entry = split (/:/, $_);
		$bin_size = $split_entry[1];
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
	
while (!defined $bin_size) {
	print "\nPlease enter the interval size for binning of the reads:";
	$bin_size = <STDIN>;
	chomp $bin_size;
	}

# open filehandles
open LIBFILE, $lib_file or die $!;
my @split_bowtie = split /\./, $bowtie_file;
$out_file = $split_bowtie[0];
chomp $out_file;
$out_file = $out_file . "_bin" . $bin_size;
while (-e $out_file . ".binned") {
	++ $counter;
	if ($counter >2) {
		my @out_file_split = split (/\_/, $out_file);
		pop @out_file_split;
		$out_file = join "_", @out_file_split;
		}
	$out_file = $out_file . "_" . $counter;
	next;
	}


open INFILE, $bowtie_file or die $!;


while (<LIBFILE>) 
{
	chomp $_;
	if ($_ =~ />/) 
	{
		$next_name = $_;
		if ($sequence_counter ==0) {$current_name = $next_name};
		++ $sequence_counter;
		$_="";
		}
	if (defined $next_name and $sequence_counter > 1) {
		@splitname = split />/, $current_name;
		@splitname2 = split / /, $splitname[1];
		chomp $splitname2[0];
		$seq_name[$sequence_counter-1] = $splitname2[0];
		my $first_ref = \$seq_name[$sequence_counter-1];
		$name_hash{$$first_ref} = $sequence_counter-1;
		$seq_length[$sequence_counter-1] = length($in_line);
		
		$total_bins[$sequence_counter-1] =  int($seq_length[$sequence_counter-1]/$bin_size);
		
		$i=0;
		while ($i <= $total_bins[$sequence_counter-1]) {
			$pos_count[$i][$sequence_counter-1]= 0; 	
			$neg_count[$i][$sequence_counter-1]= 0; 	
			++ $i;
		}
		$current_name = $next_name;
		undef $next_name;
		undef $in_line;
		
	
	}
	$in_line = $in_line . $_;
}

# make sure that single sequence files are also initialized

@splitname = split />/, $current_name;
@splitname2 = split / /, $splitname[1];
chomp $splitname2[0];
$seq_name[$sequence_counter] = $splitname2[0];
$name_hash{$seq_name[$sequence_counter]} = $sequence_counter;
$seq_length[$sequence_counter] = length($in_line);

$total_bins[$sequence_counter] =  int($seq_length[$sequence_counter]/$bin_size);

$i=0;
while ($i <= $total_bins[$sequence_counter]) {
	$pos_count[$i][$sequence_counter]= 0; 	
	$neg_count[$i][$sequence_counter]= 0; 	
	++ $i;
	}
	
		
		

close LIBFILE;





# process INFILE
while (<INFILE>) 
	{
	@items = split(/\t/, $_);
		$bin_number = int($items[3]/$bin_size);
	if ($items[1] =~ /\-/)
		{
		++ $neg_count[$bin_number][$name_hash{$items[2]}];
		}
	if ($items[1] =~ /\+/)
		{
		++ $pos_count[$bin_number][$name_hash{$items[2]}];
		}
     }
close(INFILE);


$i = 1;
while ($i <= $sequence_counter) {
	my $outfile_name="$out_file"."_"."$seq_name[$i]\.binned";
	open OUTFILE, "> $outfile_name";
	# print OUTFILE "sequence name: \t", $seq_name[$i], "\t sequence length: \t", $seq_length[$i], "\n";
	$bin_number = 0;
	print OUTFILE $seq_name[$i], "\t";
	my $position = $bin_size * $bin_number;
	print OUTFILE $position, "\t";
	print OUTFILE $pos_count[$bin_number][$i], "\t";	
	print OUTFILE "-", $neg_count[$bin_number][$i], "\t";
	print OUTFILE ($pos_count[$bin_number][$i]+$neg_count[$bin_number][$i]), "\n";	
	++ $bin_number;
	while ($bin_number <= $total_bins[$i])
		{
		print OUTFILE $seq_name[$i], "\t";
		my $position = $bin_size * $bin_number;
		print OUTFILE $position, "\t";
		print OUTFILE $pos_count[$bin_number][$i], "\t";	
		print OUTFILE "-", $neg_count[$bin_number][$i], "\t";
		print OUTFILE ($pos_count[$bin_number][$i]+$neg_count[$bin_number][$i]), "\n";	
		++ $bin_number;
		}
	# print OUTFILE "binning window size: \t", $bin_size,"\n\n";
	++$i;
	}
close (OUTFILE);





 







