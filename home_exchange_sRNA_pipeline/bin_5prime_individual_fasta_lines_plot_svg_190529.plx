#!/usr/bin/perl -w
# quantify_in_bins.plx

use Chart::Gnuplot;

# use warnings;
use strict;

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
my $genome_matching=1;


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
	if ($_ =~ /-norm:/) {
		my @split_entry = split (/:/, $_);
		$genome_matching = $split_entry[1];
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
# $out_file = $out_file . "_bin" . $bin_size;



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
	my $read_length = length($items[4])/$bin_size; # to convert position: for 5'-end, the antisense read will need this, for 3'-position the sense read
	if ($items[1] =~ /\-/)
		{
		++ $neg_count[$bin_number+ $read_length][$name_hash{$items[2]}];
		
		}
	if ($items[1] =~ /\+/)
		{
		++ $pos_count[$bin_number][$name_hash{$items[2]}];
		}
     }
close(INFILE);


$i = 0;
while ($i < $sequence_counter) {
	# replace inappropriate characters in sequence_name	
	$seq_name[$i] =~ s/\\/_/g;
	$seq_name[$i] =~ s/\|/_/g;
        $seq_name[$i] =~ tr/ +/_/;
	my $outfile_name="$out_file"."_"."$seq_name[$i]\.binned";
	open OUTFILE, "> $outfile_name";
	my $outfile_ppm="$out_file"."_"."$seq_name[$i]\.binppm";
	open OUTFILEPPM,"> $outfile_ppm";
	print OUTFILEPPM "genome-matching reads in library: \t", $genome_matching,"\n";
	my @x_val;
	my @y_val_sense;
	my @y_val_antisense;
	my @y_val_sense_ppm;
	my @y_val_antisense_ppm;
	my @data_plot;

	$bin_number = 0;
	print OUTFILE $seq_name[$i], "\t";
	my $position = $bin_size * $bin_number;
	$x_val[$bin_number]=$position;
	$y_val_sense[$bin_number]=$pos_count[$bin_number][$i];
	$y_val_sense_ppm[$bin_number]=($pos_count[$bin_number][$i] * 1000000 / $genome_matching);	
	$y_val_antisense[$bin_number]=-($neg_count[$bin_number][$i]);
	$y_val_antisense_ppm[$bin_number]=-($neg_count[$bin_number][$i] * 1000000 / $genome_matching);	
	print OUTFILE $position, "\t";
	print OUTFILE $pos_count[$bin_number][$i], "\t";	
	print OUTFILE "-", $neg_count[$bin_number][$i], "\t";
	print OUTFILE ($pos_count[$bin_number][$i]+$neg_count[$bin_number][$i]), "\n";	

	# outfile ppm
	print OUTFILEPPM $seq_name[$i], "\t";

	print OUTFILEPPM $position, "\t";
	my $pos_ppm = ($pos_count[$bin_number][$i]*1000000/$genome_matching);
	print OUTFILEPPM $pos_ppm, "\t";	
	my $neg_ppm = -($neg_count[$bin_number][$i]*1000000/$genome_matching);
	print OUTFILEPPM $neg_ppm, "\t";
	print OUTFILEPPM ($pos_ppm-$neg_ppm), "\n";

	++ $bin_number;
	while ($bin_number <= $total_bins[$i])
		{
		print OUTFILE $seq_name[$i], "\t";
		my $position = $bin_size * $bin_number;
		$x_val[$bin_number]=$position;
		$y_val_sense[$bin_number]=$pos_count[$bin_number][$i];
		$y_val_sense_ppm[$bin_number]=($pos_count[$bin_number][$i] * 1000000 / $genome_matching);	
		$y_val_antisense[$bin_number]=-($neg_count[$bin_number][$i]);
		$y_val_antisense_ppm[$bin_number]=-($neg_count[$bin_number][$i] * 1000000 / $genome_matching);			
		print OUTFILE $position, "\t";
		print OUTFILE $pos_count[$bin_number][$i], "\t";	
		print OUTFILE "-", $neg_count[$bin_number][$i], "\t";
		print OUTFILE ($pos_count[$bin_number][$i]+$neg_count[$bin_number][$i]), "\n";	
		++ $bin_number;
		
		# outfile ppm
		print OUTFILEPPM $seq_name[$i], "\t";
		$position = $bin_size * $bin_number;
		print OUTFILEPPM $position, "\t";
		my $pos_ppm = ($pos_count[$bin_number][$i] * 1000000/$genome_matching);
		print OUTFILEPPM $pos_ppm, "\t";	
		my $neg_ppm = -($neg_count[$bin_number][$i] * 1000000/$genome_matching);
		print OUTFILEPPM $neg_ppm, "\t";
		print OUTFILEPPM ($pos_ppm-$neg_ppm), "\n";

		}
	
	# create graph instance
	my $graphfile_name="$out_file"."_"."$seq_name[$i]\.svg";
	my $graph = Chart::Gnuplot->new(output => $graphfile_name, title => $graphfile_name, ylabel =>"mapped reads (absolute counts)", terminal =>"svg mousing size 1920,1080");
	
	my $dataSet_sense = Chart::Gnuplot::DataSet->new(xdata =>\@x_val, ydata =>\@y_val_sense, style =>"lines", color =>"black", linetype =>1);
	my $dataSet_antisense = Chart::Gnuplot::DataSet->new(xdata =>\@x_val, ydata =>\@y_val_antisense, style =>"lines", color => "red", linetype => 1);
	$graph->plot2d($dataSet_sense, $dataSet_antisense);
	# create graph as ppm
	my $graphfile_ppm="$out_file"."_"."$seq_name[$i]_ppm\.svg";
	my $graphfile_png="$out_file"."_"."$seq_name[$i]_ppm\.png";
	my $graph_ppm = Chart::Gnuplot->new(output => $graphfile_ppm, title => $graphfile_ppm, ylabel =>"mapped reads (ppm)", terminal  => "svg mousing size 1920,1080");
	my $graph_png = Chart::Gnuplot->new(output => $graphfile_png, title => $graphfile_png, ylabel =>"mapped reads (ppm)", terminal =>"png size 1280,1024");
	my $dataSet_sense_ppm = Chart::Gnuplot::DataSet->new(xdata =>\@x_val, ydata =>\@y_val_sense_ppm, style =>"lines", color =>"black", linetype =>1);
	my $dataSet_antisense_ppm = Chart::Gnuplot::DataSet->new(xdata =>\@x_val, ydata =>\@y_val_antisense_ppm, style =>"lines", color => "red", linetype => 1);
	$graph_ppm->plot2d($dataSet_sense_ppm, $dataSet_antisense_ppm);
	$graph_png->plot2d($dataSet_sense_ppm, $dataSet_antisense_ppm);
	close (OUTFILE);	
	close (OUTFILEPPM);
	
	# next individual element
	++$i;
	}






 







