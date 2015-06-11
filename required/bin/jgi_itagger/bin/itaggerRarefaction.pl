#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use Statistics::Descriptive;

my $usage=<<'ENDHERE';
NAME:
itaggerRarefaction.pl

PURPOSE:
OTU tables rarefaction wrapper. Takes in input multiple OTU tables and
rarefy it to the lowest value of all the OTU tables. Implements the
single_rarefaction.py from Qiime.

INPUT:
--infile <string>          : OTU table in Qiime format.
                             Can be multiple infile.
--threshold <float>        : Fraction of average at which a sample
                             is considered failed. 

OUTPUT:
--outfile <string>         : Filtered otu table in Qiime format
                             Must be the same number if outfile
                             as infiles.

NOTES:
Please make sure that the two header lines of the OTU table is 
properly formatted. For instance:

#Full OTU Counts
#OTU ID 01  03  04  05  06  07  08  09  10  11  Consensus lineage	

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, @infile, @outfile, $begin, $end, $threshold);
my $verbose = 0;

GetOptions(
    'infile=s' 		=> \@infile,
	'outfile=s' 	=> \@outfile,
	'threshold=f' 	=> \$threshold,
	'begin=s' 		=> \$begin,
	'end=s' 		=> \$end,
    'verbose' 		=> \$verbose,
    'help' 			=> \$help
);
if ($help) { print $usage; exit; }

## VALIDATE
die("--infile file required\n") unless @infile;
die("--outfile outfile required\n") unless @outfile;
die("must be equal number of --infile and --outfile args...\n") if(@infile != @outfile);

$threshold = 0.05 unless($threshold);

## MAIN

# first pass through all infiles to find the lowest abundance value.
my $max = 0;
my $min = 999999999;
my @final_sums;
foreach my $infile (@infile){
	open(IN, "<".$infile) or die "Can't open infile ".$infile."\n";
	my @sums;
	my $counter = 0;	

	while(<IN>){
		chomp;
		if($_ =~ m/#/){
			next;
		}
		
		my @row = split(/\t/, $_);
		pop(@row); shift(@row);
	
		# For first row
		if($counter == 0){
			for(my $i=0; $i<@row; $i++){
				$sums[$i] = $row[$i];
			}	
				
			$counter++;
			next;
		}
	
		for(my $i=0; $i<@row; $i++){
			if(defined($sums[$i])){
				$sums[$i] = $sums[$i] + $row[$i];		
			}else{
				print "Not defined...at indice ".$i." file: ".$infile." line ".$.."\n";
				#die ?
				#$sums[$i] = 0;
				#$sums[$i] = $sums[$i] + $row[$i];		
			}
		}
	
		$counter++;
	}
	#Push sums of these samples into the final sums array.
	foreach my $sum (@sums){
	#	print $sum."\t";
		push(@final_sums, $sum);
	}
	#print "========================\n";
	close(IN);
}

#For Debug...
foreach(@final_sums){
#	print $_."\n";
}

# pick minimal value.
# Find a way to not pick min value if that min value is an outlier.
# For instance, a sample that failed and just have a total of 4 reads
# when the average of sample abundance is say 1500 reads.
@final_sums = sort @final_sums;
my $stat = Statistics::Descriptive::Sparse->new();
foreach my $value (@final_sums){
	$stat->add_data($value);
}
my $avg = $stat->mean();

# Then only consider a minimal value to be good if it is above at least 10% of the average sample abundance.
# This is to reject failed samples.
foreach my $value (@final_sums){
	$min = $value if( ($value < $min) && ($value > $avg * $threshold) );
}
print "Performing rarefaction at baseline <".$min."> for each tables.\n";		

# Then, normalize all tables upon that value.
foreach my $infile (@infile){
	my $outfile = shift(@outfile);
	# Run single_rarefactions.py 
	my $cmd = "single_rarefaction.py";
	$cmd .= " -i ".$infile;
	$cmd .= " -o ".$outfile;
	$cmd .= " -d ".$min;
	system($cmd);
}
exit;
