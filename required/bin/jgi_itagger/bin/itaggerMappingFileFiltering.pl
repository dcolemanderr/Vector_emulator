#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use iTagger::FastaDb;
use iTagger::FastqDb;
use Data::Dumper;

my $usage=<<'ENDHERE';
NAME:
itaggerFilterObsTable.pl

PURPOSE:
Will remove from an OTU table all the samples appearing in the
mapping file provided.

INPUT:
--otu_table <string>        : OTU table.
--map <string>              : Mapping file.

OUTPUT:
--outfile <string>          : OTU table.
--outfile_filtered <string> : OTU table.

NOTES
This script can be useful if you have a sample in an OTU table
that contains significantly less reads count compared to
other samples. Once these offending samples have been removed,
rarefaction can be performed without rarefying at an 
exessively low read count value.

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $otu_table, $outfile_filtered, $map, $outfile);
my $verbose = 0;

## SCRIPTS
GetOptions(
	'otu_table=s'			=> \$otu_table,
	'outfile_filtered=s'	=> \$outfile_filtered,
	'map=s'					=> \$map,
	'outfile=s'				=> \$outfile,
    'verbose' 				=> \$verbose,
    'help' 					=> \$help
);
if ($help) { print $usage; exit; }

#VALIDATE
die "--otu_table missing...\n" unless($otu_table);
die "--map missing...\n" unless($map);
die "--outfile missing...\n" unless($outfile);

## MAIN
open(IN, "<".$map) or die "Can't open file for reading ".$map."\n";
my %hash_map;
while(<IN>){
	chomp;
	next if($_ =~ m/#/);
	my @row = split(/\t/, $_);
	$hash_map{$row[0]} = "";
}
close(IN);
#print Dumper(\%hash_map);

open(OUT, ">".$outfile) or die "Can't open file ".$outfile."\n";
open(OUT_F, ">".$outfile_filtered) or die "Can't open file ".$outfile_filtered."\n";
open(IN, "<".$otu_table) or die "Can't open file ".$otu_table."\n";
my %hash_count;
while(<IN>){
	chomp;
	if($_ =~ m/#/){
		if($_ =~ /#OTU ID/){
			print OUT "#OTU ID";
			print OUT_F "#OTU ID";
			my @row = split(/\t/, $_);
			my $i = 0;
			pop(@row); shift(@row);
			foreach my $row (@row){
				if(exists($hash_map{$row})){
					$hash_count{$i} = $i;
					print OUT_F "\t".$row;
				}else{
					print OUT "\t".$row;
				}
				$i++;
			}
			print OUT "\tConsensus lineage\n";
			print OUT_F "\tConsensus lineage\n";
			next;
		}else{
			print OUT $_."\n";
			print OUT_F $_."\n";
			next;
		}
		#print Dumper(\%hash_count);
	}
	
	# Once we have the indices of the columns we want to exlude. loop through OTU
	# table and bin columns accordingly.
	my @row = split(/\t/, $_);
	my $lineage = pop(@row);
	my $id = shift(@row);
	print OUT $id;
	print OUT_F $id;

	my $i = 0;
	foreach my $row(@row){

		if(exists $hash_count{$i}){
			print OUT_F "\t".$row;
		}else{
			print OUT "\t".$row;
		}
		$i++;
	}	
	print OUT "\t".$lineage."\n";
	print OUT_F "\t".$lineage."\n";
}

exit;
