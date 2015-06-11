#!/usr/bin/env perl

#use JGI::Bio::Iterator::FastaDb;
use warnings;
use strict;
use Getopt::Long;

## USAGE

my $usage=<<'ENDHERE';
make_RDP.pl
PURPOSE:
NOTE:Should be run on gpint node, as genepool does not allow specifed amount of RAM for jave Virtual machine.
To make an RDP database from a correctly formatted fasta file. 
Sequences should contain no illegal characters, and headers should have the following format:
>MySeqID	Root;k__Kingdom;p__Phylum;c__Class;o__Order;f__family;g__genus;
The Kingdom,Phylum,Class,Order,Family,and Genus should be single words, no spaces, and only word and numeric characters.
They should NOT contain the following characters '(,),-,",',*,' or other non-word characters.

INPUT:
    --infile_fasta <fasta> : complete file
OUTPUT:
    --outfile_name <string> : outfile name for RDP database
OPTIONS:
AUTHOR/SUPPORT:
Devin Coleman-Derr dacoleman-derr@lbl.gov
ENDHERE

## OPTIONS

my ($help, $infile_fasta, $outfile_name);
GetOptions(
    'infile_fasta=s' => \$infile_fasta,
    'outfile_name=s' => \$outfile_name,
    'h|help' => \$help
);

## VALIDATE 

if ($help) { print $usage; exit; }
die("Fasta input file required\n") unless $infile_fasta;
die("Output filename required\n") unless $outfile_name;

## MAIN

open(FASTA, $infile_fasta);
while (<FASTA>) {
 
    if($_ =~ /^>/) {
    	
	# Checks Fasta header for correct format and illegal characters.
	die("The Sequence ID or classification in $_ are not formatted correctly. >MySeqID	Root;k___Kingdom;p___Phylum;c___Class;o___Order;f___family;g___genus;")
		unless $_ =~ /^>[\w\d]*\sRoot;k___[\w\d]*;p___[\w\d]*;c___[\w\d]*;o___[\w\d]*;f___[\w\d]*;g___[\w\d]*$/;
    } else {
	
	#Checks Fasta sequences for illegal characters.
    	die("One of the fasta Sequences : $_ contains illegal characters.") unless $_ =~ /^[a-zA-Z]*$/;

    }

}

#Defines outfile names based on specified outfile name root.
my ($outfile_fasta) = "$outfile_name".".fasta";
my ($outfile_model,) = "$outfile_name".".txt";
my ($outfile_lineages) = "$outfile_name".".tab";

#Builds the tab, lineage, and fasta file for use in the RDPClassifier build process.
system( "/global/homes/t/tremblay/tools/common_tools/generate_RDP_taxonomy_from_fasta_ITS.pl --fasta $infile_fasta  --tax_level 6 --outfile_model $outfile_model --outfile_fasta $outfile_fasta --outfile_lineages $outfile_lineages --log log.txt");

#Builds the RNAClassifier files, including bergey
system( "java -Xmx9g -cp /house/homedirs/t/tremblay/misc_tools/rdp_classifier_2.3/rdp_classifier-2.3.jar edu/msu/cme/rdp/classifier/train/ClassifierTraineeMaker $outfile_model $outfile_fasta 1 version1 $outfile_name ./");

