#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
#use Bio::SeqIO;
use Bio::DB::Query::GenBank;
use Bio::DB::GenBank;
#use FastaDb;

my $usage=<<'ENDHERE';
NAME:
    get_genbank_annotation.pl
PURPOSE:
	To blast a primer sequence (short seq) against a reference sequence. This is to get hits coordinates.
INPUT:
	--infile <fasta_file> : fasta file containing only one primer seq
OUTPUT:
	--outfile <tab_file> : file containing one best hit/ref_db seq.
NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
     Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, $infile, $outfile);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'primer=s' => \$infile,
	'outfile=s' => \$outfile,
    'verbose' => \$verbose,
    'help' => \$help
);
if ($help) { print $usage; exit; }

#VALIDATE

#MAIN=====================================================================================================================================================================

my $query = "327360444";
my $query_obj = Bio::DB::Query::GenBank->new(-db => 'nucleotide', -query => $query );

my $gb_obj = Bio::DB::GenBank->new;
 
my $stream_obj = $gb_obj->get_Stream_by_query($query_obj);
 

while (my $seq_object = $stream_obj->next_seq) {  
	   
	   my $species_object = $seq_object->species;
	   my $species_string = $species_object->node_name;
	    
	   # Perlish
	   #my $species_string = $seq_object->species->node_name;
	   # either way, $species_string is "Homo sapiens"
	    
	   # get all taxa from the ORGANISM section in an array
		my @classification = $seq_object->species->classification;
		foreach $_ (@classification){
			print $_."\n";
		}
	   
}

exit;
