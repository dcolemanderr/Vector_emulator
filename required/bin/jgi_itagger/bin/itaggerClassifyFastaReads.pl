#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use File::Basename;
use File::Which;
use File::Spec::Functions qw(rel2abs);

my $usage=<<'ENDHERE';
NAME:
itaggerClassifyFastaReads.pl

PURPOSE:
Classify fasta sequences and provides taxonomy.
No QC is done.

INPUT:
--infile <string>   : fasta file

OUTPUT:
--outdir <string>  : outdir where results files are written

NOTES:

BUGS/LIMITATIONS:

AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov
ENDHERE

## OPTIONS
my ($help, @infile, $outdir);
my $verbose = 0;

## SCRIPTS
GetOptions(
    'infile=s' 	=> \@infile,
	'outdir=s' 	=> \$outdir,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help
);
if ($help) { print $usage; exit; }

#VALIDATE
die("--infile file required\n") unless @infile;
die("--outdir outdir required\n") unless $outdir;

my $bin_dir_name = dirname(rel2abs($0));
$ENV{PATH} = $bin_dir_name.":".$ENV{PATH};

################
##### MAIN #####
################

my $rdp_parallel = which('itaggerParallelRDP.pl');
if(!defined($rdp_parallel)) {
    die("itaggerParallelRDP.pl is not on the path. Quitting.");
}else{
}

my $rdp_wrapper = which('itaggerRDPWrapper.pl');
if(!defined($rdp_parallel)) {
    die("itaggerParallelRDP.pl is not on the path. Quitting.");
}else{
}

my $dummy_seqobs_from_fasta = which('itaggerDummyObsTable.pl');
if(!defined($dummy_seqobs_from_fasta)) {
    die("itaggerDummyObsTable.pl is not on the path. Quitting.");
}else{
}

my $add_tax = which('itaggerAddTaxToSeqobs.pl');
if(!defined($add_tax)) {
    die("itaggerAddTaxToSeqobs.pl is not on the path. Quitting.");
}else{
}

my $summarize_taxa = which('summarize_taxa.py');
if(!defined($summarize_taxa)) {
    die("summarize_taxa.py is not on the path. Quitting.");
}else{
}

my $graph_phylum = which('itaggerPhylumBarplot.pl');
if(!defined($graph_phylum)) {
    die("itaggerPhylumBarplot.pl is not on the path. Quitting.");
}else{
}

my $rdp_training_set_id = "GG_full_length_mito_chloro_euk";
my $compute_cluster_group = "prok-meco.p";
my $minWords = 120;
my $num_threads = 4;
my $start_at = 0;

mkdir $outdir unless -d $outdir;

my @obs_list;
my @rdp_list;
my @otu_tables;
my @prefixes;
my @tax_absolute_L1;
my @tax_absolute_L2;
my @tax_absolute_L3;
my @tax_absolute_L4;
my @tax_absolute_L5;
my @tax_absolute_L6;
my @tax_relative_L1;
my @tax_relative_L2;
my @tax_relative_L3;
my @tax_relative_L4;
my @tax_relative_L5;
my @tax_relative_L6;

foreach(@infile){
	
	my $prefix = basename($_);
	push(@prefixes, $prefix);
	## RDP WRAPPER
#	my $cmd = $rdp_parallel;
#	$cmd .= " --infile ".$_;
#	$cmd .= " --RDP_training_set ".$rdp_training_set_id;
#	$cmd .= " --outfile ".$outdir."/rdp_".$prefix.".tab";
#	$cmd .= " --minWords ".$minWords;
#	$cmd .= " --num_threads ".$num_threads;
#	wl($log_r, "Running RDP parallel wrapper : ".$cmd);
#	system($cmd) if($start_at <= 1);
#	$? != 0 ? die "command failed: $!\n" : print STDOUT "RDP wrapper (--no_clustering option) successfuly executed\n" if($verbose);

	my $cmd = $rdp_wrapper;
	$cmd .= " --compute_cluster_group ".$compute_cluster_group;
    $cmd .= " --fasta ".$_;
	$cmd .= " --query_number 300"; #Could be higher (i.e. 200 or 300)...
	$cmd .= " --RDP_training_set ".$rdp_training_set_id;
	$cmd .= " --outdir ".$outdir."/rdp_".$prefix."/";
	$cmd .= " --outfile ".$outdir."/rdp_".$prefix.".tab";
	$cmd .= " --log ".$outdir."/rdp_log.txt";
	$cmd .= " --verbose";
	system($cmd) if($start_at <= 1);
	$? != 0 ? die "command failed: $!\n" :  print "Command ".$cmd." successfuly executed\n" if($verbose);
	
	## GENERATE FAKE SEQOBS TABLES
	$cmd = $dummy_seqobs_from_fasta;
	$cmd .= " --fastq ".$outdir."/".$prefix.".fastq";
	$cmd .= " --outfile ".$outdir."/obs_".$prefix.".tab";
	system($cmd) if($start_at <= 2);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Generation of a dummy seqobs table (--no_clustering option) successfuly executed\n" if($verbose);

	## MERGE RDP AND SEQOBS > OTU TABLE ONE OTU TABLE FOR EACH SAMPLE; Add rdp taxonomy to SeqObs table
	$cmd = $add_tax;
	$cmd .= " --seqobs ".$outdir."/obs_".$prefix.".tab";
	$cmd .= " --rdp ".$outdir."/rdp_".$prefix.".tab";
	$cmd .= " --cutoff 0.50";
	$cmd .= " --outfile ".$outdir."/otu_table_".$prefix.".tab";
	$cmd .= " --outfile_failed ".$outdir."/otu_table_".$prefix.".tab";
	$cmd .= " --tax_level best ";
	system($cmd) if($start_at <= 3);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "OTU tables successfuly generated (--no_clustering option)\n";			

	## SUMMARIZE TAXA
    #Summarize taxonomy with absolute abundance
	for(my $i=1; $i<7; $i++){
		$cmd = $summarize_taxa;
		$cmd .= " -i ".$outdir."/otu_table_".$prefix.".tab";
		$cmd .= " -L ".$i;
		$cmd .= " -o ".$outdir."/absolute/";
		$cmd .= " -a";
		system($cmd) if($start_at <= 4);
		$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L".$i." (absolute abundance) successfuly generated (--no_clustering option) \n";
		push(@tax_absolute_L1, $outdir."/absolute/otu_table_".$prefix."_L".$i.".txt" ) if($i == 1);	
		push(@tax_absolute_L2, $outdir."/absolute/otu_table_".$prefix."_L".$i.".txt" ) if($i == 2);	
		push(@tax_absolute_L3, $outdir."/absolute/otu_table_".$prefix."_L".$i.".txt" ) if($i == 3);	
		push(@tax_absolute_L4, $outdir."/absolute/otu_table_".$prefix."_L".$i.".txt" ) if($i == 4);	
		push(@tax_absolute_L5, $outdir."/absolute/otu_table_".$prefix."_L".$i.".txt" ) if($i == 5);	
		push(@tax_absolute_L6, $outdir."/absolute/otu_table_".$prefix."_L".$i.".txt" ) if($i == 6);	
	}
	push(@obs_list, $outdir."/obs_".$prefix.".tab");
	push(@rdp_list, $outdir."/rdp_".$prefix.".tab");
	push(@otu_tables, $outdir."/otu_table".$prefix.".tab");
}

if($start_at <= 41){
	my $j=0;
	
	## Have to fix an error in the tool that generate graphs. There is a problem when there is only one classification lineage present.
	#$cmd = $graph_phylum;
	#foreach(@tax_absolute_L1){
	#	$cmd .= " --infile ".$_;
	#	$cmd .= " --names ".$prefixes[$j];
	#	$j++;
	#}
	#$cmd .= " --outfile_graph ".$singleReadsTax_dir."/singleReads_taxonomy_L1.pdf";
	#$cmd .= " --outfile_table ".$singleReadsTax_dir."/singleReads_taxonomy_L1.tab";
	#print $cmd."\n";
	#system($cmd);
	#$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L1 compilation (barplots and tables) (--no_clustering option) successfully generated.\n";

	my $cmd = $graph_phylum;
	$j=0;
	foreach(@tax_absolute_L2){
		$cmd .= " --infile ".$_;
		$cmd .= " --names ".$prefixes[$j];
		$j++;
	}
	$cmd .= " --outfile_graph ".$outdir."/singleReads_taxonomy_L2.pdf";
	$cmd .= " --outfile_table ".$outdir."/singleReads_taxonomy_L2.tab";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L2 compilation (barplots and tables) (--no_clustering option) successfully generated.\n";
	
	$cmd = $graph_phylum;
	$j=0;
	foreach(@tax_absolute_L3){
		$cmd .= " --infile ".$_;
		$cmd .= " --names ".$prefixes[$j];
		$j++;
	}
	$cmd .= " --outfile_graph ".$outdir."/singleReads_taxonomy_L3.pdf";
	$cmd .= " --outfile_table ".$outdir."/singleReads_taxonomy_L3.tab";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L3 compilation (barplots and tables) (--no_clustering option) successfully generated.\n";

	$cmd = $graph_phylum;
	$j=0;
	foreach(@tax_absolute_L4){
		$cmd .= " --infile ".$_;
		$cmd .= " --names ".$prefixes[$j];
		$j++;
	}
	$cmd .= " --outfile_graph ".$outdir."/singleReads_taxonomy_L4.pdf";
	$cmd .= " --outfile_table ".$outdir."/singleReads_taxonomy_L4.tab";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L4 compilation (barplots and tables) (--no_clustering option) successfully generated.\n";

	$cmd = $graph_phylum;
	$j=0;
	foreach(@tax_absolute_L5){
		$cmd .= " --infile ".$_;
		$cmd .= " --names ".$prefixes[$j];
		$j++;
	}
	$cmd .= " --outfile_graph ".$outdir."/singleReads_taxonomy_L5.pdf";
	$cmd .= " --outfile_table ".$outdir."/singleReads_taxonomy_L5.tab";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L5 compilation (barplots and tables) (--no_clustering option) successfully generated.\n";

	$cmd = $graph_phylum;
	$j=0;
	foreach(@tax_absolute_L6){
		$cmd .= " --infile ".$_;
		$cmd .= " --names ".$prefixes[$j];
		$j++;
	}
	$cmd .= " --outfile_graph ".$outdir."/singleReads_taxonomy_L6.pdf";
	$cmd .= " --outfile_table ".$outdir."/singleReads_taxonomy_L6.tab";
	system($cmd);
	$? != 0 ? die "command failed: $!\n" : print STDOUT "Taxonomic summary L6 compilation (barplots and tables) (--no_clustering option) successfully generated.\n";
}


exit;
