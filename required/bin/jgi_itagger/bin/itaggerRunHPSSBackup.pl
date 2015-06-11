#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use iTagger::FastaDb;
use iTagger::FastqDb;
use Data::Dumper;
use File::Find;

my $usage=<<'ENDHERE';
NAME:
hpss_itags.pl

PURPOSE:
To back Itag runs.
***********************************************
*Only the JGI analyst who is responsible of ***
*Itags production runs should run this script.*
***********************************************

INPUT:
--indir <string>    : Directory where itags runs are to be backed-up.
--delete            : Set if Itag runs directories that have been
                      backed-up are to be deleted.
--htar              : Set this flag if htar is to be used for back-up.
--hsi               : Set this flag if hsi is to be used for back-up.
                      *SPECIFY AT LEAST --htar of --hsi.
--num_threads <int> : Number of threads to use. Default :1.

OUTPUT:
--outdir <string> : Directory on the HPSS system where you want to 
                    backup your Itag runs. 

NOTES:
This script is intended for the JGI Itags production data analysis operators.
You have to specify one in-directory and one out-directory. In the in-dir, 
the script will search for all directory having the following nomenclature:
<project_id>_<itag_run>. Example: 1234122_M10011.1 --indir has to be a path
on the /global/blabla/... system. --outdir has to be on the HPSS backup
server (A:/blabla/blabla/) Without the A:.

hpss_itags.pl <outfile_hsi_commands> <outfile_htar_commands> --indir <string> --outdir <string>
OR the bash wrapper script:
run_hpss_backup.pl --indir <string> --outdir <string>

BUGS/LIMITATIONS:
Itag runs are expected to grow in size in the future. So htar might not
be appropriate to backup files exceeding 64 GB... At the moment each
fastqs are compressed to .gz which helps. Will have to update 
the backup script accordingly in the futyre to deal with bigger files...
Maybe do a mix of hsi and htar?
 
AUTHOR/SUPPORT:
Julien Tremblay - jtremblay@lbl.gov

ENDHERE

## OPTIONS
my ($help, $indir, $outdir, $hsi, $htar, $delete, $num_threads);
my $verbose = 0;

GetOptions(
    'indir=s' 		=> \$indir,
	'outdir=s'		=> \$outdir,
	'delete'		=> \$delete,
	'hsi'			=> \$hsi,
	'htar'			=> \$htar,
	'num_threads=i' => \$num_threads,
    'verbose' 		=> \$verbose,
    'help' 			=> \$help
);
if ($help) { print $usage; exit; }

## Validate
die "Usage: hpss_itags.pl <outfile> --indir <string> --outdir <string>\n" if($ARGV[0] eq "--indir");
die "Usage: hpss_itags.pl <outfile> --indir <string> --outdir <string>\n" if($ARGV[0] eq "--outdir");
die "--outdir arg missing\n" unless($outdir);
die "--indir arg missing\n" unless($indir);
die "Please specify either --hsi OR --htar\n" if(!$hsi and !$htar);

$num_threads = 1 unless($num_threads);

my $commands_hsi = $ARGV[0];
my $commands_htar = $ARGV[1];

## It is intended to recursively through a directory tree
## and compress all *.fastq in *.fastq.gz
sub eachFile{
	my $filename = $_;
	my $fullpath = $File::Find::name;
	#remember that File::Find changes your CWD, 
	#so you can call open with just $_

	if (-e $filename) { 

		if(substr($filename, -6) eq ".fastq"){
			print STDOUT "Compressing ".$fullpath." into .gz archive...\n";
			system("pigz -p ".$num_threads." ".$fullpath);
			$? != 0 ? die "command failed: $!\n" : print STDERR "Successfuly compressed ".$fullpath." into .gz archive...\n";
		}
	}
}

## MAIN

# Scan in-directory and build hash structure
my %hash;
my %hash_M; # Just store MXXXXX master directory names into that hash.

# Compress .fastq into .gz
find (\&eachFile, $indir);

# Then build directory structure into hash.
opendir (DIR, $indir) or die $!;
while (my $file = readdir(DIR)) {

	next unless (-d "$indir/$file");
	next if($file eq  ".." or $file eq ".");

	# Print current dir	
	print $indir."/".$file."\n" if($verbose);
	next if($file =~ m/_BACKEDUP/); # skip if directory has already been backed-up.
	if($file =~ m/(M\d+)/){
		$hash_M{$1} = $indir."/".$1;
		print $indir."/".$1."\n";
	}

	# And print its sub dirs.
	if(-d "$indir/$file"){
		opendir (DIR2, $indir."/".$file) or die $!;
		while (my $file_2 = readdir(DIR2)) {

			next unless (-d "$indir/$file/$file_2");
			next if($file_2 eq  ".." or $file_2 eq ".");

			print $indir."/".$file."/".$file_2."\n" if($verbose);

			# Parse file_2 nomenclature
			if($file_2 =~ m/(\d+)_(M\d+\.\d+)/){
				my $project_id = $1;
				my $itag_run_id = $2;	
				$hash{$itag_run_id}{project_id} 			= $project_id;
				$hash{$itag_run_id}{itag_run_id} 			= $itag_run_id;
				$hash{$itag_run_id}{Mxxxx_dir}	 			= $indir."/".$file;
				$hash{$itag_run_id}{in_path} 				= $indir."/".$file."/".$file_2;
				$hash{$itag_run_id}{out_path_project_id} 	= $outdir."/".$project_id;
				$hash{$itag_run_id}{out_path_itag_run_id} 	= $outdir."/".$project_id."/".$itag_run_id;
			}elsif($file_2 =~ m/(\d+)\.(M\d+\.\d+)/){
				my $project_id = $1;
				my $itag_run_id = $2;	
				$hash{$itag_run_id}{project_id} 			= $project_id;
				$hash{$itag_run_id}{itag_run_id} 			= $itag_run_id;
				$hash{$itag_run_id}{Mxxxx_dir}	 			= $indir."/".$file;
				$hash{$itag_run_id}{in_path} 				= $indir."/".$file."/".$file_2;
				$hash{$itag_run_id}{out_path_project_id} 	= $outdir."/".$project_id;
				$hash{$itag_run_id}{out_path_itag_run_id} 	= $outdir."/".$project_id."/".$itag_run_id;
			}	
		}
		closedir(DIR2);	
	}
}
closedir(DIR);

# print hash for debug
print "HASH DUMPING:\n" if($verbose);
print Dumper(\%hash) if($verbose);

if($hsi){
	open(HSI, ">".$commands_hsi) or die "Can't open file ".$commands_hsi;

	print HSI "\#!\/bin\/bash\n\n";
	print HSI "hsi -h archive.nersc.gov \"";

	# First, make directory if they don't exist
	foreach my $itag_run_id (keys %hash) {

		# Make project_id dir. with hsi, the mkdir command does not erase dir if it aleady exists.
		print HSI "mkdir ".$outdir."/".$hash{$itag_run_id}{project_id}."; ";

		# Make the itag_run_id 
		print HSI "mkdir ".$outdir."/".$hash{$itag_run_id}{project_id}."/".$hash{$itag_run_id}{itag_run_id}."; ";

		# Change local directory - No need to that in the end... left in comment just in case.
		#print HSI "lcd ".$hash{$project_id}{in_path}."; ";

		# Change HSI directory
		print HSI "cd ".$hash{$itag_run_id}{out_path_itag_run_id}."; ";

		# Copy everything in that local dir here. (use cput).
		print HSI "cput -R ".$hash{$itag_run_id}{in_path}."; ";
		print HSI "cput * ".$hash{$itag_run_id}{in_path}."/; ";
	}
	print HSI "\"";

	close(HSI);
}elsif($htar){

	# If htar, use hsi to first create directories.
	# Then use htar to copy the files into that dir

	print STDOUT "htarring directory...\n";	

	open(HSI, ">".$commands_hsi) or die "Can't open file ".$commands_hsi;
	open(HTAR, ">".$commands_htar) or die "Can't open file ".$commands_htar;

	# When bash script will run, exit if htar failed.
	print HTAR "\#!\/bin\/bash\n\n";
	print HTAR "set -e\n\n";

	print HSI "\#!\/bin\/bash\n\n";
	print HSI "hsi -h archive.nersc.gov \"";

	# First, make directory if they don't exist
	foreach my $itag_run_id (keys %hash) {

		# Make project_id dir. with hsi, the mkdir command does not erase dir if it aleady exists.
		print HSI "mkdir ".$outdir."/".$hash{$itag_run_id}{project_id}."; ";

		# make the itag_run_id. For now, for htar, don't create that dir. 
		#print HSI "mkdir ".$outdir."/".$project_id."/".$hash{$project_id}{itag_run_id}."; ";

		# change HSI directory
		print HTAR "htar -cf ".$hash{$itag_run_id}{out_path_itag_run_id}.".tar ";
		print HTAR $hash{$itag_run_id}{in_path}."\n";
	}
	print HSI "\"";

	# Then flag directories that have been backed-up as such by adding the "_BACKEDUP" string at the end of directory name.
	print HTAR "\n## Directories that have been renamed and can now be deleted:##\n" if(!$delete);
	print HTAR "\n## Deleting the following directories\n" if($delete);
	foreach my $itag_run_id (keys %hash_M) {
		if($delete){
			print HTAR "rm -rf ".$hash{$itag_run_id}{Mxxxx_dir}."\n";
		}else{
			my $cmd = "mv ".$hash_M{$itag_run_id}." ".$hash_M{$itag_run_id}."_BACKEDUP\n";
			print HTAR $cmd."\n";
		}
	}

	close(HSI);
	close(HTAR);

}else{
	die "Please specify either --hsi OR --htar\n";
}
exit;
