#!/usr/bin/env perl

=pod

=head1 NAME

itaggerReadQC.pl

=head1 DESCRIPTION

This script performs all steps required to process a library of Illumina reads prior to clustering.  Input file should contain all reads corresponding to one library (e.g. one barcode of one lane) in Fastq format (may be compressed).

=head1 OPTIONS

=over 5

=item --name C<string>

Short and unique identifier for this library/experimental condition.

=item --in C<path>

Raw reads for this library (one) in Fastq format; may be compressed (usually corresponds to one barcode of one lane of Illumina data).

=item --config C<path>

Config file in Ini format.  Refer to config file documentation for a description of the fields and their default values.

=item --dir C<path>

Output folder. Filenames are auto-generated for consistency and convenience.

=back

=head1 AUTHORS

Julien Tremblay (julien.tremblay@mail.mcgill.ca),
Edward Kirton (eskirton@lbl.gov)

=head1 COPYRIGHT

Copyright (c) 2013 US DOE Joint Genome Institute.  Use freely under the same license at Perl itself. Refer to wrapped tools for their own license and copyright information.

=cut

use strict;
use warnings;
use Env qw/TMPDIR RDP_JAR_PATH/;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(rel2abs);
use iTagger::ReadQC qw(libQC);

$SIG{INT} = sub { exit };    #Handle ungraceful exits with CTRL-C.

our $VERSION = "0.2";

# OPTIONS
my ( $help, $man, $libName, $configFile, $infile, $dir);
my $verbose = 0;
GetOptions(
    'name=s'   => \$libName,
    'config=s' => \$configFile,
    'infile=s' => \$infile,
    'dir=s' => \$dir,
    'verbose' => \$verbose,
	'help' => \$help,
	'man'  => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

die("--in file required\n") unless $infile;
die("--out dir required\n") unless $dir;
die("--config file required\n") unless $configFile;
die("--name required\n") unless $libName;

$libName =~ s/\W+/_/g;

$dir = rel2abs($dir);
-d $dir or mkdir($dir) or die($!);

$infile = rel2abs($infile);
$configFile = rel2abs($configFile);

libQC($configFile, $dir, $infile, $libName, $verbose);
exit;
