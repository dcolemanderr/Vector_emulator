#!/usr/bin/env perl

=pod

=head1 NAME

itaggerRdpClassifier.pl

=head1 DESCRIPTION

Run RDP Classifier with multiple threads and filter results.

=head1 OPTIONS

=over 5

=item --infile C<path>

Input file in FASTA format.  The header must be in the format of "<int>;size=<int>;obs=[<int>:<int>]+; where the first integer is the ID, the size is the total number of reads in the cluster (i.e. sum of observations); and the obs colon-separated list is the counts per library.

=item --training C<path>

RDP classifier training file (rRNAClassifier.properties)

=item --outfile C<path>

RDP outfile

=item --pass C<path>

Passed items file

=item --fail C<path>

Failed items file

=item --min C<int>

Minimum words (default=120)

=item --level C<string>

Taxonomic level for filtering (kingdom|phylum|class|order|family|genus)

=item --cutoff C<float>

Confidence cutoff level (0..1; default=0.5)

=item --threads C<int>

Number of processor threads to use (default=8).

=item @ARGV

A list of library names is required.  The order must match the obs list in the sequence header lines.

=back

=cut

use strict;
use warnings;
use iTagger::RdpClassifier;
use Pod::Usage;
use Getopt::Long;

our $VERSION = '0.1';

# OPTIONS
my ($help, $man, $inFile, $trainingFile, $rdpOutFile, $envComOutFile, $envComFailOutFile);
my $minWords = 120;
my $level = 'phylum';
my $cutoff = 0.5;
my $threads = 8;
GetOptions(
    'infile=s' => \$inFile,
    'training=s' => \$trainingFile,
    'outfile=s' => \$rdpOutFile,
    'pass=s' => \$envComOutFile,
    'fail=s' => \$envComFailOutFile,
    'min=i' => \$minWords,
    'level=s' => \$level,
    'cutoff=f' => \$cutoff,
    'threads=i' => $threads,
    'help|?' => \$help,
    'man' => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

die("Infile required\n") unless $inFile;
die("Training file required\n") unless $trainingFile;
die("RDP outfile required\n") unless $rdpOutFile;
die("Pass file required\n") unless $envComOutFile;
die("Fail file required\n") unless $envComFailOutFile;
die("Library names required\n") unless @ARGV;

iTagger::RdpClassifier::classify($inFile, $trainingFile, $rdpOutFile, $envComOutFile, $envComFailOutFile, \@ARGV, $minWords, $level, $cutoff, $threads);
exit;

__END__
Copyright (c) 2014 US DOE Joint Genome Institute.  Use freely under the same terms as RDP Classifier itself.
