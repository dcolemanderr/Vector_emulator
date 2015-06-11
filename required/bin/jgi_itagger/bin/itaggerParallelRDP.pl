#!/usr/bin/env perl

=pod

=head1 NAME

itaggerParallelRDP.pl

=head1 DESCRIPTION

Run RDP Classifier in parallel using threads.

=head1 OPTIONS

=over 5

=item @ARGV

Specify columns (library names) in appropriate order.

=item --in C<path>

Input file in Fasta format with headers in iTagger format.

=item --train C<path>

RDP training file.

=item --rdp C<path>

RDP output in tabular format.

=item --envcom C<path>

Environmental community (abundance) file in tabular format.

=item --min C<int>

Minimum words.

=item --level C<string>

Taxonomic level at which to report classification (default="genus").

=item --cutoff C<float>

Threshold for deciding which OTUs to filter from output (0..1; default = 0.5).

=item --threads C<int>

Number of threads (default=8).

=back

=head1 AUTHORS

Wrapper by Edward Kirton and Julien Tremblay
RDP Classifier by James R. Cole, Qiong Wang, James M. Tiedje

=head1 COPYRIGHT/LICENSE

Use freely under the same license as RDP Classifier itself.

=cut

use strict;
use warnings;
use Env qw(TMPDIR RDP_JAR_PATH);
use Getopt::Long;
use Pod::Usage;
use iTagger::RdpClassifier;

our $tmpdir = -e '/scratch' ? "/scratch/$$" : "$TMPDIR/$$";
-d $tmpdir or mkdir($tmpdir) or die($!);
$SIG{INT} = sub { system("rm -rf $tmpdir") };

my ( $help, $man, $infile, $rdpOutfile, $envComOutfile, $trainFile );
my $minWords = 120;
my $threads = 8;
my $level = 'genus',
my $cutoff = 0.5;
GetOptions(
	'in=s'      => \$infile,
	'train=s'   => \$trainFile,
	'min=i'     => \$minWords,
	'rdp=s'     => \$rdpOutfile,
	'envcom=s'  => \$envComOutfile,
    'level=s'   => \$level,
    'cutoff=f'  => \$cutoff,
	'threads=i' => \$threads,
	'help'      => \$help,
	'man'       => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

die("\$RDP_JAR_PATH not defined\n") unless $RDP_JAR_PATH;
die("--in file arg required\n")     unless ($infile);
die("--rdp file arg required\n")    unless ($rdpOutfile);
die("--envcom file arg required\n") unless ($envComOutfile);
die("--train file required\n")      unless ($trainFile);
die("Library names required\n")     unless (@ARGV);
die("--minWords must be <= 200 \n") if ( $minWords > 200 );
die("--level required\n")           unless $level;
die("--cutoff must be from 0-1\n")  if $cutoff < 0 or $cutoff > 1;

iTagger::RdpClassifier::classify( $infile, $trainFile, $rdpOutfile, $envComOutfile, "$envComOutfile.FAIL", \@ARGV, $minWords, $level, $cutoff, $threads );
exit;
