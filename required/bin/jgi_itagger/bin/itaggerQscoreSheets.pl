#!/usr/bin/env perl

=pod

=head1 NAME

itaggerQscoreSheets.pl

=head1 DESCRIPTION

This script generates Qscore sheets using fastx_quality_stats tool. If reads have IDs following iTagger format (lib*integer), then each section is reported separately.

=head1 OPTIONS

=over 5

=item @ARGV

List of reads files (in FASTQ format).

=item --log C<path>

Log summary (optional)

=item --threads C<int>  

Number of threads.  Each infile will be processed in a separate thread.  Default=8.

=item --outfile C<path>

Path to output file.  If multiple infiles are provided, each will be processed separately and each section will be separated by ">>NAME" (header) and "//" (footer) strings. Optional; default=STDOUT.

    >><ID>
    .
    .
    .
    //

=back

=head1 AUTHORS

Julien Tremblay (julien.tremblay@mail.mcgill.ca).
Edward Kirton (eskirton@lbl.gov).
Quality scores by fastx_quality_stats (gordon@cshl.edu).

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Env qw/TMPDIR/;
use Parallel::ForkManager;
use File::Copy;
use IO::File;

$| = 1; # UNBUFFER OUTPUT

## OPTIONS
my ($help, $man, $outfile, $logfile);
my $num_threads = 8;
GetOptions(
	'outfile=s' 	=> \$outfile,
	'log=s' 		=> \$logfile,
    'threads=i'     => \$num_threads,
	'help'          => \$help,
	'man'           => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
die("Infiles required\n") unless @ARGV;
die("\$TMPDIR not defined\n") unless $TMPDIR;
	
# GENERATE TABLE FOR EACH INFILE
my @infiles = @ARGV;
my @tmpfiles = ();
my $pm = new Parallel::ForkManager($num_threads);
foreach my $infile (@infiles)
{
    my ($basename) = fileparse($infile, qr/\.[^.]*/);
    my $tmpfile = "$TMPDIR/$basename.tsv";
    push @tmpfiles, $tmpfile;
    my $pid = $pm->start and next;
    system("fastx_quality_stats -i $infile -o $tmpfile.tmp");
    move("$tmpfile.tmp", $tmpfile);
    $pm->finish;
}
$pm->wait_all_children;

# COMBINE TMPFILES
if ( @tmpfiles > 1 )
{
    my $out;
    if ( $outfile )
    {
        $out = new IO::File(">$outfile") or die($!);
    } else
    {
        $out = *STDOUT;
    }
    foreach my $tmpfile (@tmpfiles)
    {
        my ($name) = fileparse($tmpfile, qr/\.[^.]*/);
        print $out '>>', $name, "\n";
        my $in = new IO::File("<$tmpfile") or die($!);
        while (<$in>) { print $out $_ }
        close($in);
        print $out '//', "\n";
    }
    close($out) if $outfile;
} elsif ( $outfile )
{
    move($tmpfile, $outfile) or di($!);
} else
{
    # PIPE OUTPUT
    my $in = new IO::File("<$tmpfile") or die($!);
    while (<$in>) { print }
    close($in);
}
unlink(@tmpfiles);
exit;
