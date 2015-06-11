#!/usr/bin/env perl

=pod

=head1 NAME

itaggerBinByLength.pl

=head1 DESCRIPTION

Simple script to bin sequences in a file by a range of lengths.

=head1 OPTIONS

=over 5

=item --in C<path>

Input file in Fastq or Fasta format (optional; default=STDIN).

=item --out C<path>

Directory where to store files of sequences binned by length.

=item coords C<string>

coords of seq length to bin.  For instance: C<--coords 0-290-345-360-375-385-400>

=back

=head1 AUTHOR

Julien Tremblay
Edward Kirton

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use iTagger::FastaDb;
use iTagger::FastqDb;

## OPTIONS
my ($help, $man, $infile, $coords, $outdir);
my $verbose = 0;

GetOptions(
    'infile=s'  => \$infile,
	'coords=s'	=> \$coords,
	'outdir=s'	=> \$outdir,
    'verbose' 	=> \$verbose,
    'help' 		=> \$help,
	'man'       => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

# VALIDATE
my @coords = split(/-/, $coords);
die("Invalid --coords: $coords\n") unless @coords;

# DETERMINE FORMAT
open(IN, '<', $infile) or die($!);
my $line = <IN>;
close(IN);
my $format;
if ( $line =~ /^@/ ) { $format = 'fastq' }
elsif ( $line =~ /^>/) { $format = 'fasta' }
else { die("Invalid infile\n") }

# OUTFILES
my @fhs;
for ( my $i=0; $i<=$#coords-1; $i++)
{
    my $coord = $coords[$i];
    my $outfile = "$outdir/$coord.$format";
    open(my $fh, '>', $outfile) or die($!);
	push(@fhs, $fh);
}

# MAIN
my $db = $format eq 'fastq' ? iTagger::FastqDb->new($infile) : iTagger::FastaDb->new($infile);
while( my $curr = $db->next_seq )
{
    for (my $i=0; $i<$#coords-1; $i++)
    {
        if ( $curr->len > $coords[$i] && $curr->len <= $coords[$i+1] )
        {
            my $fh = $fhs[$i];
            print $fh $curr->output();
            last;
        }
    }
}
close $_ foreach(@fhs);
exit;
