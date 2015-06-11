#!/usr/bin/env perl

use strict;
use warnings;
use iTagger::FastqDb::Tiny;

die("Usage: $0 <itaggerOutputDir/data> > allReads.fasta\n") unless @ARGV == 1;
my $dir = shift @ARGV;
opendir(DIR, $dir) or die($!);
my @libs = grep { $_ !~ /^\./ } readdir(DIR);
closedir(DIR);

my $rec;
foreach my $lib (@libs)
{
    my $n = 0;
    my $file = "$dir/$lib/hiqual.fastq";
    my $db = new iTagger::FastqDb::Tiny($file) or die("File not found: $file\n");
    print '>', $lib, '.', ++$n, "\n", $rec->seq, "\n" while $rec = $db->next_seq;
    print STDERR "$lib : $n seqs\n";
}
