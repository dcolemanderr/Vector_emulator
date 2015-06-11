#!/usr/bin/env perl

=pod

=head1 NAME

itaggerClustering.pl

=head1 DESCRIPTION

Cluster sequences

=head1 OPTIONS

=over 5

=item --config C<path>

Configuration file

=item --dir C<path>

Base output folder.

=back

=head1 AUTHOR/SUPPORT

Edward Kirton (ESKirton@LBL.gov),
Julien Tremblay (JTremblay@lbl.gov),
usearch by R.C. Edgar.

=cut

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(rel2abs);
use iTagger::ClusterOTUs qw(clusterOTUs);

my ($help, $man, $configFile, $dir);
my $threads = 8;
my $verbose = 0;
GetOptions(
    'config=s'  => \$configFile,
    'dir=s'     => \$dir,
    'threads=i' => \$threads,
    'verbose'   => \$verbose,
	'help'      => \$help,
    'man'       => \$man
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
die("--config file required\n") unless $configFile;
die("--dir path required\n") unless $dir;
clusterOTUs($configFile, $dir, $threads, $verbose);
exit;
