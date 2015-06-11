=pod

=head1 NAME

iTagger::RdpClassifier

=head1 DESCRIPTION

This wrapper parallelizes RDP Classifier and generates and environmental
community file from RDP results and OTU Fasta with header-encoded abundances
(in iTagger format).

=head1 FUNCTIONS

=over 5

=cut

package iTagger::RdpClassifier;

use strict;
use warnings;
use Carp;
use threads;
use POSIX qw/mkfifo/;
use Env qw(TMPDIR RDP_JAR_PATH);
use iTagger::FastaDb;
use iTagger::FastqDb;

our $VERSION = 1.0;

our $tmpdir = -e '/scratch' ? '/scratch' : $TMPDIR;

=item classify

Run RDP Classifier and generate environmental community table. Default 8 threads.

=cut

sub classify
{
    my ($inFile, $trainingFile, $rdpOutFile, $envComOutFile, $envComFailOutFile, $libNames, $minWords, $level, $cutoff, $threads) = @_;
    confess("Missing parameters\n") unless $inFile and -e $inFile and $trainingFile and -e $trainingFile and $rdpOutFile and $envComOutFile and $libNames and @$libNames and $cutoff;
    confess("Invalid cutoff, $cutoff\n") if $cutoff < 0 or $cutoff > 1;
    $envComFailOutFile = '/dev/null' unless $envComFailOutFile;
    $cutoff = 0 unless defined($cutoff);
    $minWords = 120 unless $minWords;
    $threads = 8 unless $threads;
    my %levelIndex = 
    (
        'kingdom' => 0, 'phylum'  => 1, 'class'   => 2, 
        'order'   => 3, 'family'  => 4, 'genus'   => 5
    );
    if ( defined($level) )
    {
        $level = lc($level);
        confess("Invalid level $level\n") unless exists($levelIndex{$level});
    }

    # SPLIT INPUT (RDP CLASSIFIER WON'T READ FROM A FIFO)
    my @inFiles = map { "$tmpdir/$$.rdp.$_.fasta" } (1..$threads);
    my @fhs;
    my @otus = ();
    for ( my $i=0; $i<$threads; $i++)
    {
        open( my $fh, '>', $inFiles[$i] ) or confess($!);
        push @fhs, $fh;
    }
    open(my $in, '<', $inFile) or confess($!);
    my $line = <$in>;
    close($in);
    if ($line =~ /^>/) { $in = new iTagger::FastaDb($inFile) }
    elsif ($line =~ /^@/) { $in = new iTagger::FastqDb($inFile) }
    else { confess("Expected infile in Fasta or Fastq format\n") }
    my $i = -1;
    while ( my $rec = $in->next_seq )
    {
        $i = 0 if ++$i == $threads;
        my $fh = $fhs[$i];
        print $fh $rec->output;
        confess("Bad header: ".$rec->id."\n") unless $rec->id =~ /^(\d+);size=\d+;obs=([\d:]+)$/;
        my $index = $1;
        my @obs = split(/:/, $2);
        $otus[$index] = \@obs;
    }
    while ( my $fh = shift @fhs ) { close $fh }

    # RUN RDP CLASSIFIER
    my @fifos = map { "$tmpdir/$$.rdp.$_.tsv" } (1..$threads);
    foreach my $fifo (@fifos)
    {
        POSIX::mkfifo($fifo, 0600) or confess("Cannot mkfifo, $fifo: $!\n");
    }
    my @threads;
    for ( my $i=0; $i<$threads; $i++)
    {
        push @threads, new threads(\&_runRdp, $RDP_JAR_PATH, $inFiles[$i], $fifos[$i], $trainingFile, $minWords);
    }
        
    # GATHER OUTPUT
    @fhs = ();
    foreach my $fifo (@fifos)
    {
        open( my $fh, '<', $fifo ) or confess($!);
        push @fhs, $fh;
    }
    open(my $out, '>', $rdpOutFile) or confess($!);
    open(my $pass, '>', $envComOutFile) or confess($!);
    open(my $fail, '>', $envComFailOutFile) or confess($!);
    print $pass join("\t", '#OTU', @$libNames, 'Consensus lineage'), "\n";
    print $fail join("\t", '#OTU', @$libNames, 'Consensus lineage'), "\n";
    my $done = 0;
    while (!$done)
    {
        foreach my $fh (@fhs)
        {
            my $line = <$fh>;
            if( defined($line) )
            {
                print $out $line;
                my @row = split(/\t/, $line);
                my $hdr = $row[0];
                confess unless $hdr =~ /^(\d+);size=\d+;obs=[\d:]+$/;
                my $index = $1;
                my @taxon = map { $_ // 0 } @row[5,8,11,14,17,20];
                my @prob  = map { $_ // 0 } @row[7,10,13,16,19,22];
                my $ok = 1;
                if ( defined($level) )
                {
                    if ( $ok = $prob[$levelIndex{$level}] >= $cutoff ? 1 : 0 )
                    {
                        @taxon = @taxon[0..$levelIndex{$level}];
                        @prob  = @prob[0..$levelIndex{$level}];
                    } # else write entire record to failed file
                } else
                {
                    my @best = ();
                    for (my $i=0; $i<=$#taxon; $i++)
                    {
                        push @best, $taxon[$i] if $prob[$i] >= $cutoff;
                    }
                    @taxon=@best;
                }
                my $taxon = join(';', @taxon);
                my $line = join("\t", $index, @{$otus[$index]}, $taxon). "\n";
                if ( $ok ) { print $pass $line } else { print $fail $line }
            } else
            {
                $done = 1;
                last;
            }
        }
    }
    close($out);
    close($pass);
    close($fail);
    foreach my $fh (@fhs) { close($fh) }
    while ( my $thr = shift @threads ) { $thr->join }
    unlink(@inFiles, @fifos);
}

=item _runRdp

Run RDP Classifier

=cut

sub _runRdp
{
    my ($rdp, $in, $out, $train, $min) = @_;
    #system("java -Xmx1g -jar $rdp -q $in -o $out -t $train --minWords $min");
    system("java -jar $rdp -q $in -o $out -t $train --minWords $min");
    confess unless $? == 0;
}

1;

=back

=head1 AUTHORS

This wrapper by Edward Kirton, Julien Tremblay.
RDP Classifier by James R. Cole, Qiong Wang, James M. Tiedje.

=head1 COPYRIGHT/LICENSE

Copyright (c) 2013 US DOE Joint Genome Institute.  Use freely under the same license as RDP Classifier itself.

=cut
