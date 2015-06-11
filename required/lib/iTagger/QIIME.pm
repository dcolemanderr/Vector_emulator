=pod

=head1 NAME

iTagger::QIIME

=head1 DESCRIPTION

Functions for generating OTU reports, wrapping several QIIME scripts.

=head1 FUNCTIONS

=over 5

=cut

package iTagger::QIIME;

use strict;
use warnings;
use Carp;
use Env qw(TMPDIR);
use Date::Format;
use JSON;
use File::Copy;
use threads;
require Exporter;
use iTagger::Stats qw(trimmedStdDev);

our $VERSION = 1.0;
our @ISA = qw(Exporter);
our @EXPORT = qw(qiimeReport);
our $tmpDir = -e '/scratch' ? '/scratch' : $TMPDIR;

=item qiimeReport

Wraps several QIIME scripts

=cut

sub qiimeReport
{
    my ($configFile, $dir) = @_;
    confess("Missing args") unless $configFile and $dir;
    my $start = time;
    my $config = Config::Tiny->read($configFile);

    # CONVERT OTU TSV TO BIOM JSON
    my $inFile = "$dir/otu/otu.tax.filtered.tsv";
    my $biomFile = "$dir/otu/otu.biom";
    my $depth = otuToBiom($inFile, $biomFile);

    ## BEGIN THREADS
    my @threads;
    my $thr;

    # FILTER OTU TABLES AT ASSORTED LEVELS, CALC ALPHA DIVERSITY LEVELS
	($thr) = new threads(\&abundanceThreshold, $inFile, "$dir/abundance_thresholds");
    push @threads, $thr;

    # RAREFY OTU TABLES
    my $rarefiedFile = "$dir/otu/rarefied.$depth.biom";
    my $rarefiedFilteredFile = "$dir/otu/rarefied.$depth.filtered.biom";
    $thr = new threads(\&singleRarefaction, $biomFile, $rarefiedFile, $rarefiedFilteredFile, $depth);
    push @threads, $thr;
    my $rarefied1kFile = "$dir/otu/rarefied.1000.biom";
    my $rarefied1kFilteredFile = "$dir/otu/rarefied.1000.filtered.biom";
    $thr = new threads(\&singleRarefaction, $biomFile, $rarefied1kFile, $rarefied1kFilteredFile, 1000);
	push @threads, $thr;

    # SUMMARIZE TAXONOMY WITH ABSOLUTE AND RELATIVE ABUNDANCE
    my $relDir = "$dir/tax_mapping/relative";
    $thr = new threads(\&summarizeTaxa, $biomFile, $relDir);
    push @threads, $thr;
    my $absDir = "$dir/tax_mapping/absolute";
    summarizeTaxa($biomFile, $absDir, 1);

    # PLOT THE RANK-ABUNDANCE CURVE OF ALL SAMPLES
    $thr = new threads(\&plotRankAbundanceGraph, $biomFile, "$dir/otu/log_rank_abundance.pdf");
	push @threads, $thr;

    # MAKE TAXA PLOT ABSOLUTE ABUNDANCE
    my @files = map { "$dir/tax_mapping/absolute/otu_L$_.txt" } (2..6);
    my $outDir = "$dir/tax_mapping/plots";
    $thr = new threads(\&plotTaxaSummary, \@files, $outDir);
    push @threads, $thr;

    # JOIN THREADS
    while ( $thr = shift @threads) { $thr->join }

    # GENERATE A PHYLUM BARPLOT FOR A QUICK VIEW OF HOW THE DATA LOOKS
	_run("itaggerPhylumBarplot.pl --infile $dir/tax_mapping/absolute/otu_L2.txt --outFile_graph $dir/tax_mapping/taxonomy_phylum_L2.pdf --outFile_table $dir/tax_mapping/taxonomy_phylum_L2.tab --title Phylum");

    # FINISH
    my $dur = time - $start;
    print "QIIME done after $dur sec\n";
}

=item summarizeTaxa

Execute QIIME summarize_taxa.py

=cut

sub summarizeTaxa
{
    my ($inFile, $outDir, $abs) = @_;
    $abs = $abs ? '-a' : '';
    _run("summarize_taxa.py -i $inFile -o $outDir $abs");
}

=item plotRankAbundanceGraph

Execute QIIME plot_rank_abundance_graph.py

=cut

sub plotRankAbundanceGraph
{
    my ($inFile, $outFile) = @_;
	_run("plot_rank_abundance_graph.py -n -s '*' -i $inFile -o $outFile");
    unlink("${outFile}_log.txt");
}

=item plotTaxaSummary

Execute QIIME plot_taxa_summary.py

=cut

sub plotTaxaSummary
{
    my ($inFiles, $outDir) = @_;
    _run("plot_taxa_summary.py -i ".join(',', @$inFiles)." -o $outDir -c bar");
}

=item _run

Run an external executable, logconfess upon nonzero exit status.

=cut

sub _run
{
    my ($cmd) = @_;
    my $output = `$cmd 2>&1`;
    confess("FAILURE RUNNING $cmd : $output") unless $? == 0;
}

=item otuToBiom

Output biom json file.

=cut

sub otuToBiom
{
    my ($inFile, $outFile) = @_;
    confess("Missing args") unless $inFile and $outFile;

    # INIT JSON RECORD
    my $biom = 
    {
        id => undef,
        format => "Biological Observation Matrix 0.9.1-dev",
        format_url => "http://biom-format.org/documentation/format_versions/biom-1.0.html",
        type => "OTU table",
        generated_by => "iTagger revision $VERSION",
        date => time2str('%Y-%m-%dT%X', time),
        matrix_type => "dense",
        matrix_element_type => "int"
    };

    # LIBRARIES
    open ( my $in, '<', $inFile) or confess($!);
    my $hdr = <$in>;
    my @hdr = split(/\t/, $hdr);
    my @libs = @hdr[1..($#hdr-1)];
    for (my $i=0; $i<=$#libs; $i++)
    {
        $libs[$i] = { id => $libs[$i], metadata => {} };
    }
    $biom->{columns} = \@libs;
    my $numLibs = scalar(@libs);

    # OTUS JSON
    my @rows = ();
    my @data = ();
    my @libSizes = ( (0) x $numLibs );
    while (<$in>)
    {
        chomp;
        my @obs = split(/\t/);
        my $id = shift @obs;
        my $tax = pop @obs;
        my @tax = split(/;/, $tax);
        push @rows, { id => $id, metadata => { taxonomy => \@tax } };
        push @data, \@obs;
        for ( my $i=0; $i<=$#obs; $i++ ) { $libSizes[$i] += $obs[$i] }
    }
    close($in);
    my $numOtus = scalar(@rows);
    $biom->{rows} = \@rows;
    $biom->{data} = \@data;
    $biom->{shape} = [ $numOtus, $numLibs ];

    # WRITE
    open ( my $out, '>', $outFile) or confess($!);
    print $out to_json($biom);
    close($out);

    # CALCULATE CUTOFF OF RAREFACTION OF TRIMMED MEAN
    # (DISCARDING SMALLEST 5% AND LARGEST 5% OF VALUES)
    # (IF MORE THAN 2 LIBRARIES)
    return 0 unless @libSizes > 2;
    my ($trStDev, $mean) = trimmedStdDev(\@libSizes, 0.05);
    confess("Unable to calculate mean and stdev of trimmed lib sizes list") unless $mean;
    my $cutoff = 0;
    my $i = 2;
    while ( $cutoff <= 0 )
    {
        $cutoff = int($mean - $i * $trStDev + 0.5);
        $i -= 0.5;
    }
	return $cutoff;
}

=item singleRarefaction

Execute QIIME single_rarefaction.py

=cut

sub singleRarefaction
{
    my ($inFile, $outFile, $outFile2x2, $depth) = @_;
	_run("single_rarefaction.py -i $inFile -o $outFile -d $depth");
    filterOtuTable($outFile, $outFile2x2, 2, 2);
}

=item filterOtuTable

Filter OTUs that don't have an abundance of at least <threshold> at least <frequency> times.  For instance, keep OTUs that have at least 10 reads in 2 samples.  The input OTU table is expected to have the OTU ID in the first column and it's taxonomic classification in the last column.

=cut

sub filterOtuTable
{
    my ($inFile, $outFile, $threshold, $frequency) = @_;
    confess("Missing args") unless $inFile and $outFile and $threshold and $frequency;
    unless ( -e $inFile and -s $inFile )
    {
        print "WARNING: OTU table ($inFile) does not exist or is empty.\n";
        return;
    }
	open(my $in, '<', $inFile) or confess($!);
	open(my $out, '>', $outFile) or confess($!);
    my $line = <$in>;
    print $out $line;
	while( my $line = <$in> )
    {
        my @row = split(/\t/, $line);
        my $count = 0;
        for (my $i=1; $i<$#row; $i++)
        {
            if ( $row[$i] >= $threshold )
            {
                if ( ++$count >= $frequency )
                {
                    print $out $line;
                    last;
                }
            }

        }
	}
	close($in);
	close($out);
}

=item abundanceThreshold

From one OTU table, generate new OTU tables having OTUs above certain abundance threshold only. Abundance threshold calculated will be 10e-6, 10e-5, 10e-4, 10e-3, 10e-2 and 10e-1.

Then, generate a table having different alpha divesity values for each 
abundance threshold values.

=cut

sub abundanceThreshold
{
    my ($inFile, $outDir) = @_;
    confess("Missing args") unless $inFile and $outDir;

    # LOAD OTUS, CALC LIB SIZES
    my @otus = ();
    my @sizes = ();
    my $total = 0;
    open(my $in, '<', $inFile) or confess($!);
    my $hdr = <$in>;
    confess("Invalid OTU file") unless $hdr =~ /^#/;
    while (<$in>)
    {
        push @otus, $_;
        my @row = split(/\t/);
        my $sum = 0;
        for ( my $i=1; $i<$#row; $i++ ) { $sum += $row[$i] }
        push @sizes, $sum;
        $total += $sum;
    }
    close($in);

    # OUTPUT FILTERED OTU TABLES
    # keep only OTUs having abundance higher than xx%
    -d $outDir or mkdir($outDir) or confess($!);
    my @thresholds = (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1);
    my @otuFiles = map { "$outDir/otu.$_.tsv" } @thresholds;
    my @biomFiles = map { "$outDir/otu.$_.biom" } @thresholds;
    my @alphaFiles = map { "$outDir/otu.$_.alpha.tsv" } @thresholds;
    my @fhs;
    foreach my $otuFile (@otuFiles)
    {
        open(my $fh, '>', $otuFile) or confess($!);
        print $fh $hdr;
        push @fhs, $fh;
    }
    while (@otus)
    {
        my $otu = shift @otus;
        my $size = shift @sizes;
        my $fraction = $size / $total;
        for ( my $i=0; $i<=$#thresholds; $i++ )
        {
            next unless $fraction > $thresholds[$i];
            my $fh = $fhs[$i];
            print $fh $otu;
        }
    }
    foreach my $fh (@fhs) { close($fh) }
    for (my $i=0; $i<=$#thresholds; $i++)
    {
        ## Check if OTU table is empty.
        open(my $in, '<', $otuFiles[$i]) or confess($!);
        my $ok = -1;
        while (<$in>) { last if ++$ok }
        close($in);
        if ( $ok )
        {
            otuToBiom($otuFiles[$i], $biomFiles[$i]);
        } else
        {
            unlink($otuFiles[$i]);
            $otuFiles[$i] = $biomFiles[$i] = $alphaFiles[$i] = undef;
        }
    }

    # CALC ALPHA DIVERSITY
    my @threads;
    for ( my $i=0; $i<=$#thresholds; $i++ )
    {
        next unless defined($biomFiles[$i]);
        my ($thr) = new threads(\&alphaDiversity, $biomFiles[$i], $alphaFiles[$i]);
        push @threads, $thr;
    }
    foreach my $thr (@threads) { $thr->join }
}

=item alphaDiversity

Execute QIIME alpha_diversity.py and reformat header line.

=cut

sub alphaDiversity
{
    my ($inFile, $outFile) = @_;
    _run("alpha_diversity.py -i $inFile -o $outFile -m chao1,observed_species,shannon,simpson");
    # REFORMAT HEADER LINE
    open(my $in, '<', $outFile) or confess("Unable to open $outFile: $!");
    open(my $out, '>', "$outFile.tmp") or confess($!);
    print $out "#Sample\tchao1\tobserved_species\tshannon\tsimpson\n";
    my $head = <$in>;
    while (<$in>) { print $out $_ }
    close($in);
    close($out);
    move("$outFile.tmp", $outFile);
}

=item _set

Return value from config or default value.  If default not defined, it must be defined in the config file, otherwise logconfess.

=cut

sub _set
{
    my ($config, $part, $key, $default) = @_;
    confess("Missing args") unless $config and $part and $key;
    if ( !exists($config->{$part}) )
    {
        $config->{$part} = {};
        return $config->{$part}->{$key} = $default if defined($default);
        confess("Config file missing section, $part");
    } elsif ( exists($config->{$part}->{$key}) )
    {
        return $config->{$part}->{$key}
    } elsif ( defined($default) )
    {
        return $config->{$part}->{$key} = $default;
    } else
    {
        confess("Config missing required $part/$key");
    }
}

# NYI

=item calcMinDepths

If a fraction is provided, it will return the minimum OTU depth which will contain that fraction of reads.  Otherwise, it will return a list of minimum depths which capture 10%, 20%, ..., 90% of reads.

=cut

sub calcMinDepths
{
    my ($inFile, $pctAR) = @_;
    confess("Missing args") unless $inFile and @$pctAR;
    # CALC TOTAL SIZE
    open( my $in, '<', $inFile ) or confess("Can't open file, $inFile: $!");
    my $hdr = <$in>;
    my @sizes = ();
    my $totSize = 0;
    while (<$in>)
    {
        my @row = split(/\t/);
        my $size = 0;
        for ( my $i = 1; $i < $#row; $i++ ) { $size += $row[$i] // 0 }
        push @sizes, $size;
        $totSize += $size;
    }
    close($in);
    confess("Empty infile") unless $totSize;
    # CALC MIN DEPTHS
    my $sum = 0;
    my $minDepth;
    my @minDepths;
    my @x = map { $_/100 } @$pctAR;
    foreach my $x (@x)
    {
        my $target = $totSize * $x;
        while ($sum < $target) { $sum += $minDepth = shift @sizes }
        push @minDepths, $minDepth;
    }
    return @minDepths;
}

1;

=back

=head1 AUTHORS

Julien Tremblay, Edward Kirton

=head1 COPYRIGHT

Copyright (c) 2013 US DOE Joint Genome Institute.  Use freely under the same license as QIIME itself.

=cut

__END__
