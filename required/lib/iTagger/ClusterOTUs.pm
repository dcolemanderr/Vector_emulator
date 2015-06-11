=pod

=head1 NAME

iTagger::ClusterOTUs

=head1 DESCRIPTION

OTU clustering and taxonomic classification, wrapping usearch and RDP classifier.  To increase speed and allow processing of larger number of clusters than will fit in RAM (particularly for 32bit usearch), alternates between single-threaded clustering and multi-threaded searching. Requires 4-16 threads.

=head1 FUNCTIONS

=over 5

=cut

package iTagger::ClusterOTUs;

use strict;
use warnings;
use Carp;
use Env qw(TMPDIR RDP_JAR_PATH);
use constant { MAX_RADIUS => 3, MIN_ITERATION_SIZE => 10_000, MAX_ITERATION_SIZE => 150_000, MIN_LIB_SIZE => 1_000, MIN_SEQ_SIZE_TO_CLUSTER => 2, SEQS_PER_SINGLETS_FILE => 10_000 };
require Exporter;
use File::Copy;
use Parallel::ForkManager;
use POSIX qw/mkfifo/;
use Config::Tiny;
use iTagger::FastaDb;
use iTagger::Stats qw(sum);

our $VERSION = 1.0;
our @ISA = qw(Exporter);
our @EXPORT = qw(clusterOTUs);
our $tmpDir = -e '/global/scratch/colemand' ? '/global/scratch/colemand' : $TMPDIR;
our $usearchTime = 0;

=item clusterOTUs

This is the only exported function.  Performs iterative clustering, and generates OTU Fasta and abundance tables.

=cut

sub clusterOTUs
{
    my ($configFile, $dir, $threads) = @_;
    confess("Missing args") unless $configFile and $dir;
    my $start = time;

    my $outDir = "$dir/otu";
    -d $outDir or mkdir($outDir) or confess("Unable to mkdir $outDir: $!");
    my $logFile = "$outDir/cluster.log";

    my $config = Config::Tiny->read($configFile);
    confess("Config missing USEARCH section") unless exists($config->{USEARCH});
    my $maxRadius = _set($config, 'USEARCH', 'MAX_RADIUS', MAX_RADIUS);
    confess("MAX_RADIUS=$maxRadius invalid; must be 1..5") unless $maxRadius >= 1 and $maxRadius <= 5;
    my $refDb = _set($config, 'USEARCH', 'REF_DB', '');
    my $iterationSize = _set($config, 'USEARCH', 'ITERATION_SIZE', 0);
    my $minLibSize = _set($config, 'USEARCH', 'MIN_LIB_SIZE', MIN_LIB_SIZE);
    unless ($threads)
    {
        $threads = `nproc`;
        chomp $threads;
        confess("Number of threads not defined") unless $threads;
    }

    my $singletBase = "$tmpDir/$$.singlets";
    my $tmpOtuFile = "$tmpDir/$$.otu.fasta";
    my $tmpFile1   = "$tmpDir/$$.seqs.1.fasta";
    my $tmpFile2   = "$tmpDir/$$.seqs.2.fasta";
    my $tmpNoChim  = "$tmpDir/$$.otu.nochim.fasta";
    my ($libOrder, $counter);
    my $numClusterable;
    my $numSinglets;
    my $libClusterableSizes;
    my $libSingletSizes;

    for ( my $radius = 1; $radius <= $maxRadius; $radius++ )
    {
        if ( $radius == 1 )
        {
            ($libOrder, $iterationSize, $numClusterable, $numSinglets, $libClusterableSizes, $libSingletSizes) = _input($dir, $tmpFile1, $tmpFile2, $singletBase, $iterationSize, $logFile, $minLibSize, $threads);
        } else
        {
            $counter = _splitFasta($tmpOtuFile, $tmpFile1, $tmpFile2, $iterationSize);
        }
        while (1)
        {
            my $numClusters = _cluster($tmpFile1, $tmpOtuFile, $radius, $libOrder, $threads);
            last unless -e $tmpFile2 and -s $tmpFile2;
            $counter = _search($tmpOtuFile, $tmpFile1, $tmpFile2, $radius, $threads, $iterationSize);
        }
    }
    unlink($tmpFile2) if -e $tmpFile2;

    # REPORT: COUNT CHIMERA FILTERED DURING CLUSTERING
    my $beforeRefChimeraFilter = _countChimeraFilteredClustering($tmpOtuFile, $libClusterableSizes, $logFile);

    # MAP AND COUNT SINGLETS
    if ( -e "$singletBase.1.fasta" )
    {
        ($counter, $beforeRefChimeraFilter) = _searchSinglets($tmpOtuFile, $singletBase, $maxRadius, $threads, $libOrder, $libSingletSizes, $logFile);
    }

    # CHIMERA FILTER
    my $afterRefChimeraFilter = _chimeraFilter($tmpOtuFile, $refDb, $tmpNoChim, $beforeRefChimeraFilter, $logFile);
    move($tmpNoChim, $tmpOtuFile);
    move("$tmpNoChim.obs", "$tmpOtuFile.obs");

    # RDP CLASSIFIER
    my $fastaFile = "$outDir/otu.fasta.bz2";
    my $rdpFile = "$outDir/rdp.tsv";
    my $otuTaxFile = "$outDir/otu.tax.tsv";
    my $otuNoTaxFile = "$outDir/otu.unk.tsv";
    my $libClassified = _classify($config, $tmpOtuFile, $fastaFile, $rdpFile, $otuTaxFile, $otuNoTaxFile, $libOrder, $logFile, $afterRefChimeraFilter, $threads);
    unlink($tmpOtuFile);

    # TAXONOMY FILTER
    my $taxFilteredFile = "$outDir/otu.tax.filtered.tsv";
    my $numTaxPass = _taxFilter($config, $otuTaxFile, $taxFilteredFile, $logFile, $libClassified);

    # DONE
    print "total usearch time = $usearchTime sec\n";
    my $dur = time - $start;
    print "Clustered $numClusterable clusterable, $numSinglets singlets, yielding $numTaxPass good clusters after $dur sec\n";
}

## PRIVATE FUNCTIONS

=item _input

Combine and depreplicate seq-obs files from multiple libraries, output in FASTA format, sorted by decreasing abundance, with 'size=<int>' tag in header as required by USEARCH.

=cut

# TODO USE THREADS
sub _input
{
    my ($dir, $outFile1, $outFile2, $singletBase, $iterationSize, $logFile, $minLibSize, $threads) = @_;
    confess("Missing args") unless $dir and $outFile1 and $outFile2 and $logFile and $minLibSize and $threads;

    # FIND SEQOBS FILES, CHECK LIB SIZES -- SKIP IF TOO SMALL
    my $dataDir = "$dir/data";
    opendir(DIR, $dataDir) or confess("Unable to readdir $dataDir: $!");
    my @contents = sort { $a cmp $b } grep { $_ !~ /^\./ } readdir(DIR);
    closedir(DIR);
    my @libs;
    my %libFiles;
    my %libSizes; # wrote to log below
    my @inFiles = ();
    foreach my $lib (@contents)
    {
        unless ( -d "$dataDir/$lib" )
        {
            carp("Folder does not exist: $dataDir/$lib");
            next;
        }
        my $file = "$dataDir/$lib/seqobs.tsv.bz2";
        unless ( -e $file )
        {
            carp("Infile not found: $dataDir/$lib/seqobs.tsv.bz2");
            next;
        }
        open(my $fh, "bzcat $file |") or confess("Unable to read $file: $!");
        my $hdr = <$fh>; # header
        chomp $hdr;
        confess("Invalid seqobs file header") unless $hdr =~ /^#lib=(.+);size=(\d+)$/;
        confess("Invalid seqobs file library id") unless $1 eq $lib;
        my $libSize = $2;
        if ( $libSize < $minLibSize )
        {
            print "Discard lib, $lib, as total size, $libSize < $minLibSize\n";
            close($fh);
            next;
        }
        push @libs, $lib;
        $libFiles{$lib} = $file;
        $libSizes{$lib} = $libSize;
        # CACHE FIRST RECORD
        my $line = <$fh>;
        chomp $line;
        my @row = split(/\t/, $line);
        confess("Invalid file: $file") unless @row == 2;
        my ($seq, $count) = @row;
        push @inFiles, [ $fh, $#libs, $seq, $count ];
    }
    confess("No seqobs files found") unless scalar(keys %libFiles) and scalar(@inFiles);
    my @obs0 = ( (0) x scalar(@libs) ); # empty obs list

    # DEREPLICATE, FILTER LOW-ABUNDANCE CLUSTERS, SORT.
    # CLUSTERABLE SEQUENCES WRITTEN TO .fasta AND .fasta.obs FILES,
    # WHILE LOW-ABUNDANCE SEQUENCES (E.G. SINGLETS) DO NOT HAVE .fasta.obs FILES;
    # INSTEAD A SPARSE REPRESENTATION IS USED IN THE HEADER.
    my $tmpFile = "$tmpDir/$$.seqobs.tmp";
    open(my $out, "| sort -k1,1nr -T $tmpDir -o $tmpFile") or confess("Unable to pipe sort $tmpFile: $!");
    my $counter = 0; # number of clusterable sequences (not reads)
    my $singletCounter = 0; # unique id for singlets
    my $singletFileNum = 1; # unique id for singlets filenames
    my $singletFile = "$singletBase.$singletFileNum.fasta";
    my $singletFileCounter = 0; # number of seqs in current singlets file
    open(my $low, '>', $singletFile) or confess("Unable to write $singletFile: $!");
    my ($i, $seq, $count) = _nextSeq(\@inFiles);
    my @obs = @obs0; # deep copy
    $obs[$i] = $count;
    my @libClusterableSizes;
    my @libSingletSizes;
    while (@inFiles)
    {
        my ($i2, $seq2, $count2) = _nextSeq(\@inFiles);
        if ( $seq eq $seq2 )
        {
            $obs[$i2] += $count2;
        } else
        {
            my $size = _size(\@obs);
            if ( $size < MIN_SEQ_SIZE_TO_CLUSTER )
            {
                if ( ++$singletFileCounter > SEQS_PER_SINGLETS_FILE )
                {
                    close($low);
                    $singletFileCounter = 0;
                    ++$singletFileNum;
                    $singletFile = "$singletBase.$singletFileNum.fasta";
                    open($low, '>', $singletFile) or confess("Unable to write $singletFile: $!");
                }
                print $low '>', ++$singletCounter, ';size=', $size, ';obs=', _obsToSparseString(\@obs), "\n" ,$seq, "\n";
                for ( my $i=0; $i<=$#obs; $i++ ) { $libSingletSizes[$i] += $obs[$i] }
            } else
            {
                print $out $size,"\t",$seq,"\t",join(',', @obs),"\n";
                ++$counter;
                for ( my $i=0; $i<=$#obs; $i++ ) { $libClusterableSizes[$i] += $obs[$i] }
            }
            $seq = $seq2;
            @obs = @obs0; # deep copy empty list
            $obs[$i2] = $count2;
        }
    }
    my $size = _size(\@obs);
    if ( $size < MIN_SEQ_SIZE_TO_CLUSTER )
    {
        print $low '>', ++$singletCounter, ';size=', $size, ';obs=', _obsToSparseString(\@obs), "\n" ,$seq, "\n";
        for ( my $i=0; $i<=$#obs; $i++ ) { $libSingletSizes[$i] += $obs[$i] }
    } else
    {
        print $out $size,"\t",$seq,"\t",join(',', @obs),"\n";
        ++$counter;
        for ( my $i=0; $i<=$#obs; $i++ )
        {
            $libClusterableSizes[$i] += $obs[$i];
        }
    }
    close($out);
    close($low);

    # SUMMARY STATS LOG
    my @libSizes = map { $libSizes{$_} } @libs;
    my $totLibSize = sum(\@libSizes);
    my @clusterablePct;
    my @singletsPct;
    my $totClusterable = sum(\@libClusterableSizes);
    my $totSinglet = sum(\@libSingletSizes);
    my $totClusterablePct = $totLibSize ? int($totClusterable/$totLibSize*1000+0.5)/10 . '%' : 'N/A';
    my $totSingletPct = $totLibSize ? int($totSinglet/$totLibSize*1000+0.5)/10 . '%' : 'N/A';
    for ( my $i=0; $i<=$#libSizes; $i++ )
    {
        $clusterablePct[$i] = $libSizes[$i] ? (int($libClusterableSizes[$i]/$libSizes[$i]*1000+0.5)/10) . '%' : 'N/A';
        $singletsPct[$i] = $libSizes[$i] ? (int($libSingletSizes[$i]/$libSizes[$i]*1000+0.5)/10) . '%' : 'N/A';
    }
    open(my $log, '>', $logFile) or confess("Unable to write $logFile: $!");
    print $log
        join("\t", 'LIB', 'TOTAL', map { _commify($_) } @libs),"\n",
        join("\t", 'INPUT', _commify($totLibSize), map { _commify($_) } @libSizes), "\n",
        join("\t", 'CLUSTERABLE', _commify($totClusterable), map { _commify($_) } @libClusterableSizes), "\n",
        join("\t", 'CLUSTERABLE PCT', $totClusterablePct, @clusterablePct), "\n",
        join("\t", 'SINGLETS', _commify($totSinglet), map { _commify($_) } @libSingletSizes), "\n",
        join("\t", 'SINGLETS PCT', $totSingletPct, @singletsPct), "\n";
    close($log);

    # SPLIT
    $iterationSize = int($counter/5) unless $iterationSize;
    if ( $iterationSize < MIN_ITERATION_SIZE ) { $iterationSize = MIN_ITERATION_SIZE }
    elsif ( $iterationSize > MAX_ITERATION_SIZE ) { $iterationSize = MAX_ITERATION_SIZE }
    _initialSplitFasta($tmpFile, $outFile1, $outFile2, $iterationSize);
    unlink($tmpFile);

    return (\@libs, $iterationSize, $counter, $singletCounter, \@libClusterableSizes, \@libSingletSizes);
}

=item _nextSeq

Private function returns next sorted seq, used for merging input.

=cut

sub _nextSeq
{
    my ($inFiles) = @_; # fh, index, seq, count
    # PICK NEXT
    my $i=0;
    for (my $j=1; $j<=$#$inFiles; $j++)
    {
        $i=$j if $inFiles->[$j]->[2] lt $inFiles->[$i]->[2];
    }
    my ($fh, $index, $seq, $count) = @{$inFiles->[$i]};
    # BUFFER NEXT
    my $line = <$fh>;
    if ( $line )
    {
        chomp $line;
        my ($seq, $count) = split(/\t/, $line);
        $inFiles->[$i]->[2] = $seq;
        $inFiles->[$i]->[3] = $count;
    } else
    {
        splice($inFiles, $i, 1); # EOF
    }
    return ($index, $seq, $count);
}

=item _size

Private function returns total size of @obs list

=cut

sub _size
{
    my ($obs) = @_;
    my $size = 0;
    foreach (@$obs) { $size += $_ }
    return $size;
}

=item _obsToSparseString

Given a list of observations (counts), return a string with a sparse array reprentation.

=cut

sub _obsToSparseString
{
    my ($obs) = @_;
    my @items = ();
    for ( my $i=0; $i<=$#$obs; $i++ )
    {
        my $x = $obs->[$i];
        push @items, "$i:$x" if $x;
    }
    return join(',', @items);
}

=item _initialSplitFasta

Private function to create a file of specified size (for clustering) and write the remainder to a second file (for searching).

=cut

sub _initialSplitFasta 
{
    my ($inFile, $outFile1, $outFile2, $iterationSize) = @_;
    open(my $in, '<', $inFile) or confess("Unable to read $inFile: $!");
    open(my $out, '>', $outFile1) or confess("Unable to write $outFile1: $!");
    open(my $obsOut, '>', "$outFile1.obs") or confess("Unable to write $outFile1.obs: $!");
    my $index = -1;
    while (<$in>)
    {
        chomp;
        my ($size, $seq, $obs) = split(/\t/);
        if ( ++$index == $iterationSize )
        {
            close($out);
            close($obsOut);
            open($out, '>', $outFile2) or confess("Unable to write $outFile2: $!");
            open($obsOut, '>', "$outFile2.obs") or confess("Unable to write $outFile2.obs: $!");
        }
        print $out '>',$index,';size=',$size,"\n",$seq,"\n";
        print $obsOut $index, "\t", $obs, "\n";
    }
    close($out);
    close($obsOut);
    return $index+1;
}

=item _cluster

Given Fasta file of dereplicated sequences (sorted by decreasing abundance), cluster OTUs at specified radius, map the reads to the centroids, calculate OTU sizes, and output OTU Fasta.  Delete input file.

=cut

sub _cluster
{
    my ($inFile, $outFile, $radius, $libOrder, $threads) = @_;
    confess("Missing args") unless $inFile and $outFile and $radius and $threads;
    my $identity = (100-$radius)/100;
    my @tmpFiles = ();

    # CLUSTER
    push @tmpFiles, my $centroids = "$tmpDir/$$.centroids.fasta";
    my $start = time;
    _run("usearch -cluster_otus $inFile -otu_radius_pct $radius -otus $centroids");
    $usearchTime += (time-$start);

    # MAP CLUSTERED READS TO CENTROIDS
    push @tmpFiles, my $hitFile = "$tmpDir/$$.centroids.hits";
    $start = time;
    _run("usearch -threads $threads -usearch_global $inFile -db $centroids -strand plus -id $identity -userout $hitFile -userfields query+target");
    $usearchTime += (time-$start);

    # LOAD CENTROIDS
    my $otus = _loadOtus($centroids, "$inFile.obs");

    # LOAD HITS
    my %hits = ();
    open(my $in, '<', $hitFile) or confess("Unable to read $hitFile: $!");
    while (<$in>)
    {
        chomp;
        my ($query, $subj) = split(/\t/);
        next if $query eq $subj;
        confess("Invalid query ID") unless $query =~ /^(\d+);size=\d+$/;
        my $queryId = $1;
        confess("Invalid subject ID") unless $subj =~ /^(\d+);size=\d+$/;
        my $subjId = $1;
        next if $queryId < $subjId; # don't add to smaller cluster
        confess("Invalid OTU $subjId") unless exists($otus->{$subjId});
        carp("Clustering error, OTU $queryId should not exist; it should have been added to $subjId") if exists($otus->{$queryId}); # TODO
        $hits{$queryId} = $subjId;
    }
    close($in);

    # GET QUERY OBS, ADD TO OTU
    open($in, '<', "$inFile.obs") or confess("Unable to read $inFile.obs: $!");
    while (<$in>)
    {
        chomp;
        my ($queryId, $obs) = split(/\t/);
        next unless exists($hits{$queryId});
        my $subjId = $hits{$queryId};
        my @queryObs = split(/,/, $obs);
        my $subjObs = $otus->{$subjId}->[1];
        for (my $i=0; $i<=$#$subjObs; $i++) { $subjObs->[$i] += $queryObs[$i] }
    }
    close($in);

    # CLEANUP
    unlink(@tmpFiles, $inFile, "$inFile.obs");

    # OUTPUT
    return _outputOtus($otus, $outFile);
}

=item _splitFasta

Private function to write specified number of seqs to one file and remainder to another.

=cut

sub _splitFasta
{
    my ($inFile, $outFile1, $outFile2, $iterationSize) = @_;
    unlink($outFile2, "$outFile2.obs") if -e $outFile2;
    my $in = new iTagger::FastaDb($inFile);
    open(my $obsIn, '<', "$inFile.obs") or confess("Unable to read $inFile.obs: $!");
    open(my $out, '>', $outFile1) or confess("Unable to write $outFile1: $!");
    open(my $obsOut, '>', "$outFile1.obs") or confess("Unable to write $outFile1.obs: $!");
    my $counter = 0;
    while ( my $rec = $in->next_seq )
    {
        print $out $rec->output;
        my $obs = <$obsIn>;
        print $obsOut $obs;
        if ( ++$counter == $iterationSize )
        {
            close($out);
            close($obsOut);
            open($out, '>', $outFile2) or confess("Unable to write $outFile2: $!");
            open($obsOut, '>', "$outFile2.obs") or confess("Unbale to write $outFile2.obs: $!");
        }
    }
    close($out);
    close($obsOut);
    return $counter;
}

=item _outputOtus

Order OTUs by decreasing size and output.  OTUs are renamed, starting with 0.

=cut

sub _outputOtus
{
    my ($otus, $fastaFile) = @_;
    my $obsFile = "$fastaFile.obs";

    # CALCULATE TOTAL SIZES
    foreach my $id ( keys %$otus )
    {
        my $size = 0;
        my $obs = $otus->{$id}->[1];
        foreach (@$obs) { $size += $_ }
        $otus->{$id}->[2] = $size;
    }

    # SORT BY DECREASING SIZE
    my @ids = sort { $otus->{$b}->[2] <=> $otus->{$a}->[2] } keys %$otus;

    # OUTPUT
    open(my $out1, '>', $fastaFile) or confess("Unable to write $fastaFile: $!");
    open(my $out2, '>', $obsFile) or confess("Unable to write $obsFile: $!");
    my $index = -1; # NEW ID
    while (@ids)
    {
        my $id = shift @ids;
        my ($seq, $obs, $size) = @{$otus->{$id}};
        print $out1 '>', ++$index, ';size=', $size, "\n", $seq, "\n";
        print $out2 $index, "\t", join(',', @$obs), "\n";
    }   
    close($out1);
    close($out2);
    return $index+1; # NUMBER OF OTUS
}

=item _loadOtus

Private function to load OTUs from disk into RAM.  The obs file may contain more records than the Fasta file; only the OTUs in the Fasta file are used.  The Fasta file may contain OTUs out of order.

=cut

sub _loadOtus
{
    my ($fastaFile, $obsFile) = @_;
    my %otus = (); # id => [seq, obs]
    # LOAD OTUS FROM FASTA FILE
    my $db = new iTagger::FastaDb($fastaFile);
    while ( my $rec = $db->next_seq )
    {
        my $hdr = $rec->id;
#	print "\n$hdr\n"; die;
#modified 05292015
#       confess("Invalid OTU ID: $hdr") unless $hdr =~ /^(\d+);size=\d+$/;
       	confess("Invalid OTU ID: $hdr") unless $hdr =~ /^(\d+)/;

        my $id = $1;
        $otus{$id} = [ $rec->seq ];
    }
    # GET OBS FOR OTUS FOUND ABOVE ONLY (SOME MAY HAVE BEEN FILTERED AS CHIMERIC)
    open(my $in, '<', $obsFile) or confess("Unable to read $obsFile: $!");
    while (<$in>)
    {
        chomp;
        my ($id, $obs) = split(/\t/);
        next unless exists($otus{$id});
        my @obs = split(/,/, $obs);
        $otus{$id}->[1] = \@obs;
    }
    close($in);
    return \%otus;
}

=item _sumOtus

Count the total number of sequence per libs in all OTUs.

=cut

sub _sumOtus
{
    my ($inFile) = @_;
    my @sums;
    open (my $in, '<', $inFile) or confess("Unable to read $inFile: $!");
    while (<$in>)
    {
        chomp;
        my ($id, $obs) = split(/\t/);
        my @obs = split(/,/, $obs);
        for ( my $i=0; $i<=$#obs; $i++ ) { $sums[$i] += $obs[$i] }
    }
    return \@sums;
}

=item _search

Map unclustered reads files to OTUs, calculate new OTU sizes, sort by decreasing abundance, and generate new file to be clustered (with OTUs at top of file), plus second file if number of unclustered sequences exceeds maximum number of seqs per iteration.

=cut

# SEARCH FILE2 VS OTU, DELETE OTU, AND WRITE NEW FILE1 AND FILE2, WHERE NEW OTUS ARE AT TOP OF FILE1.  FILE2 MAY BE EMPTY IF THERE ARE LESS SEQUENCES TO CLUSTER THAN ITERATION SIZE.
sub _search
{
    my ($otuFile, $file1, $file2, $radius, $threads, $iterationSize) = @_;
    confess("Missing args") unless $otuFile and $file1 and $file2 and $radius and $threads and $iterationSize;
    my $identity = (100-$radius)/100;

    # MAP FILE2 VS OTU (CENTROIDS) FILE
    my $hitFile = "$tmpDir/$$.usearch.hits";
    my $start = time;
    _run("usearch -threads $threads -usearch_global $file2 -db $otuFile -strand plus -id $identity -userout $hitFile -userfields query+target");
    $usearchTime += (time-$start);

    # LOAD OTUS
    my $otus = _loadOtus($otuFile, "$otuFile.obs");

    # LOAD HITS
    my %hits = ();
    open(my $in, '<', $hitFile) or confess("Unable to read $hitFile: $!");
    while (<$in>)
    {
        chomp;
        my ($query, $subj) = split(/\t/);
        next if $query eq $subj;
        confess("Invalid query ID: $query") unless $query =~ /^(\d+);size=\d+$/;
        my $queryId = $1;
        confess("Invalid subject ID: $subj") unless $subj =~ /^(\d+);size=\d+$/;
        my $subjId = $1;
        confess("OTU not found: $subjId") unless exists($otus->{$subjId});
        $hits{$queryId} = $subjId;
    }
    close($in);

    # GET QUERY OBS, ADD TO OTU
    open($in, '<', "$file2.obs") or confess("Unable to read $file2.obs: $!");
    while (<$in>)
    {
        chomp;
        my ($queryId, $obs) = split(/\t/);
        next unless exists($hits{$queryId});
        my $subjId = $hits{$queryId};
        my @queryObs = split(/,/, $obs);
        my $subjObs = $otus->{$subjId}->[1];
        for (my $i=0; $i<=$#$subjObs; $i++) { $subjObs->[$i] += $queryObs[$i] }
    }
    close($in);

    # OUTPUT UNCLUSTERED SEQS
    my $counter = scalar(keys %$otus); # iteration size includes otus
    my $time = time;
    my $tmpFile1 = "$tmpDir/$$.$time.tmpFile.1.fasta";
    my $tmpFile2 = "$tmpDir/$$.$time.tmpFile.2.fasta";
    my $in1 = new iTagger::FastaDb($file2);
    open(my $in2, '<', "$file2.obs") or confess("Unable to read $file2.obs: $!");
    open(my $out1, '>', $tmpFile1) or confess("Unable to write $tmpFile1: $!");
    open(my $out2, '>', "$tmpFile1.obs") or confess("Unable to write $tmpFile1.obs: $!");
    while ( my $rec1 = $in1->next_seq )
    {
        confess("Invalid header") unless $rec1->id =~ /^(\d+);size=\d+$/;
        my $id1 = $1;
        my $rec2 = <$in2>;
        my ($id2, $obs) = split(/\t/, $rec2);
        chomp $obs;
        confess("Fasta:obs mismatch ($id1:$id2)") unless $id1 == $id2;
        next if exists($hits{$id1});
        print $out1 $rec1->output;
        print $out2 $rec2;
        if ( ++$counter == $iterationSize )
        {
            close($out1);
            close($out2);
            open($out1, '>', $tmpFile2) or confess("Unable to write $tmpFile2: $!");
            open($out2, '>', "$tmpFile2.obs") or confess("Unable to write $tmpFile2.obs: $!");
        }
    }
    close($in2);
    close($out1);
    close($out2);
    unlink($file2, "$file2.obs", $hitFile);
    if ( -e $tmpFile2 )
    {
        move($tmpFile2, $file2) or confess("Unable to mv $tmpFile2 $file2: $!");
        move("$tmpFile2.obs", "$file2.obs") or confess("Unable to mv $tmpFile2.obs $file2.obs: $!");
    }

    # OUTPUT OTUS AND APPEND UNMAPPED PART 1
    _outputOtus($otus, $file1);
    _run("cat $tmpFile1 >> $file1");
    _run("cat $tmpFile1.obs >> $file1.obs");
    unlink($tmpFile1, "$tmpFile1.obs");
    return $counter;
}

=item _countChimeraFilteredClustering

Count number of sequences per library that were filtered as chimeric during clustering.

=cut

sub _countChimeraFilteredClustering
{
    my ($tmpOtuFile, $libClusterableSizes, $logFile) = @_;
    my @sizes = ( (0) x scalar(@$libClusterableSizes) );
    my $otus = _loadOtus($tmpOtuFile, "$tmpOtuFile.obs");
    foreach my $id ( keys %$otus )
    {
        my $obs = $otus->{$id}->[1];
        for ( my $i=0; $i<=$#sizes; $i++ ) { $sizes[$i] += $obs->[$i] }
    }
    my @chimeric;
    for (my $i=0; $i<=$#sizes; $i++)
    {
        $chimeric[$i] = $libClusterableSizes->[$i] - $sizes[$i];
    }
    my $totLibClusterableSize = sum($libClusterableSizes);
    my $totChimPct = $totLibClusterableSize ? (int(sum(\@chimeric)/$totLibClusterableSize*1000+0.5)/10) . '%' : 'N/A';
    my @chimericPct;
    for (my $i=0; $i<=$#chimeric; $i++)
    {
        $chimericPct[$i] = $libClusterableSizes->[$i] ? (int($chimeric[$i]/$libClusterableSizes->[$i]*1000+0.5)/10) . '%' : 'N/A';
    }
    open(my $log, '>>', $logFile) or confess("Unable to append $logFile: $!");
    print $log 
        join("\t", 'CLUST CHIMERA', _commify(sum(\@chimeric)), map { _commify($_) } @chimeric), "\n",
        join("\t", 'CLUST CHIMERA PCT', $totChimPct, @chimericPct), "\n";
    close($log);
    return \@sizes;
}

=item _searchSinglets

Map unclustered reads files to OTUs, calculate new OTU sizes, sort by decreasing abundance, and overwrite OTU file.  Deletes singlet files.

=cut

sub _searchSinglets
{
    my ($otuFile, $singletBase, $radius, $threads, $libs, $libSingletSizes, $logFile) = @_;
    confess("Missing args") unless $otuFile and $singletBase and $radius and $threads and $libs and $libSingletSizes and $logFile;
    my $identity = (100-$radius)/100;

    # LOAD OTUS
    my $otus = _loadOtus($otuFile, "$otuFile.obs");

    # MAP SINGLETS (IN UNSPECIFIED NUMBER OF FILES)
    my $fileNum = 0;
    my $numLibs = scalar(@$libs);
    my @libSingletsClustered = ( (0) x $numLibs ); # FOR SUMMARY STATS
    while (1)
    {
        ++$fileNum;
        my $singletFile = "$singletBase.$fileNum.fasta";
        last unless -e $singletFile;

        # MAP SINGLETS VS OTU (CENTROIDS) FILE
        my $hitFile = "$tmpDir/$$.singlets.hits";
        my $start = time;
        _run("usearch -threads $threads -usearch_global $singletFile -db $otuFile -strand plus -id $identity -userout $hitFile -userfields query+target");
        $usearchTime += (time-$start);

        # PARSE HITS, ADD TO OTUS
        open(my $in, '<', $hitFile) or confess("Unable to read $hitFile: $!");
        while (<$in>)
        {
            chomp;
            my ($query, $subj) = split(/\t/);
            # QUERY
            confess("Invalid query ID: $query") unless $query =~ /^(\d+);size=\d+;obs=(.+)$/;
            my $queryId = $1;
            my $queryObs = $2;
            my @queryObs = split(/,/, $queryObs);
            # TARGET/SUBJECT
            confess("Invalid subject ID: $subj") unless $subj =~ /^(\d+);size=\d+$/;
            my $subjId = $1;
            confess("OTU not found: $subjId") unless exists($otus->{$subjId});
            my $subjObs = $otus->{$subjId}->[1];
            # ADD
            foreach my $item (@queryObs)
            {
                my ($i, $x) = split(/:/, $item);
                $subjObs->[$i] += $x;
                $libSingletsClustered[$i] += $x;
            }
        }
        close($in);
        unlink($singletFile, $hitFile);
    }

    # COUNT OTU SIZES
    my @beforeRefChimeraFilter = ( (0) x scalar(@$libs) );
    foreach my $id (keys %$otus)
    {
        my $obs = $otus->{$id}->[1];
        for (my $i=0; $i<=$#$libs; $i++ ) { $beforeRefChimeraFilter[$i] += $obs->[$i] }
    }

    # LOG
    my @libSingletsClusteredPct;
    my $totSingClus = my $totSingSize = 0;
    for (my $i=0; $i<=$#libSingletsClustered; $i++)
    {
        if ( $libSingletSizes->[$i] )
        {
            $libSingletsClusteredPct[$i] = $libSingletSizes->[$i] ?  (int($libSingletsClustered[$i] / $libSingletSizes->[$i] * 1000 + 0.5)/10) . '%' : 'N/A';
            $totSingClus += $libSingletsClustered[$i];
            $totSingSize += $libSingletSizes->[$i];
        } else
        {
            $libSingletsClusteredPct[$i] = 'N/A';
        }
    }
    my $totSingMapPct = $totSingSize ? (int($totSingClus/$totSingSize*1000+0.5)/10) . '%' : 'N/A';
    open(my $log, '>>', $logFile) or confess("Unable to append $logFile: $!");
    print $log
        join("\t", 'SINGLETS MAPPED', _commify(sum(\@libSingletsClustered)), map { _commify($_) } @libSingletsClustered), "\n",
        join("\t", 'SINGLETS MAPPED PCT', $totSingMapPct, @libSingletsClusteredPct), "\n",
        join("\t", 'BEFORE REF-CHIMERA FILTER', _commify(sum(\@beforeRefChimeraFilter)), map { _commify($_) } @beforeRefChimeraFilter), "\n";
    close($log);

    my $counter = _outputOtus($otus, $otuFile);
    return ($counter, \@beforeRefChimeraFilter);
}

=item _chimeraFilter

Filter chimeric sequences by comparison to reference database.

=cut

sub _chimeraFilter
{
    my ($inFile, $refDb, $outFile, $beforeSizes, $logFile) = @_;
    confess("Missing args") unless $inFile and $outFile;

    my $afterSizes = $beforeSizes;
    if ( $refDb )
    {
        # FILTER CHIMERAS
        my $start = time;
        _run("usearch -uchime_ref $inFile -db $refDb -strand plus -nonchimeras $outFile");
        $usearchTime += (time-$start);

        # OUTPUT NEW OTUS FASTA AND OBS FILES
        my $otus = _loadOtus($outFile, "$inFile.obs");
        _outputOtus($otus, $outFile);

        # SUMMARY STATS
        $afterSizes = _sumOtus("$outFile.obs");
        my @chimeric;
        my @chimericPct;
        for ( my $i=0; $i<=$#$afterSizes; $i++ )
        {
            my $before = $beforeSizes->[$i];
            $chimeric[$i] = $before - $afterSizes->[$i];
            $chimericPct[$i] = $before ? int($chimeric[$i]/$before*1000+0.5)/10 . '%' : 'N/A';
        }
        my $totBefore = sum($beforeSizes);
        my $totChimeric = $totBefore - sum($afterSizes);
        my $totChimericPct = $totBefore ? int($totChimeric/$totBefore*1000+0.5)/10 . '%' : 'N/A';
        open(my $log, '>>', $logFile) or confess("Unable to append $logFile: $!");
        print $log
            join("\t", 'REF-CHIMERA', _commify($totChimeric), map { _commify($_) } @chimeric), "\n",
            join("\t", 'REF-CHIMERA PCT', $totChimericPct, @chimericPct), "\n";
        close($log);
    }

    # SUMMARY LOG
    my @afterSizesPct = ();
    for ( my $i=0; $i<=$#$afterSizes; $i++ )
    {
        $afterSizesPct[$i] = $beforeSizes ? int($afterSizes->[$i]/$beforeSizes->[$i]*1000 + 0.5)/10 . '%' : 'N/A';
    }
    my $totBeforeSizes = sum($beforeSizes);
    my $totAfterSizes = sum($afterSizes);
    my $totAfterSizesPct = $totBeforeSizes ? int($totAfterSizes/$totBeforeSizes*1000+0.5)/10 . '%' : 'N/A';
    open(my $log, '>>', $logFile) or confess("Unable to append $logFile: $!");
    print $log
        join("\t", 'CLUSTERED', _commify($totAfterSizes), map { _commify($_) } @$afterSizes), "\n",
        join("\t", 'CLUSTERED PCT', $totAfterSizesPct, @afterSizesPct), "\n";
    close($log);
    return $afterSizes;
}

=item classify

Run RDP Classifier in multiple threads and generate environmental community table.

=cut

sub _classify
{
    my ($config, $inFile, $outFasta, $rdpOutFile, $otuPassFile, $otuFailFile, $libNames, $logFile, $libClustered, $threads) = @_;
    confess("Missing args") unless $config and $inFile and $outFasta and $rdpOutFile and $otuPassFile and $otuFailFile and $libNames and $logFile and $threads;
    confess("\$RDP_JAR_PATH not defined") unless $RDP_JAR_PATH;
    confess("\$TMPDIR not defined") unless $TMPDIR;
    my $start = time;

    # CONFIG
    my $trainingFile = _set($config, 'RDP_CLASSIFIER', 'TRAINING_FILE');
    my $minWords = _set($config, 'RDP_CLASSIFIER', 'MIN_WORDS', 120);
    my $level = _set($config, 'RDP_CLASSIFIER', 'LEVEL', 'class');
    my $cutoff = _set($config, 'RDP_CLASSIFIER', 'CUTOFF', 0.5);
    confess("Invalid cutoff, $cutoff") if $cutoff < 0 or $cutoff > 1;
    $cutoff = 0 unless defined($cutoff);
    $minWords = 120 unless $minWords;
    my %levelIndex = 
    (
        'kingdom' => 0, 'phylum'  => 1, 'class'   => 2, 
        'order'   => 3, 'family'  => 4, 'genus'   => 5
    );
    $level = lc($level);
    confess("Invalid level $level") unless exists($levelIndex{$level});

    # SPLIT INPUT TO TMPFILES, ROUND-ROBIN STYLE.  FIFOS ARE NOT USED FOR INPUT BECAUSE RDP CLASSIFIER WON'T READ FROM A FIFO.
    my @inFiles = map { "$tmpDir/$$.rdp.in.$_.fasta" } (1..$threads);
    my @fhs;
    for ( my $i=0; $i<$threads; $i++)
    {
        open( my $fh, '>', $inFiles[$i] ) or confess("Unable to write $inFiles[$i]: $!");
        push @fhs, $fh;
    }
    open( my $out, "| bzip2 --best > $outFasta") or confess("Unable to pipe gzip $outFasta: $!");
    my $db = new iTagger::FastaDb($inFile);
    my $i = -1;
    while ( my $rec = $db->next_seq )
    {
        # remove ";size=X" from header
        my $id = $rec->id;
        confess("Invalid id: $id") unless $id =~ /^(\d+);size=\d+/;
        $rec->base($id = $1);
        print $out $rec->output; # save copy
        $i = 0 if ++$i == $threads;
        my $fh = $fhs[$i];
        print $fh $rec->output; # rdp thread
    }
    while ( my $fh = shift @fhs ) { close($fh) }
    close($out);

    # RUN RDP CLASSIFIER. FIFOS AREN'T USED FOR OUTPUT BECAUSE NOT ALL QUERIES PRODUCE ANY OUTPUT FROM RDP CLASSIFIER.
    my @outFiles = map { "$tmpDir/$$.rdp.out.$_.tsv" } (1..$threads);
    my $pm = new Parallel::ForkManager($threads);
    for ( my $i=0; $i<$threads; $i++)
    {
        my $pid = $pm->start and next;
        _runRdp($RDP_JAR_PATH, $inFiles[$i], $outFiles[$i], $trainingFile, $minWords);
        $pm->finish;
    }
    $pm->wait_all_children;
    unlink(@inFiles);

    # GATHER OUTPUT
    my @libPass = ( (0) x scalar(@$libNames) ); # FOR STATS
    my @libFail = ( (0) x scalar(@$libNames) ); # FOR STATS
    my @buffer = ();
    for ( my $i=0; $i<$threads; $i++)
    {
        open( my $fh, '<', $outFiles[$i] ) or confess("Unable to read $outFiles[$i]: $!");
        my $line = <$fh>;
        if ( defined($line) )
        {
            chomp $line;
            my @row = split(/\t/, $line);
            my $id = shift @row;
            $fhs[$i] = $fh;
            $buffer[$i] = [ $id, \@row ];
        }
    }
    open(my $in, '<', "$inFile.obs") or confess("Unable to read $inFile.obs: $!");
    open($out, '>', $rdpOutFile) or confess("Unable to write $rdpOutFile: $!"); # save a copy
    open(my $pass, '>', $otuPassFile) or confess("Unable to write $otuPassFile: $!");
    open(my $fail, '>', $otuFailFile) or confess("Unable to write $otuFailFile: $!");
    print $pass join("\t", '#OTU', @$libNames, 'Consensus lineage'), "\n";
    print $fail join("\t", '#OTU', @$libNames, 'Consensus lineage'), "\n";
    while ( my $queryLine = <$in> )
    {
        chomp $queryLine;
        my ($id, $obs) = split(/\t/, $queryLine);
        my @obs = split(/,/, $obs);

        # next RDP result (some queries don't produce results)
        my $index = undef; # which thread index contains record for OTU $id
        my $found = 0;
        my $row = undef; # RDP result arrayref
        for ( my $i=0; $i < $threads; $i++ )
        {
            next unless defined($buffer[$i]); # this thread has no more results
            next unless $buffer[$i]->[0] == $id;
            $index = $i; # query OTU found at this index
            $row = $buffer[$index]->[1];
            $found = 1;
            last;
        }
        unless ($found)
        {
            # NO RDP RESULT FOR THIS QUERY
            print $fail join("\t", $id, @obs, ''), "\n";
            for ( my $i=0; $i<=$#obs; $i++ ) { $libFail[$i] += $obs[$i] }
            next;
        }

        # UPDATE BUFFER
        my $fh = $fhs[$index];
        my $line = <$fh>;
        if ( $line )
        {
            # update buffer with next RDP result
            chomp $line;
            my @nextRow = split(/\t/, $line);
            my $nextId = shift @nextRow;
            $buffer[$index] = [ $nextId, \@nextRow ];
        } else
        {
            # no more results for this thread
            close($fh);
            $buffer[$index] = undef;
        }

        # PRINT RDP OUTPUT (merge all threads' output)
        print $out join("\t", $id, @$row), "\n";

        # EVALUATE RDP RESULT
        my @taxon = map { $_ // 0 } @$row[4,7,10,13,16,19];
        my @prob  = map { $_ // 0 } @$row[6,9,12,15,18,21];
        my $ok = $prob[$levelIndex{$level}] >= $cutoff ? 1 : 0;
        my @best = ();
        for (my $i=0; $i<=$#taxon; $i++)
        {
            if ($prob[$i] >= $cutoff) { push @best, $taxon[$i] }
            else { last }
        }
        @taxon=@best;
        my $taxon = join(';', @taxon);
        if ( $ok )
        {
            print $pass join("\t", $id, @obs, $taxon), "\n";
            for ( my $i=0; $i<=$#obs; $i++ ) { $libPass[$i] += $obs[$i] }
        }
        else
        {
            print $fail join("\t", $id, @obs, $taxon), "\n";
            for ( my $i=0; $i<=$#obs; $i++ ) { $libFail[$i] += $obs[$i] }
        }
    }
    close($in);
    close($out);
    close($pass);
    close($fail);
    unlink(@outFiles);
    open (my $log, '>>', $logFile) or confess("Unable to append $logFile: $!");
    my $numClassified = sum(\@libPass);
    my $numUnclassified = sum(\@libFail);
    my $numClustered = $numClassified + $numUnclassified;
    my $numClassifiedPct = $numClustered ? int($numClassified/$numClustered*1000+0.5)/10 . '%' : 'N/A';
    my $numUnclassifiedPct = $numClustered ? int($numUnclassified/$numClustered*1000+0.5)/10 . '%' : 'N/A';
    my @libPassPct;
    my @libFailPct;
    for (my $i=0; $i<=$#libPass; $i++)
    {
        $libPassPct[$i] = $libClustered->[$i] ? int($libPass[$i]/$libClustered->[$i]*1000+0.5)/10 . '%' : 'N/A';
        $libFailPct[$i] = $libClustered->[$i] ? int($libFail[$i]/$libClustered->[$i]*1000+0.5)/10 . '%' : 'N/A';
    }
    print $log
        join("\t", 'CLASSIFIED', _commify($numClassified), map { _commify($_) } @libPass), "\n",
        join("\t", 'CLASSIFIED PCT', $numClassifiedPct, @libPassPct), "\n",
        join("\t", 'UNCLASSIFIED', _commify($numUnclassified), map { _commify($_) } @libFail), "\n",
        join("\t", 'UNCLASSIFIED', $numUnclassifiedPct, @libFailPct), "\n";
    close($log);
    return \@libPass;

    my $dur = time - $start;
    print "RDP Classifier: $numClassified/$numClustered; $dur sec\n";
}

=item _runRdp

Run RDP Classifier

=cut

sub _runRdp
{
    my ($rdp, $in, $out, $train, $min) = @_;
    _run("java -jar $rdp -q $in -o $out -t $train --minWords $min");
}

=item _taxFilter

Filter OTU table by taxonomy

=cut

# ECCE FIX RETURN VALUE -- RETURNING TOTAL NUMBER OF CLUSTERED SEQS, NEED TOTAL NUMBER OF CLUSTERS.

sub _taxFilter
{
    my ($config, $inFile, $outFile, $logFile, $libClassified) = @_;
    confess("Missing args") unless $config and $inFile and $outFile and $logFile and $libClassified;
    my $taxFilter = _set($config, 'RDP_CLASSIFIER', 'FILTER', undef);
    my $re;
    my @filter = split(/\|/, $taxFilter);
    foreach my $filter (@filter)
    {
        confess("Invalid filter string, $filter") unless $filter =~ /^[kpcofg]__\w+/i;
    }
    $re = qr/^$taxFilter/i;

    my $numLibs = scalar(@$libClassified);
    my @libPass = my @libFail = ( (0) x $numLibs );
    open(my $in, '<', $inFile) or confess("Unable to read $inFile: $!");
    open(my $out, '>', $outFile) or confess("Unable to write $outFile: $!");
    my $line = <$in>;
    print $out $line; # header
    while ($line = <$in>)
    {
        my @row = split(/\t/, $line);
        my $id = shift @row;
        my $taxon = pop @row;
        chomp $taxon;
        if ( $taxon =~ $re )
        {
            print $out $line;
            for (my $i=0; $i<=$#libPass; $i++) { $libPass[$i] += $row[$i] }
        } else
        {
            for (my $i=0; $i<=$#libFail; $i++) { $libFail[$i] += $row[$i] }
        }
    }
    close($in);
    close($out);
    my $totPass = sum(\@libPass);
    my $totFail = sum(\@libFail);
    my $totClassified = sum($libClassified);
    my $totPassPct = $totClassified ? int($totPass/$totClassified*1000+0.5)/10 . '%' : 'N/A';
    my $totFailPct = $totClassified ? int($totFail/$totClassified*1000+0.5)/10 . '%' : 'N/A';
    my @libPassPct;
    my @libFailPct;
    for (my $i=0; $i<=$#libPass; $i++)
    {
        $libPassPct[$i] = $libClassified->[$i] ? int($libPass[$i]/$libClassified->[$i]*1000+0.5)/10 . '%' : 'N/A';
        $libFailPct[$i] = $libClassified->[$i] ? int($libFail[$i]/$libClassified->[$i]*1000+0.5)/10 . '%' : 'N/A';
    }
    open(my $log, '>>', $logFile) or confess("Unable to append $logFile: $!");
    print $log
        join("\t", 'TAX-FILTER PASS', _commify($totPass), map { _commify($_) } @libPass), "\n",
        join("\t", 'TAX-FILTER PASS PCT', $totPassPct, @libPassPct), "\n",
        join("\t", 'TAX-FILTER FAIL', _commify($totFail), map { _commify($_) } @libFail), "\n",
        join("\t", 'TAX-FILTER FAIL PCT', $totFailPct, @libFailPct), "\n";
    close($log);
    return $totPass;
}

=item _run

Private function to run external executable.

=cut

sub _run
{
    my ($cmd) = @_;
    my $output = `$cmd 2>&1`;
    confess("FAILURE RUNNING $cmd : $output") unless $? == 0;
}

=item _set

Return value from config or default value. If default not defined, it must be defined in the config file, otherwise logconfess.

=cut

sub _set
{
    my ($config, $part, $key, $default) = @_;
    confess("Missing args") unless $config and $part and $key;
    unless ( exists($config->{$part}) )
    {
        return $default if $default;
        confess("Config file missing section, $part");
    }
    if ( exists($config->{$part}->{$key}) ) { return $config->{$part}->{$key} } 
    elsif ( defined($default) ) { return $default }
    else { confess("Config missing required $part/$key") }
}

=item _commify

Add commas to integers for readability

=cut

sub _commify
{
    my $x = shift;
    return 0 unless $x;
    my $integer = reverse $x;
    $integer =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $integer;
}


1;

=back

=head1 LIMITATIONS

The clustering is single-threaded, but we have attempted to speed it up by interspersing with multi-threaded mapping steps.  Consequently, the full number of reserved cores are not used at all times, resulting in some waste of resources.

=head1 NOTES

Special handling of 454 data (i.e. homopolymer-insensitive clustering) was dropped from this version.

=head1 AUTHORS

Edward Kirton and Julien Tremblay.
USEARCH by Robert Edgar (http://drive5.com)
RDP Classifier by James R. Cole, Qiong Wang, James M. Tiedje.

=head1 COPYRIGHT

Copyright (c) 2013 US DOE Joint Genome Institute. Use freely under the same license as Perl itself.  Refer USEARCH and RDP Classifier documentation for their own copyright and license information.

=cut

__END__
