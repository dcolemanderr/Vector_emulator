#!/global/home/users/colemand/ActivePerl-5.18/bin/perl

=pod

=head1 NAME

itagger.pl

=head1 PURPOSE

Pipeline for amplicon analysis (e.g. 16S rRNA, fungal ITS)

=head1 SYNOPSIS

    itagger.pl --config C<path> --lib C<path> --dir C<path>

=head1 OPTIONS

=over 5

=item --config C<path>

Config file in Ini format; see below.

=item --lib C<path>

Table of library name and location of Fastq file.  File should contain all and only reads associated with a library (i.e. one barcode of one lane of one run).

=item --dir C<path>

Base folder for output.

=item --group C<string>

Upon successful completion, set group for all folders, files; optional.

=item --tarball C<path>

Compress output; does not delete uncompressed output dir.  Optional.

=back

=head1 CONFIGURATION FILE

All paths and parameters are included in the configuration file, in Ini format, which has sections described below and are formatted thusly:

    [AMPLICON]
    NAME=16S_V4

=head2 AMPLICON

Information about the amplicon being sequenced.

=over 5

=item NAME = C<string>

Name of the amplicon; examples: 16S_V, ITS.  Required.

=item LEN_MEAN = C<int>

Average length of the amplicon sequence, as determined from Reference database.  That is, base-pairs between primers.  Required.

=item LEN_STDEV = C<int>

Standard deviation in amplicon length.  Required.

=item TRIM5 = C<int>

The number of bp from the 5' end of the amplicon to exclude from clustering and classification.  That is, bases immediately following the forward primer.  Allows exclusion of moderately conserved sequence which are problematic for classification.  For example, used for ITS but not 16S analysis.  Optional.

=item TRIM3 = C<int>

The number of bp from the 3' end of the amplicon to trim (i.e. bases following reverse primer).  Optional.

=back

=head2 DUK

We currently use C<duk> to remove contaminant reads; namely, sequence associated with sequencing library construction (control and adapter sequences).

=over 5

=item CONTAM_DB = C<path>

Multiple reference files in Fasta format are supported. Typically these would be QC sequences (e.g. PhiX) and adapter sequences.  Providing these as separate files will result in their numbers being reported separately. Optional.

=item KMER = C<int>

Kmer size for duk search.  Higher number increass specificity at expense of sensitivity.  Optional; default=21.

=item STEP = C<int>

Step size for search.  Lower number increases sensitivity but is slower.  Optional; default=1.

=item CUTOFF = C<float>

Cutoff detection level score.  Optional; default=1.

=back

=head2 CUTADAPT

C<cutadapt> is used to remove primer sequences.  For paired reads, if the primer is not found in both foward and reverse reads, the read-pair is filtered.

=over 5

=item PRIMER1 = C<string>

Forward (read1) primer; may contain ambiguous nucleotide characters.  Must be in same strand as read1.  Required.

=item PRIMER2 = C<string>

Reverse (read2) primer; may contain ambiguous nt chars.  Must be in same strand as read2.  Required.

=item ERROR_RATE = C<float>

Maximum fraction of errors allowed in detected primer sequence.  Increasing this number will allow more matches, but may increase false positives.  As we expect and require finding primer matches, a high error rate is recommended.  Optional; default=0.4.

=item MIN_OVERAP = C<int>

The minimum length of primer that must be found in the read to be considered a good hit.  Optional; default=10.

=back

=head2 Merging Paired Reads

It is assumed the amplicon primers were designed such that the read-pairs will overlap, to provide complete coverage of the amplicon sequence.  Paired-end reads are merged (assembled) into a single sequence using either C<Flash> or C<Pandaseq>.  Pairs which are not combined are discarded.

=head3 FLASH

=over 5

=item MIN_OVERLAP = C<int>

Minimum bp overlap required to merge paired reads.  Optional; default=10.

=item MAX_MISMATCH = C<float>

Maximum fraction of errors allowed in the overlap region.  Optional; default=0.3.

=back

=head3 PANDASEQ

=over 5

=item THRESHOLD = C<float>

Matching threshold (0-1); optional (default=0.5).

=item MIN_OVERLAP = C<int>

Minimum overlap in base pairs; optional (default=10).

=back

=head2 AMPQC

Quality filtering of amplicons prior to clustering and classification.

=over 5

=item EXP_ERR = C<float>

Filter amplicons whose total expected number of errors exceeds this threshold.  Optional; default=1.0.

=back

=head2 USEARCH

All libraries' equences are dereplicated and sorted by decreasing total normalized size before being clustered using C<USEARCH>.  We cluster iteratively, starting from a low centroid radius.  The clustering step is partially parallelized for speed.

=over 5

=item MIN_LIB_SIZE = C<int>

Discard libraries with less than this number of sequences.  Optional; default=1000.

=item MAX_RADIUS = C<int>

Sequences with less than this distance (percent sequence mismatch) from the centroid are added to the cluster.  Optional; default=3.

=item REF_DB = C<path>

Reference database of known amplicon sequences, used for chimera checking. The database should be of high quality and not include chimeric sequences.  It is recommended that the reference sequences are trimmed to the amplicon sequence (i.e. electronic PCR) and dereplicated.  They may also be clustered to reduce the database size.  Required.

=back

=head2 OTU_FILTER

=over 5

=item MIN_NORM_SIZE = C<int>

An OTU must have C<MIN_NUM_LIBS> different libraries with at least C<MIN_NORM_SIZE> total normalized reads, otherwise the OTU is filtered.  Optional; default=10.

=item MIN_NUM_LIBS = C<int>

See C<MIN_NORM_SIZE>.  Optional; default=2.

=back

=head2 RDP_CLASSIFIER

Parameters for RDP Classifier and for filtering clusters from further analysis.  For example, to discard chloroplast sequences from 16S analysis.

=over 5

=item TRAINING_FILE = C<path>

Amplicon-specific RDP Classifier training file.  Required.

=item LEVEL = C<string>

An OTU must be classified at this level with a confidence above cutoff (below) in order to pass.  Optional; default="class".

=item CUTOFF = C<float>

Filter OTUs which do not meet minimum confidence level at the specified taxonomic classification (above).  Optional default=0.5.

=item FILTER = C<string>

To filter OTUs which do not match that described in search string.  String follows the pattern: X__NAME, where X indicated taxonomic level (k=kingdom, p=phylum, c=class, o=order, f=family, g=genus).  Multiple search strings should be separated by pipe ("|") characters with no spaces.  Use quotation marks if your search string contains any spaces.  Examples: k__bacteria|k__archaea, k__fungi.  Case insensitive.  Optional.

=back

=cut

=head2 QIIME

Analysis of OTU tables using C<QIIME>.

=over 5

=item RAREFY_MIN = C<float>

Cutoff for single-rarefaction.  Optional; by default use 10% of normalized mean.

=back

=head2 DIVERSITY I<(beta-testing at this time)>

Diversity plots and estimates using C<R> scripts.  Run only if C<FACTORS> is defined.

=over 5

=item FACTORS = C<PATH>

Table of libraries X factors.  Library IDs must match exactly to those used in the OTU table (i.e. C<DATA> section above).

=back

=head1 AUTHORS

Julien Tremblay (julien.tremblay@mail.mcgill.ca),
Edward Kirton (ESKirton@LBL.gov)

=head1 COPYRIGHT

Copyright (c) 2013 United States Department of Energy Joint Genome Institute.  Use freely under the same license as Perl itself.  Refer to wrapped tools for their own licenses and copyright information.

=cut

print "Running\n\n";

use strict;
use warnings;
use Env qw/TMPDIR RDP_JAR_PATH/;
use Getopt::Long;
use Pod::Usage;
use File::Spec::Functions qw(rel2abs);
use Config::Tiny;
use FindBin qw($Bin);
use Parallel::ForkManager;
use iTagger::ReadQC qw(libQC);
use iTagger::ReadQC::Meta qw(summary);
use iTagger::ClusterOTUs qw(clusterOTUs);
use iTagger::QIIME qw(qiimeReport);

our $VERSION = "1.1";

my ($help, $man, $configFile, $libsFile, $dir, $group, $tarball);
my $threads = `nproc`;
chomp $threads;
GetOptions(
    'config=s' => \$configFile,
    'libs=s' => \$libsFile,
    'dir=s' => \$dir,
    'group=s' => \$group,
    'threads=i' => \$threads,
    'tarball=s' => \$tarball,
    'help' => \$help,
    'man' => \$man ) or pod2usage(2);
pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;
die("--threads required\n") unless $threads;

# LOGGER
print "Begin iTagger version $VERSION with $threads threads\n";

# ENV
die("\$TMPDIR not defined") unless defined($TMPDIR);
die("\$RDP_JAR_PATH not defined") unless defined($RDP_JAR_PATH);

# OUTPUT DIR
die("--dir path required") unless $dir;
$dir = rel2abs($dir);
-d $dir or mkdir($dir) or die("Unable to mkdir $dir: $!");

# CONFIG FILE
die("--config file required") unless $configFile;
die("--config file not found: $configFile") unless -f $configFile;
$configFile = rel2abs($configFile);
my $config = Config::Tiny->read($configFile);
foreach (qw/AMPLICON AMPQC USEARCH RDP_CLASSIFIER/)
{
    die("Config file missing section: $_") unless exists($config->{$_});
}
$config->{_}->{VERSION} = $VERSION;
$config->write("$dir/config.ini");

# LIBS
die("--libs file required") unless $libsFile;
die("--libs file not found: $libsFile") unless -f $libsFile;
$libsFile = rel2abs($libsFile);
my %libs; # libName => [ files ]
my %files; # file => undef
open( my $in, '<', $libsFile) or die($!);
while (<$in>)
{
    chomp;
    my @row = split(/\t/);
    die("Invalid libs file") unless @row == 2;
    my ($name, $file) = @row;
    die("Duplicate reads file, $file") if exists($files{$file});
    $files{$file} = undef;
    $name =~ s/\W/_/g;
    $libs{$name} = [] unless exists($libs{$name});
    push @{$libs{$name}}, $file;
}
close($in);

# CREATE LIST OF INPUT, IDENTIFY LIBS WITH REPS AND MAKE NAMES UNIQ
my @input; # [ file, libName ]
open(my $out, '>', "$dir/libs.tsv") or die($!);
foreach my $name ( sort keys %libs )
{
    my $numReps = scalar(@{$libs{$name}});
    if ( $numReps == 1 )
    {
        my $file = $libs{$name}->[0]; 
        push @input, [ $name, $file ];
        print $out "$name\t$file\n";
    } else
    {
        for ( my $i=1; $i<=$numReps; $i++ )
        {
            my $file = $libs{$name}->[$i-1]; 
            push @input, [ "$name.$i", $file ];
            print $out "$name.$i\t$file\n";
        }
    }
}
close($out);

# READ QC
my $pm = new Parallel::ForkManager($threads);
-d "$dir/data" or mkdir("$dir/data") or die($!);
my $qc = 0;
while ( @input )
{
    my $ar = shift @input;
    my ($name, $file) = @$ar;
    next if -e "$dir/data/$name/readQC.log";
    unless ( -e $file )
    {
        warn("Lib $name file not found: $file");
        next;
    }
    ++$qc;
    $pm->start and next;
    my $qc = new iTagger::ReadQC($configFile, "$dir/data", $file, $name);
    $qc->libQC();
    $pm->finish;
}
$pm->wait_all_children;
my $report = new iTagger::ReadQC::Meta($configFile, "$dir/data");
$report->summary();

# CLUSTERING
my $otuFile = "$dir/otu/envcom.pass.tsv";
if ( $qc or ! -e $otuFile )
{
    clusterOTUs($configFile, $dir, $threads);
}

# QIIME
qiimeReport($configFile, $dir);

# README
open($out, '>', "$dir/README.txt") or die($!);
while (<DATA>) { print $out $_ }
close(DATA);
print $out 'duk: ', `duk | grep Version | tail -1 | cut -f 1 -d ':'`;
print $out 'cutadapt: ', `cutadapt --version`;
print $out 'flash: ', `which flash`;
print $out `pandaseq -v 2>&1`;
print $out `usearch --version`;
print $out "RDP Classifier: $RDP_JAR_PATH\n";
print $out 'QIIME: ', `which alpha_diversity.py`;
close($out);

# FIX PERMS
if ( $group )
{
    system("chgrp -R $group $dir");
    system("chmod -R a+rX $dir");
    system("chmod -R ug+rwX $dir");
}

# TARBALL
if ( $tarball )
{
    system("tar -cvzf $tarball $dir");
}
exit;
__DATA__
iTagger v1.1 METHODS

The iTagger amplicon analysis pipeline uses several publicly available tools to analyze amplicon libraries, such as 16S rRNA or fungal ITS variable regions for phylogenetic analysis.  All libraries to be compared should be identically constructed, sequenced, and analyzed.

I. INPUT:

(1) Configuration file in INI format with parameters and paths of reference databases.
    OUTPUT: The config file is copied to <OUTDIR>/config.ini
(2) Libraries tabular file indicates library/condition name and path to Fastq file.

II. READ QC:

Each library's Fastq file is processed as described below and the results are saved in the data/ folder, which also includes a readQC.log file which indicates the read-pairs at each step and the percentage of pairs which pass each step (i.e. percentages are per-step, not of total input).
    OUTPUT: Saved in <OUTDIR>/data/<LIB_NAME>/
            Summarized in readQC.log

(1) CONTAM FILTER: Filter one or more contaminants using Duk (e.g. PhiX control, sequencing library adapter dimers, human contaminants, etc.).  For paired reads, the entire pair is filtered if one end-read has a high-scoring hit.
    OUTPUT: duk.log
(2) PRIMER TRIM: PCR primers (of the conserved region) are trimmed away.  For paired-end reads, both forward and reverse primers must be found, otherwise the entire pair is filtered.
(3) HARD TRIM: For particular libraries (e.g. fungal ITS), it is useful to trim a predefined number of bases from the 5' and 3' ends of the sequence.  For fungal ITS, we observed better RDP Classifier results after trimming conserved regions.  We do not recommend hard trimming for 16S sequences.
(4) ITERATIVE PAIR MERGING: Reads are trimmed as a pair, removing the last base from whichever end has the highest expected error in a window 5bp wide.  Reads are trimmed from mean + 3 standard deviations to mean - 3 standard deviations in 0.5 standard deviation steps.  After each trimming step, pairs are merged into single sequences with either Flash or Pandaseq.  Pairs which are not merged continue to the next round of trimming.  Paired reads which are not combined are discarded.
    OUTPUT: Not-combined reads, nc.fastq.bz2
(5) EXPECTED-ERROR FILTER: Merged reads are filtered if they have an expected number of errors which exceeds the threshold.  The config file indicates the maximum number of expected errors per 100bp. Note that Flash and Pandaseq produce different quality scores in the overlap-assembled regions.
    OUTPUT: Filtered extended reads, ex.fastq.bz2
            Quality report, qualStat.pdf, qualStat.tsv
(7) DEREPLICATE: Count the number of times each sequence is observed and output in tabular (seq-obs) format, ordered by sequence.
    OUTPUT: seqobs.tsv.bz2

III. CLUSTERING:

USEARCH is used for clustering, although there is a provision for iterative clustering which can (a) provide faster processing and (b) allow processing of larger files than can normally be processed (particularly with the 32bit version).  RDP Classifier is used for taxonomic classification of the resultant cluster centroid sequences and it's accuracy is dependent upon providing a well-curated RDP reference database.
    OUTPUT: Saved in <OUTDIR>/otu
            Summarized in cluster.log

(1) MERGE LIBRARIES: The seq-obs files for all libraries are merged, dereplicated, and sorted by decreasing abundance.  Low-abundance sequences are separated and excluded from clustering, step 2 (although they will be mapped and counted in step 3).
(2) ITERATIVE CLUSTER OTUS: Refer to the USEARCH documentation for the algorithm description.  Our use of USEARCH differs slightly from that described in the USEARCH documentation in that we iterate between single-threaded clustering and multi-threaded searching in order to reduce run-time. We also use .obs files for tracking cluster members, so a final mapping and counting step (as described in the USEARCH docs) is unnecessary.  Clustering is done iteratively starting at 99% identity, and decreasing by 1% identity until the level described in the config file is reached (e.g. 97% for 16S, 95% for Fungal ITS).
(3) MAP LOW-ABUNDANCE SEQUENCES: Rare sequences, which cannot form their own clusters, are mapped to the cluster centroid sequences and counted.
(4) REFERENCE DB CHIMERA FILTER: Centroid sequences are compared to the reference database and likely chimeric sequences discarded, using UCHIME.
    OUTPUT: Final cluster centroids, otu.fasta.bz2
(5) CLASSIFICATION: Assign taxonomic classification to each cluster using RDP Classifier.  The config file indicates a taxonomic level (e.g. family) and confidence level (e.g. 0.5) which is used to decide which classifications are useful.  Clusters which can be acceptably classified are output to otu.tax.tsv, while the others are written to otu.unk.tsv.
    OUTPUT: RDP output, rdp.tsv
            Classified OTUs, otu.tax.tsv
            Unclassified OTUs, otu.unk.tsv
(6) TAX FILTER: Clusters with classifications which do not match those indicated in the config file are discarded and the desired clusters are written to the otu.tax.filtered.tsv file.
    OUTPUT: Final OTU table, otu.tax.filtered.tsv

IV. TAXONOMIC ANALYSIS:

QIIME is used to manipulate the final OTUs file.  A few graphs are generated plus some rarefied tables which may be useful for subsequent analysis.

(1) GENERATE BIOM: The BIOM JSON file is generated from the OTU tabular file.
    OUTPUT: otu.biom
(2) ABUNDANCE THRESHOLD: Filter OTUs are assorted levels and calculate alpha diversities.
    OUTPUT: Several files in <OUTDIR>/abundance_thresholds/
(3) SINGLE RAREFACTION: This is done at 1000 and at a level calculated from the trimmed mean and standard deviation of the library sizes (10% trimmed, and calc. mean - i * stdev; while i is the highest number from 2 to 0.5, step 0.5, until the cutoff is above 0).
    OUTPUT: <OUTDIR>/otu/rarefied.1000.biom, rarefied.1000.filtered.biom
            <OUTDIR>/otu/rarefied.<X>.biom, rarefied.<X>.filtered.biom
(4) SUMMARIZE TAXONOMY: with both relative and absolute abundance
    OUTPUT: <OUTDIR>/tax_mapping/relative/
            <OUTDIR>/tax_mapping/absolute/
(5) PLOT RANK-ABUNDANCE: Generate PDF rank-abundance graph of all samples.
    OUTPUT: <OUTDIR>/otu/log_rank_abundance.pdf
(6) PLOT TAXA SUMMARY: Make taxa plot of absolute abundance
        OUTPUT: <OUTDIR>/tax_mapping/plots/
(7) PHYLUM BARPLOT: generate a phylum-level barplot using absolute abundance for a quick overview of the data.
        OUTPUT: <OUTDIR>/tax_mapping/taxonomy_phylum_L2.tab

* * *

SUMMARY OF OUTPUT FILES:

(1) <OUTDIR>/

    config.ini : parameter settings used
    libs.txt : libraries/samples and their sequence files
    README.txt : description of method

(2) <OUTDIR>/data/

    readQC.log : summary of read preprocessing, one column per library/sample
        Fields:
        "Input" : starting number of read-pairs
        "Filtered" : read-pairs remaining after named filter (e.g. PhiX)
        "Primer-trimmed" : read-pairs remaining after primer trimming
        "Extended" : number of read-pairs successfully merged into consensus
        "Length-filtered" : reads after short reads removed
        "High quality" : reads which pass expected-error filter
        "Unique" : number of unique sequences (includes singlets)

(3) <OUTDIR>/data/<LIB_NAME>/

    readQC.log : summary statistics for this library/sample
    noContam.fastq : reads after contaminants (e.g. adapter) have been filtered
    duk.log : contaminant filter log

(4) <OUTDIR>/OTU

    cluster.log : summary of clustering statistics -- total and for each library
        Fields:
        "Input" : number of dereplicated sequences provided
        "Clusterable" : number of sequences used in clustering (excludes singlets)
        "Singlets" : singlets cannot form clusters but can join clusters
        "Singlets Mapped" : number of singlets assigned to clusters
        "Ref-chimera" : clusters discarded as chimeric by comparison to ref db
        "Clustered" : total number of reads assigned to remaining clusters
        "Classified" : number of reads which belong to classifiable clusters
        "Unclassified" : number of reads belonging to unclassifiable clusters
        "Tax-filter pass" : reads associated with clusters of specified taxa
    otu.fasta.bz2 : sequences of cluster representatives
    rdp.tsv : RDP Classifier output
    otu.tax.tsv : abundance table with assigned taxonomies
    otu.unk.tsv : abundance table for unclassified OTUs only
    otu.tax.filtered.tsv : abundance table for OTUs of desired clades only
    otu.biom : abundance table in biom json format (e.g. for QIIME)
    log_rank_abundance.pdf : rank-abundance curve of all libs, generated by QIIME
    rarefied.1000.biom : biom json file with only 1000 reads sampled per library, generated by QIIME's single_rarefaction.py
    rarefied.1000.filtered.biom : above, with low-abundance OTUs removed
    rarefied.<X>.biom : single rarefaction with sample size calculated from sample stats (see methods above)
    rarefied.<X>.filtered.biom : above, with low-abundance OTUs removed

(5) <OUTDIR>/tax_mapping

    taxonomy_phylum_L2.pdf : phylum barplot for a quit view of the data
    taxonomy_phylum_L2.tab : supporting datafile for above
    taxonomy_phylum_L2.tab.perc : supporting datafile for above

(6) <OUTDIR>/tax_mapping/plots/

    bar_charts.html : HTML taxonomy report generated by QIIME
    charts/ css/ js/ raw_data/ : supporting files, included in above

(7) <OUTDIR>/tax_mapping/absolute/

    otu_L<X>.txt : absolute abundance taxa plot at multiple taxonomic levels, generated by QIIME's plot_taxa_summary.py
    otu_L<X>.biom : biom json format (e.g. for QIIME)

(8) <OUTDIR>/tax_mapping/relative/

    otu_L<X>.txt : relative abundance taxa plot at multiple taxonomic levels, generated by QIIME's plot_taxa_summary.py
    otu_L<X>.biom : biom json format (e.g. for QIIME)

* * *

iTagger was written by Julien Tremblay (julien.tremblay@mail.mcgill.ca) and Edward Kirton (ESKirton@LBL.gov) and is Copyright (c) 2013 by the US DOE Joint Genome Institute but is freely available for use without any warranty under the same license as Perl itself.  v1.1 was released 12/12/2013.  Refer to wrapped tools for their author and license information.

* * *

External executables used:

