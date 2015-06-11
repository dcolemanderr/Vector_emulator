=pod

=head1 NAME

iTagger::ReadQC

=head1 DESCRIPTION

This module contains functions for preprocessing reads prior to clustering.  It should receive a FASTQ file of raw sequence data corresponding to one barcode of one lane (ie. all seqs belong to same library or experimental condition).

=head1 METHODS

=cut

package iTagger::ReadQC;

use strict;
use warnings;
use Carp;
use Env qw(TMPDIR ITAGGER_LOG_CONFIG);
use File::Copy;
use Config::Tiny;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(fileparse);
use constant { MIN_AMPLICON_LEN => 50, MIN_FINAL_LEN => 32, TRIM_NUM_STDEVS => 2.5 };
require Exporter;
use iTagger::FastqDb;
use iTagger::FastqDb::Tiny;
use iTagger::Stats;

our $VERSION = 1.2;
our @ISA = qw(Exporter);
our @EXPORT = qw(new libQC);

=head2 CONSTRUCTOR AND DESTRUCTOR

=over 5

=item new

Constructor.  Config is either a hashref or the path to an .ini format file.

=cut

sub new
{
    my ($class, $config, $dir, $inFile, $libName) = @_;

    # INIT
    my $startTime = time;
    my $this =
    {
        name    => $libName,
        infile  => $inFile,
        start   => $startTime,
        logfile => undef,
        logfh   => undef
    };
    bless $this, $class;

    # PARAMS
    confess("Missing args") unless $config and $dir and $inFile and $libName;

    # OUTPUT DIR
    -d $dir or mkdir($dir) or confess("Unable to mkdir $dir: $!");
    $dir .= "/$libName" unless $dir =~ /$libName$/;
    -d $dir or mkdir($dir) or confess("Unable to mkdir $dir: $!");
    $this->{dir} = $dir;
    $this->run("mkdir -p $dir");

    # TMPDIR
    my $tmpdir = -e '/scratch' ? '/scratch' : $TMPDIR;
    $this->{tmpdir} = $tmpdir = "$tmpdir/$$.$startTime";
    -d $tmpdir or mkdir($tmpdir) or confess("Unable to mkdir $tmpdir: $!");

    # CONFIG
    if ( ref($config) eq '' ) { $config = Config::Tiny->read($config) }
    elsif ( ref($config) ne 'HASH' ) { confess("Missing config file/hash") }
    foreach my $section (qw/AMPLICON AMPQC/)
    {
        confess("Config file missing required section: $section") unless exists($config->{$section});
    }
    if ( exists($config->{FLASH}) and exists($config->{PANDASEQ}) )
    {
        confess("Invalid config file: FLASH xor PANDASEQ");
    }
    $this->{config} = $config;
    
    # LOGFILE
    my $logFile = $this->{logfile} = "$dir/readQC.log";

    # PAIRED OR NOT?
    $this->{paired} = ( exists($config->{FLASH}) or exists($config->{PANDASEQ}) ) ? 1 : 0;

    # INIT LOGFILE (OVERWRITES EXISTING)
    my $tmpLogFile = "$logFile.tmp";
    open(my $log, '>', $tmpLogFile) or confess($!);
    print $log "#STEP\tCOUNT\tPCT OF PREVIOUS STEP\n";
    $this->{logfh} = $log;

    return $this;
}

=item DESTROY

Renames logfile before exiting.

=cut

sub DESTROY
{
    my $this = shift;

    # CLEANUP TMPDIR
    $this->run("rm -rf $this->{tmpdir}");

    # RENAME LOG FILE
    my $fh = $this->{logfh};
    close($fh) if defined($fh);
    my $logfile = $this->{logfile};
    if ( $logfile and -e "$logfile.tmp" )
    {
        move("$logfile.tmp", $logfile) or confess($!);
    }
}

=back

=head2 DATA PROCESSING

=over 5

=item libQC

Main function to QC a library of reads.  All reads should belong to the same experimental condition (library).  Parameters come from required configuration file.

=cut

sub libQC
{
    my ($this) = @_;

    # INIT VARS
    my $dir = $this->{dir};
    my $tmpdir = $this->{tmpdir};
    my $config = $this->{config};

    # INPUT: VALIDATE, UNCOMPRESS IF NECESSARY, CONVERT TO SANGER-ENCODING IF NECESSARY, AND COUNT.
    # TODO: IF UNCOMPRESSED AND ALREADY SANGER-ENCODED, SIMPLY COUNT AND SYMLINK, DO NOT COPY.
    # TODO: DELETE ITAGGER::FASTQDB AND USE ITAGGER::FASTQDB::TINY EVERYWHERE
    unless ( -e $this->{infile} )
    {
        carp("Lib $this->{name} input file not found: $this->{infile}");
        return;
    }
    unless ( -s $this->{infile} )
    {
        $this->log('Input', 0);
        carp("Lib $this->{name} input file is empty: $this->{infile}");
        return;
    }
    my $inFile = "$tmpdir/input.fastq";
    my $numInput = 0;
    my $db = new iTagger::FastqDb($this->{infile});
    open(my $out, '>', $inFile) or confess($!);
    if ( $this->{paired} )
    {
        while ( my $pair = $db->next_pair )
        {
            ++$numInput;
            print $out $pair->[0]->output, $pair->[1]->output;
        }
    } else
    {
        while ( my $seq = $db->next_seq )
        {
            ++$numInput;
            print $out $seq->output;
        }
    }
    close($out);
    $this->log('Input', $numInput);
    return unless $numInput;

    # FILTER CONTAMINANTS.  THERE MAY BE ZERO OR MORE CONTAM DB FILES.
    my $noContamFile;
    my $numNoContam;
    if ( exists($config->{DUK}->{CONTAM_DB}) )
    {
        $noContamFile = "$dir/noContam.fastq";
        my $numNoContamAR;
        my $contamDbNamesAR;
        ($numNoContamAR, $contamDbNamesAR) = $this->duk($inFile, $noContamFile, "$dir/duk.log");
        my $n0 = $numInput;
        for ( my $i=0; $i<= $#$contamDbNamesAR; $i++ )
        {
            my $dbName = $contamDbNamesAR->[$i];
            $numNoContam = $numNoContamAR->[$i];
            $this->log("Filtered $dbName", $numNoContam, $n0);
            $n0 = $numNoContam;
        }
    } else
    {
        $noContamFile = "$tmpdir/noContam.fastq";
        $numNoContam = $numInput;
        symlink($inFile, $noContamFile) or confess($!);
    }
    return unless $numNoContam;

    # TRIM PRIMERS USING CUTADAPT
    my $primerTrimPassFile = "$tmpdir/primerTrim.pass.fastq";
    my $primerTrimFailFile = undef; #"$dir/primerTrim.fail.fastq.bz2";
    my $numPrimerTrim;

    if ( exists($this->{config}->{CUTADAPT}) and (
        ( exists($this->{config}->{CUTADAPT}->{PRIMER1}) and $this->{config}->{CUTADAPT}->{PRIMER1} ) or
        ( exists($this->{config}->{CUTADAPT}->{PRIMER2}) and $this->{config}->{CUTADAPT}->{PRIMER2} ) ) )
    {
        ($numPrimerTrim) = $this->cutadapt($noContamFile, $primerTrimPassFile, $primerTrimFailFile);
        $this->log('Primer-trimmed', $numPrimerTrim, $numNoContam);
        return unless $numPrimerTrim;
    } else
    {
        $numPrimerTrim = $numNoContam;
        $this->log('Primer-trimmed', $numPrimerTrim, $numNoContam);
        symlink($noContamFile, $primerTrimPassFile) or confess("Symlink failure");
    }

    # TODO MOVE BELOW TO SEPARATE METHODS

    #########################
    ## UNPAIRED READS

    unless ( $this->{paired} )
    {
        # OPTIONAL HARD TRIM AND LENGTH FILTER
        my $lengthFilteredFile = "$this->{tmpdir}/lengthFiltered.fastq";
        my ($numLengthFiltered) = $this->lengthFilter($primerTrimPassFile, $lengthFilteredFile);
        $this->log('Length-filtered', $numLengthFiltered, $numPrimerTrim);
        return unless $numLengthFiltered;

        # EXPECTED ERROR FILTER
        my $hiQualFile = "$dir/qc.hiqual.fastq";
        my $numHiQual = $this->expErrFilter($lengthFilteredFile, $hiQualFile);
        $this->log("Quality-filtered", $numHiQual, $numLengthFiltered);
        return unless $numHiQual;

        # DEREPLICATE
        my $numUniq = $this->seqObs($hiQualFile, "$dir/qualStat.tsv", "$dir/qualStat.pdf", 'High quality sequences', "$dir/seqobs.tsv.bz2");
        $this->log("High quality, unique", $numUniq, $numHiQual);

        # FINISH
        my $dur = int((time - $this->{start})/60 * 10 + 0.5)/10;
        print "QC $this->{name} ($numInput unpaired reads) done after $dur min\n";
        return;
    }

    #########################
    ## PAIRED READS

    # MERGE PAIRED READS
    my $exFile = "$dir/ex.fastq.bz2";
    my $ncFile = "$dir/nc.fastq.bz2";
    my $numEx;
    my $numNc;
    if ( exists($config->{FLASH}) )
    {
        # MERGE PAIRED READS USING FLASH
        if ( -s $primerTrimPassFile )
        {
            ($numEx, $numNc) = $this->iterativeFlash($primerTrimPassFile, $exFile, $ncFile);
        } else
        {
            logcarp("THERE ARE NO SEQS IN PRIMERTRIMPASSFILE\n"); # ECCE
        }
    } else
    {
        # TRIMS PRIMERS AND MERGE PAIRED READS USING PANDASEQ
        ($numEx, $numNc) = $this->iterativePandaseq($primerTrimPassFile, $exFile, $ncFile, "$dir/pandaseq.log.bz2");
    }
    $this->log('Extended', $numEx, $numPrimerTrim);
    return unless $numEx;

    # OPTIONAL HARD TRIM AND LENGTH FILTER
    my $lengthFilteredFile = "$this->{tmpdir}/lengthFiltered.fastq";
    my ($numLengthFiltered) = $this->lengthFilter($exFile, $lengthFilteredFile);
    $this->log('Length-filtered', $numLengthFiltered, $numEx);
    return unless $numLengthFiltered;

    # EXPECTED ERROR FILTER
    my $hiQualFile = "$dir/hiqual.fastq";
    my $numHiQual = $this->expErrFilter($lengthFilteredFile, $hiQualFile);
    $this->log("High quality", $numHiQual, $numLengthFiltered);
    return unless $numHiQual;

    # DEREPLICATE
    my $numExUniq = $this->seqObs($hiQualFile, "$dir/qualStat.tsv", "$dir/qualStat.pdf", 'Extended and filtered sequences', "$dir/seqobs.tsv.bz2");
    $this->log("Unique", $numExUniq, $numHiQual);

    # FINISH
    my $dur = int((time - $this->{start})/60 * 10 + 0.5)/10;
    print "QC $this->{name} ($numInput paired reads) done after $dur min\n";
}

=item log

Append message to log summary stats file.

=cut

sub log
{
    my ($this, $label, $num, $denom) = @_;
    return unless $this->{logfile};
    $num = 0 unless defined($num);
    my $pct = $denom ? (int($num/$denom*1000+0.5)/10).'%' : '';
    my $numStr = reverse $num;
    $numStr =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    $numStr = scalar reverse $numStr;
    my $fh = $this->{logfh};
    print $fh $label, "\t", $numStr, "\t", $pct, "\n";
}

=item duk

Filter contaminant reads.  There may be one or more reference db files.  Each reference db is filtered separately (in order) so that the number of hits per db is known.

=cut

sub duk
{
    my ($this, $inFile, $outFile, $logFile) = @_;

    # INIT VARS
    my $config = $this->{config};
    return unless exists $config->{DUK}->{CONTAM_DB};
    my @contamDbs = split(/,/, $config->{DUK}->{CONTAM_DB});
    foreach (qw/KMER STEP CUTOFF/)
    {
        confess("Missing duk parameter $_") unless exists($config->{DUK}->{$_});
    }
    my $numOutput;
    my @numOutput = ();
    my @dbNames = ();
    my $tmpOutFile = "$this->{tmpdir}/duk.out.fastq";
    my $tmpLogFile   = "$this->{tmpdir}/duk.log";
    my $tmpInFile = "$this->{tmpdir}/duk.in.fastq";

    # SETUP
    symlink($inFile, $tmpInFile) or confess("Error symlinking $inFile: $!");

    # FILTER
    for ( my $i=0; $i<=$#contamDbs; $i++ )
    {
        # CONTAM DB
        my $contamDb = $contamDbs[$i];
        next unless $contamDb;
        confess("Contam db does not exist: $contamDb") unless -e $contamDb;
        my ($dbName, $dbDir, $dbSuffix) = fileparse($contamDb, qr/\.[^.]*/);
        push @dbNames, $dbName;

        # INPUT
        if ( $i > 0 )
        {
            move($tmpOutFile, $tmpInFile) or confess("Error moving $tmpOutFile: $!");
        }
        $numOutput = $this->_duk($contamDb, $tmpInFile, $tmpOutFile, $tmpLogFile);
        push @numOutput, $numOutput;
        unlink($tmpInFile);

        # DUK LOG
        if ( $i == 0 )
        {
            move($tmpLogFile, $logFile) or confess("Error moving $tmpLogFile: $!");
        } else
        {
            # APPEND LOG
            open(my $in, '<', $tmpLogFile) or confess("Error reading logfile: $!");
            open(my $out, '>>', $logFile) or confess("Error appending logfile, $logFile: $!");
            while (<$in>) { print $out $_ }
            close($in);
            close($out);
            unlink($tmpLogFile);
        }
        last unless -s $tmpOutFile;
    }
    if ( -s $tmpOutFile )
    {
        move($tmpOutFile, $outFile) or confess("Error moving $tmpOutFile: $!");
    } else
    {
        open(my $out, '>', $outFile) or confess($!);
        close($out);
    }
    return (\@numOutput, \@dbNames);
}

=item _duk

Run duk contaminant filter with a single db of contaminant sequences.  If paired reads then remove any singletons generated.

=cut

sub _duk
{
    my ($this, $dbFile, $inFile, $outFile, $logFile) = @_;
    confess("Missing args") unless $dbFile and $inFile and $outFile and $logFile;

    # CONFIG
    my $kmer   = $this->{config}->{DUK}->{KMER};
    my $step   = $this->{config}->{DUK}->{STEP};
    my $cutoff = $this->{config}->{DUK}->{CUTOFF};

    # RUN DUK
    my $tmpFile = "$this->{tmpdir}/duk.fastq";
    $this->run("duk -o $logFile -n $tmpFile -m /dev/null -k $kmer -s $step -c $cutoff $dbFile $inFile");

    # PARSE LOG TO DETERMINE NUMBER OF READS INPUT AND MATCHED
    my $numInput;
    my $numMatched;
    open(my $in, '<', $logFile) or confess("Error reading log, $logFile: $!");
    while (<$in>)
    {
        if (/^# Total number of reads:\s+(\d+)$/)
        {
            $numInput = $1;
        } elsif (/^# Total number of matched reads:\s+(\d+)$/)
        {
            $numMatched = $1;
            last;
        }
    }
    close($in);
    my $numOutput = $numInput - $numMatched;

    # DONE IF UNPAIRED READS
    unless ( $this->{paired} )
    {
        move($tmpFile, $outFile) or confess("Error moving $tmpFile: $!");
        return $numOutput;
    }

    # DISCARD SINGLETONS FOR PAIRED READS
    $numOutput /= 2 if $numOutput; # covert reads to pairs
    if ( -e $tmpFile )
    {
        # FASTQDB WILL DISCARD UNPAIRED READS WHEN DB'S PAIRED PARAM IS TRUE
        my $db = new iTagger::FastqDb::Tiny($tmpFile, { paired => 1 });
        open(my $out, '>', $outFile) or confess("Error writing $outFile: $!");
        while ( my $pair = $db->next_pair )
        {
            print $out $pair->[0]->output, $pair->[1]->output;
        }
        close($out);
        unlink($tmpFile);
    }
    return $numOutput;
}

=item cutadapt

Run cutadapt to remove primer sequences.

=cut

sub cutadapt
{
    my $this = shift;
    return $this->_cutadaptPaired(@_) if $this->{paired};
    my $strandInspecific = $this->_set('CUTADAPT', 'STRAND_INSPECIFIC', 0);
    return $this->_cutadaptUnpairedStrandInspecific(@_) if $strandInspecific;
    return $this->_cutadaptUnpaired(@_);
}

=item _cutadaptUnpaired

Trim primer sequences.  Unless all specified primers are found, the read is filtered.

=cut

sub _cutadaptUnpaired
{
    my ($this, $inFile, $outFile, $failFile) = @_;
    confess("Missing args") unless $inFile and $outFile;

    # VALIDATE PARAMETERS
    my $params = $this->{config}->{CUTADAPT};
    foreach (qw/ERROR_RATE MIN_OVERLAP PRIMER1 PRIMER2/)
    {
        confess("Missing cutadapt parameter, $_") unless exists($params->{$_});
    }
    my $errRate = $params->{ERROR_RATE};
    my $minOverlap = $params->{MIN_OVERLAP};

    # DEFINE ARGUMENTS
    my $primer1 = $params->{PRIMER1} ? '-g '.$params->{PRIMER1} : '';
    my $primer2 = $params->{PRIMER2} ? '-a '.$params->{PRIMER2} : '';
    confess("CUTADAPT requires PRIMER1 or PRIMER2") unless $primer1 or $primer2;

    # RUN CUTADAPT
    $this->run("cutadapt --trimmed-only -e $errRate -O $minOverlap $primer1 $primer2 -o $outFile $inFile > /dev/null");

    # GET FAILED READS
    my $numPass = my $numFail = 0;
    my $db = new iTagger::FastqDb::Tiny($inFile);
    my $ok = new iTagger::FastqDb::Tiny($outFile);
    my $failOut = $failFile ? "| bzip2 --best > $failFile" : '>/dev/null';
    open(my $out, $failOut) or confess($!);
    my $okSeq = $ok->next_seq;
    while ( my $seq = $db->next_seq )
    {
        if ( defined($okSeq) and $seq->id eq $okSeq->id )
        {
            ++$numPass;
            $okSeq = $ok->next_seq;
        } else
        {
            ++$numFail;
            print $out $seq->output;
        }
    }
    close($out);
    return ($numPass, $numFail);
}

=item _cutadaptUnpairedStrandInspecific

Trim primer sequences. If primers are not found, reverse-complement and try again.  Filter reads in which both primers not found in either strand.  This was written for PacBio sequences.

=cut

sub _cutadaptUnpairedStrandInspecific
{
    my ($this, $inFile, $outFile, $failFile) = @_;
    confess("Missing args") unless $inFile and $outFile;

    # VALIDATE PARAMETERS
    my $params = $this->{config}->{CUTADAPT};
    foreach (qw/ERROR_RATE MIN_OVERLAP PRIMER1 PRIMER2/)
    {
        confess("Missing cutadapt parameter, $_") unless exists($params->{$_});
    }
    my $errRate = $params->{ERROR_RATE};
    my $minOverlap = $params->{MIN_OVERLAP};
    my $primer1 = $params->{PRIMER1} ? '-g '.$params->{PRIMER1} : '';
    my $primer2 = $params->{PRIMER2} ? '-a '.$params->{PRIMER2} : '';
    confess("CUTADAPT requires PRIMER1 or PRIMER2") unless $primer1 or $primer2;

    # RUN CUTADAPT (ROUND 1)
    $this->run("cutadapt --trimmed-only -e $errRate -O $minOverlap $primer1 $primer2 -o $outFile $inFile > /dev/null");

    # GET FAILED READS AND REVERSE-COMPLEMENT
    my $numPass = 0;
    my $db = new iTagger::FastqDb::Tiny($inFile);
    my $ok = new iTagger::FastqDb::Tiny($outFile);
    my $tmpInFile = "$this->{tmpdir}/cutadapt.in.fastq";
    open(my $fail, '>', $tmpInFile) or confess($!);
    my $okSeq = $ok->next_seq;
    while ( my $seq = $db->next_seq )
    {
        if ( defined($okSeq) and $seq->id eq $okSeq->id )
        {
            ++$numPass;
            $okSeq = $ok->next_seq;
        } else
        {
            $seq->revcomp;
            print $fail $seq->output;
        }
    }
    close($fail);

    # RUN CUTADAPT (ROUND 2, WITH REVERSE COMPLEMENTED READS)
    my $tmpOutFile = "$this->{tmpdir}/cutadapt.out.fastq";
    $this->run("cutadapt --trimmed-only -e $errRate -O $minOverlap $primer1 $primer2 -o $tmpOutFile $tmpInFile > /dev/null");

    # GET FAILED READS
    my $numFail = 0;
    $db = new iTagger::FastqDb::Tiny($tmpInFile);
    $ok = new iTagger::FastqDb::Tiny($tmpOutFile);
    my $failOut = $failFile ? "| bzip2 --best > $failFile" : '>/dev/null';
    open($fail, $failOut) or confess("Error writing cutadapt failfail: $!");
    open(my $pass, '>>', $outFile) or confess("Error appending cutadapt outfile: $!");
    $okSeq = $ok->next_seq;
    while ( my $seq = $db->next_seq )
    {
        if ( defined($okSeq) and $seq->id eq $okSeq->id )
        {
            ++$numPass;
            print $pass $seq->output;
            $okSeq = $ok->next_seq;
        } else
        {
            ++$numFail;
            $seq->revcomp; # REVERT TO ORIGINAL (UNNECESSARY BUT MAY AVOID CONFUSION)
            print $fail $seq->output;
        }
    }
    close($fail);
    close($pass);
    unlink($tmpInFile,$tmpOutFile);
    return ($numPass, $numFail);
}

=item _cutadaptPaired

Run cutadapt on paired sequences.  A read-pair is filtered if both primers are not found.

=cut

sub _cutadaptPaired
{
    my ($this, $inFile, $outFile, $failFile) = @_;
    return unless exists($this->{config}->{CUTADAPT});
    confess("Missing args") unless $inFile and $outFile;

    # PARAMETERS
    my $params = $this->{config}->{CUTADAPT};
    my $primer1 = $params->{PRIMER1};
    my $primer2 = $params->{PRIMER2};
    confess("CUTADAPT: PRIMER1 and/or PRIMER2 required") unless $primer1 or $primer2;
    foreach (qw/ERROR_RATE MIN_OVERLAP/)
    {
        confess("Missing cutadapt parameter, $_") unless exists($params->{$_});
    }
    my $errRate = $params->{ERROR_RATE};
    my $minOverlap = $params->{MIN_OVERLAP};

    # TEMP FILES
    my @tmpFiles;
    push @tmpFiles, my $inFile1 = "$this->{tmpdir}/cutadapt.in1.fastq";
    push @tmpFiles, my $inFile2 = "$this->{tmpdir}/cutadapt.in2.fastq";
    push @tmpFiles, my $outFile1 = "$this->{tmpdir}/cutadapt.out1.fastq";
    push @tmpFiles, my $outFile2 = "$this->{tmpdir}/cutadapt.out2.fastq";

    # SPLIT INFILE
    my $db = new iTagger::FastqDb::Tiny($inFile, { paired=>1 });
    open(my $fh1, '>', $inFile1) or confess($!);
    open(my $fh2, '>', $inFile2) or confess($!);
    while ( my $pair = $db->next_pair )
    {
        print $fh1 $pair->[0]->output;
        print $fh2 $pair->[1]->output;
    }
    close($fh1);
    close($fh2);

    # RUN CUTADAPT
    if ( $primer1 )
    {
        $this->run("cutadapt --trimmed-only -e $errRate -O $minOverlap -g $primer1 -o $outFile1 $inFile1 > /dev/null");
    } else
    {
        symlink($inFile1, $outFile1) or logconfess($!);
    }
    if ( $primer2 )
    {
        $this->run("cutadapt --trimmed-only -e $errRate -O $minOverlap -g $primer2 -o $outFile2 $inFile2 > /dev/null");
    } else
    {
        symlink($inFile2, $outFile2) or logconfess($!);
    }

    # MERGE CUTADAPT OUTPUT -- PAIR FAILS UNLESS BOTH PRIMERS FOUND
    open(my $out, '>', $outFile) or confess($!);
    my $fail = undef;
    if ( $failFile )
    {
        open($fail, "| bzip2 --best > $failFile") or confess($!);
    }
    my $in1 = new iTagger::FastqDb::Tiny($inFile1);
    my $in2 = new iTagger::FastqDb::Tiny($inFile2);
    my $trim1 = new iTagger::FastqDb::Tiny($outFile1);
    my $trim2 = new iTagger::FastqDb::Tiny($outFile2);
    my $r1 = $trim1->next_seq;
    my $r2 = $trim2->next_seq;
    my $numAdapter = my $numLengthFiltered = 0;
    while ( my $untrimmed1 = $in1->next_seq )
    {
        my $untrimmed2 = $in2->next_seq;
        my $id = $untrimmed1->base;
        if ( $r1 and $r2 and $id eq $r1->base and $id eq $r2->base )
        {
            ++$numAdapter;
            if ( $r1->filtered or $r2->filtered )
            {
                ++$numLengthFiltered;
                if ( $failFile )
                {
                    print $fail $untrimmed1->output, $untrimmed2->output;
                }
            } else
            {
                print $out $r1->output, $r2->output;
            }
            $r1 = $trim1->next_seq;
            $r2 = $trim2->next_seq;
        } else
        {
            if ( $failFile )
            {
                print $fail $untrimmed1->output, $untrimmed2->output;
            }
            $r1 = $trim1->next_seq if $r1 and $r1->base eq $id;
            $r2 = $trim2->next_seq if $r2 and $r2->base eq $id;
        }
    }
    close($out);
    close($fail) if $failFile;
    unlink(@tmpFiles);
    return ($numAdapter, $numLengthFiltered);
}

=item _qualTrim

Trim entire read-pair to desired size by sliding window mean expected error.  A read will not be trimmed beyond it's minimum read length. This method is used by the iterative merging of overlapping paired reads methods. Failed reads (e.g. too short) are written to a separate file, untrimmed.

=cut

sub _qualTrim
{
	my ($this, $inFile, $outFiles, $failFile, $targetLen, $minReadLen) = @_;
	confess("Missing args") unless $inFile and $outFiles and @$outFiles and $failFile and $targetLen;

    # INIT VARS
    $minReadLen = 20 unless defined($minReadLen);
    confess("Invalid target length=$targetLen and mininum read length=$minReadLen") unless $targetLen >= ( 2*$minReadLen );
	my $numPass = my $numFail = 0;

    # OPEN OUTFILES
    open(my $fout, '>', $failFile) or confess($!);
    my ($out, $out1, $out2, $interleaved);
    if ( @$outFiles == 1 )
    {
        $interleaved = 1;
        open($out, '>', $outFiles->[0]) or confess($!);
    } elsif ( @$outFiles == 2 )
    {
        $interleaved = 0;
        open($out1, '>', $outFiles->[0]) or confess($!);
        open($out2, '>', $outFiles->[1]) or confess($!);
    } else
    {
        confess("qualTrim requires either one or two outfiles");
    }

    # TRIM READS
	my $in = new iTagger::FastqDb::Tiny( $inFile, { paired => 1 } );
	while ( my $pair = $in->next_pair )
	{
		my ( $r1, $r2 ) = @$pair;
		my $tot = $r1->len + $r2->len;
		if ( $tot < $targetLen or $r1->len < $minReadLen or $r2->len < $minReadLen or $r1->filtered or $r2->filtered )
		{
            # TOO SHORT TO BEGIN WITH
            ++$numFail;
            print $fout $r1->output, $r2->output;
		}
		elsif ( $tot == $targetLen )
		{
            # HAPPENS TO BE DESIRED LENGTH TO BEGIN WITH
			++$numPass;
            if ( $interleaved )
            {

                print $out $r1->output, $r2->output;
            } else
            {
                print $out1 $r1->output;
                print $out2 $r2->output;
            }
		} else
        {
            # COMMON CASE -- TRIM UNTIL DESIRED LENGTH
            # ITERATIVELY TRIMS END BASE WITH GREATEST MEAN EXPECTED ERROR OF
            # LAST (3') 5 BP
            my $q1  = $r1->qual_arrayref;
            my $q2  = $r2->qual_arrayref;
            my @ee1 = map { 10**( ( $_ * -1 ) / 10 ) } reverse @$q1;
            my @ee2 = map { 10**( ( $_ * -1 ) / 10 ) } reverse @$q2;
            my $e1  = $ee1[0] + $ee1[1] + $ee1[2] + $ee1[3] + $ee1[4];
            my $e2  = $ee2[0] + $ee2[1] + $ee2[2] + $ee2[3] + $ee2[4];
            while ( $tot > $targetLen )
            {
                if ( ( $e1 > $e2 and scalar(@ee1) > $minReadLen ) or scalar(@ee2) == $minReadLen )
                {
                    $e1 -= shift @ee1;
                    $e1 += $ee1[4];
                }
                else
                {
                    $e2 -= shift @ee2;
                    $e2 += $ee2[4];
                }
                --$tot;
            }
            ++$numPass;
            $r1->trim_to( scalar(@ee1) );
            $r2->trim_to( scalar(@ee2) );
            if ( $interleaved )
            {
                print $out $r1->output, $r2->output;
            } else
            {
                print $out1 $r1->output;
                print $out2 $r2->output;
            }
        }
	}
    close($fout);
    if ( $interleaved )
    {
        close($out);
    } else
    {
        close($out1);
        close($out2);
    }
	return ($numPass, $numFail);
}

=item iterativeFlash

Attempt to merge overlapping read-pairs.  Iteratively use flash on increasingly trimmed read-pairs.

=cut

sub iterativeFlash
{
    my ($this, $inFile, $exFile, $ncFile) = @_;
    confess("Missing args") unless $inFile and $exFile and $ncFile;

    # CONFIG
    my $maxMismatch = $this->_set('FLASH', 'MAX_MISMATCH', 0.3);
    my $meanFragLen = $this->_set('AMPLICON', 'LEN_MEAN');
    my $stdevFragLen = $this->_set('AMPLICON', 'LEN_STDEV');
    my $minOverlap = $this->_set('FLASH', 'MIN_OVERLAP', 10);

    # INIT VARS
    my $outDir = "$this->{tmpdir}/flash";
    -d $outDir or mkdir($outDir) or confess($!);
    $inFile = rel2abs($inFile);
    $exFile = rel2abs($exFile);
    $ncFile = rel2abs($ncFile);
    unlink($exFile) if -e $exFile;
    my $tmpInFile = "$this->{tmpdir}/flash.input.fastq";
    my $tmpTrimmedFile = "$this->{tmpdir}/flash.trimmed.fastq";
    my $tmpFailFile = "$this->{tmpdir}/flash.trimFail.fastq";
    symlink($inFile, $tmpInFile) or confess($!);
    open(my $out, "| bzip2 --best > $exFile") or confess($!);
    my $numEx = 0;

    # MERGE, STARTING WITH LONG (MINIMALLY TRIMMED) SEQS
    for ( my $numStdDev = 3; $numStdDev >= -3; $numStdDev -= 0.5 )
    {
        # QUALITY END-TRIM PAIRS
        my $targetLen = int($meanFragLen + $numStdDev * $stdevFragLen + $minOverlap);
        last if $targetLen < MIN_AMPLICON_LEN;
        my ($numQualTrimPass, $numQualTrimFail) = $this->_qualTrim($tmpInFile, [$tmpTrimmedFile], $tmpFailFile, $targetLen, $minOverlap);

        # MERGE WITH FLASH
        my $readLen = int($targetLen/2);
        if ( -s $tmpTrimmedFile )
        {
            $this->run("flash -I -p 33 -d $outDir -m $minOverlap -M $readLen -x $maxMismatch -r $readLen -f $meanFragLen -s $stdevFragLen -t 1 $tmpTrimmedFile");
        }
        unlink($tmpInFile, $tmpTrimmedFile);

        # FIX IDS AND APPEND EXTENDED/MERGED READS
        my $file = "$outDir/out.extendedFrags.fastq";
        if ( -e $file and -s $file )
        {
            my $db = new iTagger::FastqDb::Tiny($file);
            while ( my $rec = $db->next_seq )
            {
                ++$numEx;
                $rec->pair('');
                print $out $rec->output;
            }
        }

        # PREPARE NOT COMBINED
        if ( -e "$outDir/out.notCombined.fastq" )
        {
            move("$outDir/out.notCombined.fastq", $tmpInFile) or confess($!);
        }
        # APPEND TRIM FAILED (I.E. TOO SHORT)
        $this->run("cat $tmpFailFile >> $tmpInFile") if -s $tmpFailFile;
        
        # CLEANUP
        system("rm -rf $outDir") if -e $outDir;
        last unless -s $tmpInFile;
    }
    close($out);

    # GET FULL-LENGTH NOT-COMBINED READS
    my $numNc = 0;
    if ( -s $tmpInFile )
    {
        # GET IDS OF NOT-COMBINED READS
        my %ids;
        my $db = new iTagger::FastqDb::Tiny($tmpInFile, { paired=>1 });
        while ( my $pair = $db->next_pair )
        {
            my $id = $pair->[0]->base;
            $ids{$id} = undef;
        }

        # GET SEQS
        $db = new iTagger::FastqDb::Tiny($inFile, { paired=>1 });
        open(my $out, "| bzip2 --best > $ncFile") or confess($!);
        while ( my $pair = $db->next_pair )
        {
            my $id = $pair->[0]->base;
            next unless exists($ids{$id});
            delete($ids{$id});
            print $out $pair->[0]->output, $pair->[1]->output;
            ++$numNc;
        }
        close($out);
    }
    elsif ( -e $ncFile ) { unlink($ncFile) }
    unlink($tmpInFile);
    return ($numEx, $numNc);
}

=item iterativePandaseq

Merge read-pairs using pandaseq on increasingly trimmed read-pairs.  Pandaseq is an alternative to Flash and the method here is similar.  We don't use the primer trimming capabilities of pandaseq so we can report the number of seqs with/without both primers and we wish for not-combined reads to be primer-trimmed as well (Pandaseq filters reads when primers not found).

=cut

sub iterativePandaseq
{
    my ($this, $inFile, $exFile, $ncFile, $logFile) = @_;
    confess("Missing args") unless $inFile and $exFile and $ncFile and $logFile;

    # CONFIG
    my $threshold = $this->_set('PANDASEQ', 'THRESHOLD', 0.5);
    my $minOverlap = $this->_set('PANDSEQ', 'MIN_OVERLAP', 10);
    my $meanFragLen = $this->_set('AMPLICON', 'LEN_MEAN');
    my $stdevFragLen = $this->_set('AMPLICON', 'LEN_STDEV');

    # INIT VARS
    my $minLen = int($meanFragLen - 3 * $stdevFragLen);
    my $maxLen = int($meanFragLen + 3 * $stdevFragLen);

    # INIT FILES
    my $tmpInFile = "$this->{tmpdir}/pandaseq.in.fastq";
    my $tmpFailFile = "$this->{tmpdir}/pandaseq.trim.fail.fastq";
    my $tmpFile1 = "$this->{tmpdir}/pandaseq.1.fastq";
    my $tmpFile2 = "$this->{tmpdir}/pandaseq.2.fastq";
    my $tmpLogFile = "$this->{tmpdir}/pandaseq.log";
    unlink($exFile) if -e $exFile;
    open(my $log, "| bzip2 --best > $logFile") or confess($!);
    open(my $exFh, "| bzip2 --best >  $exFile") or confess($!);
    symlink($inFile, $tmpInFile) or confess($!);

    # MERGE
    my $numEx = 0;
    for ( my $numStdDev = 3; $numStdDev >= -3; $numStdDev -= 0.5 )
    {
        # QUALITY END-TRIM PAIRS
        my $targetLen = int($meanFragLen + $numStdDev * $stdevFragLen + $minOverlap);
        last if $targetLen < MIN_AMPLICON_LEN;

        my ($numQualTrimPass, $numQualTrimFail) = $this->_qualTrim($tmpInFile, [$tmpFile1, $tmpFile2], $tmpFailFile, $targetLen, $minOverlap);
        unlink($tmpInFile);
        next unless -s $tmpFile1;

        # PANDASEQ
        my $tmpExFile = "$this->{tmpdir}/pandaseq.ex.fastq";
        $this->run("pandaseq -t $threshold -f $tmpFile1 -r $tmpFile2 -F -k 20 -T 1 -w $tmpExFile -g $tmpLogFile -l $minLen -L $maxLen -N -o $minOverlap");

        # APPEND LOG
        open(my $in, '<', $tmpLogFile) or confess($!);
        while (<$in>) { print $log $_ }
        close($in);
        unlink($tmpLogFile);

        # EXTENDED READS: REFORMAT IDS, AND APPEND OUTFILE
        my $db = new iTagger::FastqDb::Tiny($tmpExFile, {paired=>0}) or confess($!);
        my $tmpNumEx = 0; # num this round only (used for development)
        my %ex = ();
        while ( my $rec = $db->next_seq )
        {
            ++$tmpNumEx;
            my $id = $rec->base;
            confess("Unable to parse id: $id") unless $id =~ /^(\S+)\:([ATCG]+)/;
            $rec->base($id=$1);
            $rec->barcode($2);
            $ex{$id} = undef;
            print $exFh $rec->output;
        }
        unlink($tmpExFile);
        $numEx += $tmpNumEx;

        # PREPARE NOT-COMBINED READS
        open(my $out, '>', $tmpInFile) or confess($!);
        my $db1 = new iTagger::FastqDb::Tiny($tmpFile1);
        my $db2 = new iTagger::FastqDb::Tiny($tmpFile2);
        my $tmpNumNc = 0;
        while ( my $read1 = $db1->next_seq )
        {
            my $id = $read1->base;
            my $read2 = $db2->next_seq;
            unless ( defined($read2) )
            {
                my $name = $this->{name};
                confess("Pandaseq: LIB $name missing read 2 for ".$read1->base);
            }
            unless ( $read2->base eq $id )
            {
                my $name = $this->{name};
                confess("Pandaseq: LIB $name, bad not-combined pair: $id, ".$read2->base);
            }
            if ( exists($ex{$id}) )
            {
                delete($ex{$id});
            } else
            {
                print $out $read1->output, $read2->output;
                ++$tmpNumNc;
            }
        }
        close($out);
        # APPEND QUAL FILTERED (I.E. TOO SHORT)
        $this->run("cat $tmpFailFile >> $tmpInFile");
        unlink($tmpFailFile) or confess($!);
    }
    close($log);
    close($exFh);
    unlink($tmpFile1, $tmpFile2);

    # GET FULL-LENGTH NOT-COMBINED SEQS
    # GET IDS
    my $db = new iTagger::FastqDb::Tiny($tmpInFile, { paired=>1 });
    my %nc = ();
    my $numNc = 0;
    while ( my $pair = $db->next_pair )
    {
        my $id = $pair->[0]->base;
        $nc{$id} = undef;
        ++$numNc;
    }
    unlink($tmpInFile);
    # GET SEQS
    open(my $out, "| bzip2 --best > $ncFile") or confess($!);
    $db = new iTagger::FastqDb::Tiny($inFile, { paired=>1 });
    while ( my $pair = $db->next_pair )
    {
        my $id = $pair->[0]->base;
        next unless exists($nc{$id});
        print $out $pair->[0]->output, $pair->[1]->output;
        delete($nc{$id});
    }
    close($out);
    return ($numEx, $numNc);
}

=item lengthFilter

Optionally, the sequence amplicon may be hard-trimmed at 5' and/or 3' ends.  Filter short sequences.

=cut

sub lengthFilter
{
    my ($this, $inFile, $outFile) = @_;
    confess("lengthFilter: infile, outfile required") unless $inFile and $outFile;
    my $trim5 = $this->_set('AMPLICON', 'TRIM5', 0);
    my $trim3 = $this->_set('AMPLICON', 'TRIM3', 0);
    my $minFragLen = $this->_set('AMPLICON', 'LEN_MIN', 0); # undocumented option for testing only
    my $maxFragLen = $this->_set('AMPLICON', 'LEN_MAX', 0); # undocumented option for testing only
    my $meanFragLen = $this->_set('AMPLICON', 'LEN_MEAN');
    my $stdevFragLen = $this->_set('AMPLICON', 'LEN_STDEV');
    my $minLen = MIN_FINAL_LEN;
    if ( $minFragLen )
    {
        $minLen = $minFragLen;
    } else
    {
        $minLen = $meanFragLen - TRIM_NUM_STDEVS * $stdevFragLen;
    }
    $minLen = MIN_FINAL_LEN if $minLen < MIN_FINAL_LEN;
    my $maxLen = $maxFragLen ? $maxFragLen : ( $meanFragLen + TRIM_NUM_STDEVS * $stdevFragLen );

    my $numPass = my $numFail = 0;
    my $db = new iTagger::FastqDb::Tiny($inFile);
    open(my $out, '>', $outFile) or confess($!);
    while ( my $seq = $db->next_seq )
    {
        $seq->trim5_bp($trim5) if $trim5;
        $seq->trim3_bp($trim3) if $trim3;
        if ( $seq->filtered or $seq->len < $minLen or $seq->len > $maxLen )
        {
            ++$numFail;
        } else
        {
            ++$numPass;
            print $out $seq->output;
        }
    }
    close($out);
    return ($numPass, $numFail);
}

=item stitch

For read-pairs that could not be merged/extended, trim read1 and read2 to specified lengths (filtering if either too short), reverse-complement read2, and concatenate.  Optionally, a spacer sequence (e.g. 'NNNNNN') may be inserted between the read-ends.  Returns number of read-pairs successfully stitched together.

=cut

sub stitch
{
    my ($this, $inFile, $outFile) = @_;
    confess("Missing args") unless $inFile and $outFile;
    return unless exists($this->{config}->{STITCH});

    # CONFIG
    confess("STITCH READLEN1 required") unless exists($this->{config}->{STITCH}->{READLEN1});
    confess("STITCH READLEN2 required") unless exists($this->{config}->{STITCH}->{READLEN2});
    my $len1 = $this->{config}->{STITCH}->{READLEN1};
    confess("STITCH: Invalid READLEN1 = $len1") unless $len1 >= 20;
    my $len2 = $this->{config}->{STITCH}->{READLEN2};
    confess("STITCH: Invalid READLEN2 = $len2") unless $len2 >= 20;
    my $spacer = exists($this->{config}->{STITCH}->{SPACER}) ? $this->{config}->{STITCH}->{SPACER} : '';
    my $spacerQ = '!' x length($spacer);

    # TRIM
    my $db = new iTagger::FastqDb::Tiny($inFile, { paired => 1 });
    open(my $out, "| bzip2 --best > $outFile") or confess($!);
    my $pass = 0;
    while ( my $pair = $db->next_pair )
    {
        my ($r1, $r2) = @$pair;
        if ( $r1->len >= $len1 and $r2->len >= $len2 )
        {
            ++$pass;
            $r1->trim_to($len1);
            $r2->trim_to($len2);
            $r2->revcomp;
            $r1->pair('');
            my $stitched = new iTagger::FastqDb::Fastq($r1->header, $r1->seq.$spacer.$r2->seq, $r1->qual.$spacerQ.$r2->qual);
            print $out $stitched->output;
        }
    }
    close($out);
    return $pass;
}

=item qualReport

Generate report on base qualities (for unpaired sequences), using Fastx toolkit's fastx_quality_stats.py, and generate plots.

=cut

sub qualReport
{
    my ($this, $inFile, $tsvFile, $pdfFile, $label) = @_;
    confess("Missing args") unless $inFile and $tsvFile and $pdfFile and $label;
    $this->run("bzcat $inFile | fastx_quality_stats -Q33 -o $tsvFile");
    if ( -s $tsvFile )
    {
        $this->run("itaggerQscorePlots.pl --infile_1 $tsvFile --name $label --pdf $pdfFile --display 1 --single");
    }
}

=item qualReportPaired

Generate report on base qualities for paired sequences, using Fastx toolkit's fastx_quality_stats.py, and generate plots.

=cut

sub qualReportPaired
{
    my ($this, $inFile1, $inFile2, $tsvFile1, $tsvFile2, $pdfFile, $label) = @_;
    confess("Missing args") unless $inFile1 and $inFile2 and $tsvFile1 and $tsvFile2 and $pdfFile and $label;
    $this->run("fastx_quality_stats -Q33 -i $inFile1 -o $tsvFile1");
    $this->run("fastx_quality_stats -Q33 -i $inFile2 -o $tsvFile2");
    $this->run("itaggerQscorePlots.pl --infile_1 $tsvFile1 --infile_2 $tsvFile2 --name $label --pdf $pdfFile --display 1 --paired");
}

=item expErrFilter

Discard sequences with expected error exceeding threshold.

=cut

sub expErrFilter
{
    my ($this, $inFile, $outFile) = @_;

    # DEFINE STRINGENCY
    my %params = ();
    my $ee = $this->_set('AMPQC', 'EXP_ERR', 0);
    my $ee_per_kb = $this->_set('AMPQC', 'EXP_ERR_PER_KB', 20);
    if ( $ee ) { $params{exp_err} = $ee }
    else { $params{exp_err_per_kb} = $ee_per_kb }

    # FILTER
    open(my $out, '>', $outFile) or confess($!);
    my $in = new iTagger::FastqDb::Tiny($inFile, \%params);
    while ( my $rec = $in->next_seq ) { print $out $rec->output }
    $in->close_file;
    close($out);
    return $in->summary->{pass};
}

=item seqObs

Generate report on base qualitites, dereplicate sequences, and output in seq-obs tabular format.

=cut

sub seqObs
{
    my ($this, $inFile, $tsvFile, $pdfFile, $label, $seqObsFile) = @_;
    confess("Missing args") unless $inFile and $tsvFile and $pdfFile and $label and $seqObsFile;

    # READ SEQS
    open(my $out, "| fastx_quality_stats -Q33 -o $tsvFile") or confess($!);
    my %seqObs = (); # there's plenty of RAM for 1 lane
    my $totSize =0;
    my $in = new iTagger::FastqDb::Tiny($inFile);
    while ( my $rec = $in->next_seq )
    {
        print $out $rec->output;
        ++$totSize;
        my $seq = uc($rec->seq);
        if ( exists($seqObs{$seq}) ) { $seqObs{$seq} += 1 } 
        else { $seqObs{$seq} = 1 }
    }
    close($out);
    if ( -s $tsvFile )
    {
        $this->run("itaggerQscorePlots.pl --infile_1 $tsvFile --name $label --pdf $pdfFile --display 1 --single");
    }

    # SEQOBS
    open($out, "| bzip2 --best > $seqObsFile") or confess($!);
    print $out "#lib=$this->{name};size=$totSize\n";
    my @seqs = sort { $a cmp $b } keys %seqObs;
    my $numUniq = 0;
    while ( my $seq = shift @seqs )
    {
        ++$numUniq;
        my $size = $seqObs{$seq};
        delete $seqObs{$seq};
        print $out "$seq\t", $size, "\n";
    }
    close($out);
    return $numUniq;
}

=item $this->run

Run external executable, logconfess on nonzero exit status.

=cut

sub run
{
    my ($this, $cmd) = @_;
    my $output = `$cmd 2>&1`;
    confess("FAILURE RUNNING $cmd: $output") unless $? == 0;
}

=item $this->_set

Return value from config or default value.  If default not defined, it must be defined in the config, otherwise logconfess with error message.

=cut

sub _set
{
    my ($this, $section, $key, $default) = @_;
    confess("Missing args") unless $section and $key;
    my $config = $this->{config};
    if ( !exists($config->{$section}) )
    {
        $config->{$section} = {};
        return $config->{$section}->{$key} = $default if defined($default);
        confess("Config file missing section, $section");
    } elsif ( exists($config->{$section}->{$key}) )
    {
        return $config->{$section}->{$key}
    } elsif ( defined($default) )
    {
        return $config->{$section}->{$key} = $default;
    } else
    {
        confess("Config missing required $section/$key");
    }
}

1;

=back

=head1 AUTHORS

Edward Kirton, Julien Tremblay

=head1 COPYRIGHT

Copyright (c) 2013 US DOE Joint Genome Institute.  Use freely under the same license as Perl itself.  Refer to duk and flash documentation for their own copyright/license information.

=cut
