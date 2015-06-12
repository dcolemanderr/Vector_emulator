=pod

=head1 NAME

FastqDb - Simple object for handling Fastq files

=head1 SYNOPSIS

    # simple script to filter low-qual read pairs
    use iTagger::FastqDb::Tiny;
    my $file = shift @ARGV;
    my $db = new FastqFile::Tiny($file, {paired=>1, exp_err=>0.5});
    while (my $pair = $db->next_pair) {
        my ($r1,$r2) = @$pair;
        print $r1-output, $r2->output;
    }
    my ($pass, $fail) = $db->pass_fail;
    print "pass = $pass; fail = $fail\n":

=head1 DESCRIPTION

This module provides an iterator for Fastq files, plus basic QC functions

=head1 METHODS

=over 5

=cut

package iTagger::FastqDb;

use warnings;
use strict;
use Carp;
use IO::File;
use File::Spec;
use constant {SAMPLE_SIZE => 10_000};
use Env qw/TMPDIR/;
use iTagger::FastqDb::Fastq;

our $VERSION = 1.0;

=item new $file

Constructor accepts a path to an existing Fastq file (otherwise reads from standard input).

The quality format will be automatically detected; the output always uses sanger-scaled quality scores.

=cut

sub new
{
	my ($class, $file, $params) = @_;
    if (defined($file))
    {
        if ($file eq '-')
        {
            $file=undef;
        } elsif ($file !~ /\|\s*$/)
        {
            $file = readlink($file) if -l $file;
            $file = File::Spec->rel2abs($file);
            confess("Input file, $file, not found\n") unless -e $file;
        }
	}
	my $this = {
		file       => $file,
		line       => 0,
		paired     => undef,
		counter    => { pass  => 0, error => 0 },
		buffer => [],           # buffer for sampled reads
		sample => {
			readlen => undef,
			size    => undef,
			file    => undef
		}};
	bless $this, $class;
	$this->_init_qc_params($params);
	$this->_open_file;
	$this->_init;
	return $this;
}

=item DESTROY

Close file and cleanup tmp file before going away.

=cut

sub DESTROY
{
	my $this = shift;
	$this->close_file;
	unlink($this->{sample}->{file}) if exists($this->{sample}->{file}) and defined($this->{sample}->{file}) and -e $this->{sample}->{file};
}

=item _open_file

Open file, possibly compressed, or stdin.

=cut

sub _open_file
{
	my $this = shift;
	$this->close_file if exists($this->{fhr});
	my $file = $this->{file};
	my $fh   = new IO::File;
	if (!$file)
	{
		$fh = *STDIN;
	} elsif ($file =~ /\.gz$/)
	{
		open($fh, "gunzip -c $file |") or confess("Unable to open gzipped infile, $file: $!\n");
	} elsif ($file =~ /\.bz2$/)
	{
		open($fh, "bzcat $file |") or confess("Unable to open bzipped infile, $file: $!\n");
	} elsif (-e $file)
	{
		open($fh, '<', $file) or confess("Unable to open infile, $file: $!\n");
	} else
	{
		open($fh, $file) or confess("Unable to open infile, $file: $!\n");
	}
	$this->{fhr}  = \$fh;
}

=item close_file

Close file

=cut

sub close_file
{
	my $this = shift;
	return unless exists($this->{fhr});
	if (defined($this->{fhr}))
	{
		my $fh = ${$this->{fhr}};
		close($fh) if $this->{file};
	}
	delete $this->{fhr};
}

=item _init_qc_params

Init QC params hash with default values.

=cut

sub _init_qc_params
{
	my ($this, $params) = @_;
	$params = {} unless defined($params);
    # DEFAULT PARAMS
	my %qc_params =
    (
        no_qc          => 0,
        ignore_pairing => 0,
		informat       => undef,
		qual_to_sanger => undef,
		roche_mid_len  => undef,
        trim3_bp       => 0,        # trim this many bp from 3' end
        trim5_bp       => 0,        # trim this many bp from 5' end
		trim3          => 0,        # trim uncalled bases from 3' end
		min_len        => 1,        # minimum length of good sequences
		max_n          => undef,    # maxmimum number of Ns allowed in a good sequence
		low            => undef,    # low complexity threshold (0<x<1)
		low_q          => undef,    # Q-score considered low for Q-count filter
		max_num_low_q  => undef,    # Q-count filters maximum number of low-qual bases
        trim3_exp_err  => undef,    # trim 3' once 5' bp exceed max err
        exp_err        => undef,    # filter read if expected error exceed threshold (absolute)
        exp_err_per_kb => undef,    # filter read if expected error exceed threshold (relative)
        trim_to        => undef,    # trim read to exactly this many bp
		learn_barcodes     => 0,        # learn barcodes from sample if true
		barcodes_file      => undef,    # path to barcode sequences in Fasta or Tsv format
		barcodes_seq_col   => 0,        # column of Tsv file containing sequences (first col = 1)
		barcodes_label_col => undef,    # column of Tsv file containing labels (first col = 1)
		barcodes           => undef,    # { barcode => counter }
		barcode_variants   => undef,    # { 1bp error barcode => real barcode }
		barcode_labels     => undef,    # { barcode => label }
	);
	foreach my $key (keys %$params)
    {
        if ( $key eq 'paired' ) { $this->{paired} = $params->{$key} ? 1:0 }
        else { $qc_params{$key} = $params->{$key} }
    }
	$this->{qc_params} = \%qc_params;
}

=item paired

Returns true if paired, false otherwise.

=cut

sub paired { return shift->{paired} }

########################################
## SAMPLE DATA FOR PARAMETER TRAINING ##
########################################
#
# The file is randomly sampled for parameter training purposes:
#
# - whether reads are paired or not (to output empty reads if paired)
# - whether input is in older FastqIllumina format (and must be converted to Sanger-scale)
# - the molecular barcodes used (if barcoded)

=item sample_file

Returns path to the sampled data in Fastq format or undef if input was not sampled.

=cut

sub sample_file
{
	my $this = shift;
	my $file = $this->{sample}->{file};
	if (defined($file) and -e $file and -s $file)
	{
		return $file;
	} else
	{
		return $this->_write_sample;
	}
}

=item sample_arrayref

Returns a reference to an array of the buffered sample.

=cut

sub sample_arrayref
{
	my $this = shift;
	my @sample = @{$this->{buffer}};
    return \@sample;
}

=item sample_size

Returns the number of sequences sampled.

=cut

sub sample_size { return shift->{sample}->{size} }

=item _init

Buffer SAMPLE_SIZE good reads for parameter training, determine quality score encoding, if paired/unpaired, and possibly learn molecular barcodes.

=cut

sub _init
{
	my $this = shift;
	$this->{buffer} = [];
	my $counter          = 0;       # number of good reads
    my $appears_pacbio   = 0;
	my $appears_sanger   = 0;
	my $appears_illumina = 0;
	$this->{paired}   = undef;
    if ( $this->{qc_params}->{ignore_pairing} )
    {
        $this->{paired} = 0;
    }

	# COPY QC PARAMS TO SAMPLE WITHOUT QC AS QUAL-ENCODING NEEDS TO BE DETECTED FIRST
	my %qc_params = %{$this->{qc_params}};
	$this->{qc_params}->{no_qc} = 1;

	# READ FILE, DETERMINE PAIRING AND QUALITY-SCORE ENCODING (NO QC ON READS)
    my @buffer = ();
	while (my $read = $this->_next_seq)
	{
        push @buffer, $read;
		my @qual = split(//, $read->qual);
		foreach my $q (@qual)
		{
			my $ascii = ord($q);
            if    ($ascii < 33)  { confess("Invalid quality character\n") }
            elsif ($ascii > 113) { confess("Invalid quality character\n") }
            elsif ($ascii > 104) { ++$appears_pacbio }
			elsif ($ascii < 64)  { ++$appears_sanger }
			elsif ($ascii > 74)  { ++$appears_illumina }
		}
        unless ($this->{qc_params}->{ignore_pairing})
        {
            if (defined($this->{paired}))
            {
                my $read_paired = $read->pair ? 1:0;
                unless ( $read_paired == $this->{paired} )
                {
                    confess("Infile contains a mix of paired and unpaired reads\n");
                }
            } else
            {
                $this->{paired} = $read->pair ? 1 : 0;
            }
        }
		last if ++$counter >= SAMPLE_SIZE;
	}
    unless ($counter)
    {
        $this->{sample}->{size} = $this->{sample}->{readlen} = $this->{sample}->{len_stddev} = undef;
        unlink($this->{sample}->{file}) if defined($this->{sample}->{file}) and -e $this->{sample}->{file};
	    $this->{buffer} = [];
        $this->close_file;
        return;
    }

	# RESET QC PARAMS
	$this->{qc_params} = \%qc_params;

	# SAVE RESULT
	my $convert_qual = 0;
    if ( $appears_pacbio ) { $convert_qual = 2 }
    elsif ( $appears_sanger < $appears_illumina ) { $convert_qual = 1 }
	$this->{qc_params}->{qual_to_sanger} = $convert_qual;
	$this->{sample}->{size}              = $counter;

	# PROCESS BUFFERED READS (CONVERT QUAL, TRIM, AND FILTER)
	for (my $i=0; $i<=$#buffer; $i++)
	{
		$buffer[$i]->convert_qual_to_sanger if $convert_qual;
		$buffer[$i]->qc($this->{qc_params});
	}
	$this->_init_barcodes(\@buffer);

	# COUNT NUMBER OF GOOD READS IN BUFFER AND DETERMINE LENGTH DISTRIBUTION
	my @read_lengths = ();
	$counter = 0;    # NUMBER OF GOOD READS
	if ($this->{paired})
	{
		for (my $i = 0; $i <= $#buffer; $i += 2)
		{
			next if $buffer[$i]->filtered;
			next if $buffer[ $i + 1 ]->filtered;
			push @read_lengths, $buffer[$i]->len;
			push @read_lengths, $buffer[ $i + 1 ]->len;
			$counter += 2;
		}
	} else
	{
		for (my $i = 0; $i <= $#buffer; $i++)
		{
			next if $buffer[$i]->filtered;
			push @read_lengths, $buffer[$i]->len;
			++$counter;
		}
	}

	# READ MORE SEQUENCES AS SOME WERE FILTERED
	if ($counter < SAMPLE_SIZE)
	{
		if ($this->{paired}) 
		{
			my $read1;
			my $read2;
			while ($read1 = $this->_next_seq)
			{
				$read2 = $this->_next_seq;
				confess("Unpaired read: " . $read1->id . "\n") unless defined($read2);
				push @buffer, ($read1, $read2);
				unless ($read1->filtered or $read2->filtered)
				{
					$counter += 2;
                    push @read_lengths, $read1->len;
                    push @read_lengths, $read2->len;
					last if $counter >= SAMPLE_SIZE;
				}
			}
		} else
		{
			my $read;
			while ($read = $this->_next_seq)
			{
				push @buffer, $read;
				unless ($read->filtered)
				{
                    push @read_lengths, $read->len;
					last if ++$counter >= SAMPLE_SIZE;
				}
			}
		}
	}
    unless ($counter)
    {
        $this->{sample}->{size} = $this->{sample}->{readlen} = $this->{sample}->{len_stddev} = undef;
        unlink($this->{sample}->{file}) if defined($this->{sample}->{file}) and -e $this->{sample}->{file};
	    $this->{buffer} = [];
        $this->close_file;
        return;
    }
    confess("Invalid read lengths\n") unless @read_lengths;
	$this->{sample}->{size} = $counter;
	my ($mean, $stddev) = _mean_stddev(\@read_lengths);
	$this->{sample}->{readlen}    = int($mean + 0.5);
	$this->{sample}->{len_stddev} = int($stddev * 10 + 0.5) / 10;
	unlink($this->{sample}->{file}) if defined($this->{sample}->{file}) and -e $this->{sample}->{file};
    $this->{buffer} = \@buffer;
}

=item _write_sample

Write sample to disk

=cut

sub _write_sample
{
	my $this = shift;
	return undef unless defined($this->{sample}->{size}) and $this->{sample}->{size};
	my $outfile = $this->{sample}->{file} = $TMPDIR ? "$TMPDIR/$$.fastq" : "./$$.fastq";
	open(OUT, ">$outfile") or confess($!);
	for (my $i = 0; $i < $#{$this->{buffer}}; $i++)
	{
        print OUT $this->{buffer}->[$i]->output unless $this->{buffer}->[$i]->filtered;
	}
	close(OUT);
	return $outfile;
}

=item next_seq

An iterator for the Fastq database.  Returns the next valid read object.
If the reads are paired and one read is filtered, the entire pair is filtered. The &next_pair method is also useful for returning both pairs at once.

=cut

sub next_seq
{
	my ($this) = @_;
    return $this->_next_single unless $this->{paired};
    if ( exists( $this->{next_read2} ) )
    {
        my $read2 = $this->{next_read2};
        delete($this->{next_read2});
        return $read2;
    }
    else
    {
        my $pair = $this->next_pair;
        return undef unless defined($pair);
        my ( $read1, $read2 ) = @$pair;
        $this->{next_read2} = $read2;
        return $read1;
    }
}

=item _next_single

Get next single read from buffer or input.

=cut

sub _next_single
{
	my ($this) = @_;
	my $rec;
	while (1)
	{
        $rec = $this->_next_seq;
		last unless defined($rec);
		if ($rec->filtered)
		{
			$this->_tally($rec->filtered);
		} else
        {
            $this->{counter}->{pass} += 1;
            $this->{barcodes}->{$rec->barcode} += 1 if defined($this->{barcodes});
            last;
        }
	}
	return $rec;
}

=item next_pair

Returns a pair of good reads.  If one read was filtered, the entire pair is filtered.

=cut

sub next_pair
{
	my ($this) = @_;
	confess("This method cannot be used with unpaired input\n") unless $this->{paired};
	my $rec1;
	my $rec2;
    do
    {
        do
        {
            # GET READ1
            while ( !$rec1 )
            {
                $rec1 = $this->_next_seq;
                return undef unless defined($rec1);
                if ( $rec1->pair == 2 )
                {
                    $this->_tally('singleton');
                    $rec1 = undef;
                }
            };

            # GET READ2
            $rec2 = $this->_next_seq;
            if ( !defined($rec2) )
            {
                $this->_tally('singleton');
                return undef;
            }
            elsif ( $rec2->pair == 1 )
            {
                $this->_tally('singleton');
                $rec1 = $rec2; 
                $rec2 = undef;
            }
            elsif ( $rec1->base ne $rec2->base )
            {
                $this->_tally('singleton');
                $this->_tally('singleton');
                $rec1 = $rec2 = undef;
            }
        } until ( $rec1 and $rec2 );

        # CHECK IF EITHER READ WAS FILTERED
        my $filtered = 0;
        if ( $rec1->filtered )
        {
            $this->_tally( $rec1->filtered );
            ++$filtered;
        }
        if ( $rec2->filtered )
        {
            $this->_tally( $rec2->filtered );
            ++$filtered;
        }
        if ( $filtered )
        {
            $this->_tally('singleton') if $filtered == 1;
            $rec1 = $rec2 = undef;
        } else
        {
            $this->{counter}->{pass} += 2;
            if ( defined($this->{barcodes}) and $rec1->barcode )
            {
                if ( $rec1->barcode ne $rec2->barcode )
                {
                    if ( exists($this->{barcodes}->{$rec1->barcode}) )
                    {
                        $rec2->barcode($rec1->barcode);
                    } elsif ( exists($this->{barcodes}->{$rec2->barcode}) )
                    {
                        $rec1->barcode($rec2->barcode);
                    } else
                    {
                        # neither barcode good, let's at least make them the same
                        $rec2->barcode($rec1->barcode);
                    }
                }
                $this->{barcodes}->{$rec1->barcode} += 1;
            }
        }
    } until ( $rec1 and $rec2 );
	return [ $rec1, $rec2  ];
}

=item _tally

Update counter for filtered reads

=cut

sub _tally
{
	my ($this, $reason) = @_;
	if (exists($this->{counter}->{$reason}))
	{
		$this->{counter}->{$reason} += 1;
	} else
	{
		$this->{counter}->{$reason} = 1;
	}
}

=item _next_seq

Get next seq from buffer or file

=cut

sub _next_seq
{
    my ($this) = @_;
	my $rec = scalar(@{$this->{buffer}}) ? shift @{$this->{buffer}} : undef;
    return $rec if $rec;
	return undef unless exists($this->{fhr}) and defined($this->{fhr});
	my $fh = ${$this->{fhr}};
	my ($hdr, $seq, $sep, $qual);
	my $error = 0;
	while (<$fh>)
	{
		chomp;
		if (!defined($hdr))
		{
			if (/^@/)
			{
				$hdr   = $_;
				$error = 0;
			} else
			{
				$this->{counter}->{error} += 1;
			}
		} elsif (!defined($seq))
		{
			if ($_ eq '' or /^[ATCGN]+$/i)
			{
				$seq = $_;
			} else
			{
				$hdr   = undef;
				$error = 1;
				$this->{counter}->{error} += 1;
			}
		} elsif (!defined($sep))
		{
			if (/^\+/)
			{
				$sep = $_;
			} else
			{
				$hdr = $seq = undef;
				$error = 1;
				$this->{counter}->{error} += 1;
			}
		} elsif (!defined($qual))
		{
			if (length($_) == length($seq))
			{
				$qual = $_;
				last;
			} else
			{
				$hdr = $seq = $sep = undef;
				$error = 1;
				$this->{counter}->{error} += 1;
			}
		}
	}
	unless (defined($hdr))
	{
		$this->close_file;
		return undef;
	}
	return new iTagger::FastqDb::Fastq($hdr, $seq, $qual, $this->{qc_params}, $this->{barcodes}, $this->{barcode_variants});
}

=item summary

Return summary of good and filtered reads.

=cut

sub summary
{
	my $this    = shift;
    my %summary = ();
    foreach ( keys %{$this->{counter}} )
    {
        $summary{$_} = $this->{counter}->{$_} // 0;
    }
	return \%summary;

}

=item pass_fail

Return number of reads with pass, fail QC.

=cut

sub pass_fail
{
    my $this = shift;
    my $pass = $this->{counter}->{pass};
    my $fail = 0;
    foreach ( %{$this->{counter}} )
    {
        $fail += $this->{counter}->{$_} // 0 unless $_ eq 'pass';
    }
    return ($pass, $fail);
}

=item barcodes

Return a list of barcodes.

=cut

sub barcodes
{
    my $this=shift;
    return keys %{$this->{barcodes}};
}

=item barcode_counts

Returns a hash of barcodes and counts.  Returns undef if reads were not barcoded. Normally this is only run after iterating through all reads in the database.

=cut

sub barcode_counts
{
	my $this = shift;
	return undef unless defined($this->{barcodes});
	my %barcodes = %{$this->{barcodes}};    # deep copy
	return \%barcodes;
}

=item barcode_labels

Returns a hash of barcodes and labels.  Returns undef if reads were not barcoded or labels not defined.

=cut

sub barcode_labels
{
	my $this = shift;
	return undef unless defined($this->{barcode_labels});
	my %barcode_labels = %{$this->{barcode_labels}};    # deep copy
	return \%barcode_labels;
}

=item _init_barcodes

Prepare to handle barcoded reads

=cut

sub _init_barcodes
{
	my ($this, $buffer) = @_;
	return unless $this->{qc_params}->{learn_barcodes} or (exists($this->{qc_params}->{barcodes_file}) and defined($this->{qc_params}->{barcodes_file}));
	if ($this->{qc_params}->{learn_barcodes}) { $this->_learn_barcodes($buffer) }
	else { $this->_load_barcodes }
	# PROCESS BUFFERED SEQS
    for ( my $i=0; $i<=$#$buffer; $i++ )
	{
		$buffer->[$i]->check_barcode($this->{barcodes}, $this->{barcode_variants});
	}
}

=item _load_barcodes

Load the barcodes provided in a tabular or Fasta file

=cut

sub _load_barcodes
{
	my $this = shift;
	my $file = $this->{qc_params}->{barcodes_file} or confess("No barcodes file provided\n");
    my $seq_col   = $this->{qc_params}->{barcodes_seq_col};
	my $label_col = undef;
	if ( exists($this->{qc_params}->{barcodes_label_col}) and defined($this->{qc_params}->{barcodes_label_col}) )
    {
        $label_col = $this->{qc_params}->{barcodes_label_col};
    }
	my $barcode = '';
    my $label = '';
	my %barcodes         = ();
	my %barcode_variants = ();
	my %barcode_labels   = ();
    my @data = ();
	open(IN, "<$file") or confess("Unable to open molecular barcodes file, $file: $!\n");
    while (<IN>)
    {
        chomp;
        push @data, $_;
    }
	close(IN);

	# DETERMINE FORMAT AND PARSE ACCORDINGLY
	if ($data[0] =~ /^>\w+/)
	{
		# FASTA
		while (@data)
		{
			my $line = shift @data;
			if ($line =~ /^>(.+)$/)
			{
				$label = $1;
			} elsif ($line =~ /^[ATCG]+$/i)
			{
				$barcode = uc($line);
				$barcodes{$barcode} = 0;
				my @variants = _generate_barcode_variants($barcode);
				foreach my $variant (@variants)
				{
					$barcode_variants{$variant} = $barcode;
				}
				$barcode_labels{$barcode} = $label;
			} elsif ($line)
			{
				confess("Error reading molecular barcodes file\n");
			}
		}
	} else
	{
		# TABULAR
		while (@data)
		{
			my $line = shift @data;
			my @row = split(/\t/, $line);
			$barcode = uc($row[$seq_col]);
			confess("Invalid barcode, $barcode\n") unless $barcode =~ /^[ATCGN]+$/;
			$barcodes{$barcode} = 0;
			my @variants = _generate_barcode_variants($barcode);
			foreach my $variant (@variants)
			{
				$barcode_variants{$variant} = $barcode;
			}
			if (defined($label_col))
			{
				$barcode_labels{$barcode} = $row[$label_col];
			}
		}
	}
	$this->{barcodes}         = \%barcodes;
	$this->{barcode_variants} = \%barcode_variants;
	$this->{barcode_labels}   = \%barcode_labels;
}

=item _learn_barcodes

Discover barcodes using abundant sequences and assuming all real barcodes have counts within three standard deviations of the mean abundance.

=cut

sub _learn_barcodes
{
	my ($this, $buffer) = @_;
	my %putative_barcodes = ();    # seq => count
	my %barcodes          = ();    # sequence of real barcode or variant => number of times observed
	my %barcode_variants  = ();    # variant barcode sequence => real barcode sequence
	my $barcode;
	my $n;
	my %counts = ();

	# COUNT ALL BARCODE SEQUENCES IN SAMPLE
	for (my $i = 0; $i <= $#$buffer; $i++)
	{
		my $read = $buffer->[$i];
		if ($barcode = $read->barcode)
		{
			if (exists($putative_barcodes{$barcode}))
			{
				$putative_barcodes{$barcode} += 1;
			} else
			{
				$putative_barcodes{$barcode} = 1;
			}
		}
	}

	# SORT PUTATIVE BARCODES IN DESCENDING ORDER OF COUNTS
	foreach $barcode (keys %putative_barcodes)
	{
		$n = $putative_barcodes{$barcode};
		if (exists($counts{$n})) { push @{$counts{$n}}, $barcode }
        else { $counts{$n} = [$barcode] }
	}
	my @counts = sort { $b <=> $a } keys %counts;

	# PICK REAL BARCODES (ASSUMES ABUNDANCE OF ALL REAL BARCODES IS WITHIN 3SD OF MEAN)
	my @real_barcode_counts = ();
	my $mean;
	my $stddev;
	while ($n = shift @counts)
	{
		if (@real_barcode_counts == 0) { 1 } 
        elsif (@real_barcode_counts < 3)
		{
			($mean, $stddev) = _mean_stddev(\@real_barcode_counts);
			last if $n * 100 < $mean;
		} else
		{
			($mean, $stddev) = _mean_stddev(\@real_barcode_counts);
			last if $n < ($mean - 3 * $stddev);
		}

		# REAL BARCODE(S)
		foreach $barcode (@{$counts{$n}})
		{
			next unless exists($putative_barcodes{$barcode});
			push @real_barcode_counts, $n;
			$barcodes{$barcode} = 0;
			# DELETE VARIANTS
			my @variants = _generate_barcode_variants($barcode);
			foreach my $variant (@variants)
			{
                delete($putative_barcodes{$variant}) if (exists($putative_barcodes{$variant}));
				$barcode_variants{$variant} = $barcode;
			}
		}
		delete($counts{$n});
	}

	# GENERATE LABELS HASH (LABELS ARE JUST THE SEQ IN THIS CASE)
	my %labels = ();
	foreach $barcode (keys %barcodes) { $labels{$barcode} = $barcode }

	# SAVE RESULTS
	$this->{barcodes}         = \%barcodes;
	$this->{barcode_variants} = \%barcode_variants;
	$this->{barcode_labels}   = \%labels;

	# WRITE BARCODES TO FILE (TXT/TSV FORMAT)
	if ($this->{qc_params}->{barcodes_file})
	{
		open(OUT, '>' . $this->{qc_params}->{barcodes_file}) or confess("Unable to write barcodes to file, " . $this->{qc_params}->{barcodes_file} . ": $!\n");
		foreach my $barcode (keys %barcodes) { print OUT "$barcode\n" }
		close(OUT);
	}
}

=item _mean_stddev

Return mean and standard deviation of list

=cut

sub _mean_stddev
{
	my $ar  = shift;
	my $n   = scalar(@$ar);
    return (undef, undef) unless $n;
	my $tot = 0;
	foreach my $x (@$ar) { $tot += $x }
	my $mean = $tot / $n;
	my $d    = 0;
	foreach my $x (@$ar) { $d += ($mean - $x)**2 }
	my $stddev = sqrt($d / $n);
	return ($mean, $stddev);
}

=item _generate_barcode_variants

Generate 1bp error variants of the barcode

=cut

sub _generate_barcode_variants
{
	my $seq      = shift;
	my @seq      = split(//, $seq);
	my @variants = ();
	for (my $i = 0; $i <= $#seq; $i++)
	{
		foreach my $nt (qw/A T C G N/)
		{
			next if $nt eq $seq[$i];
			my @variant = @seq;
			$variant[$i] = $nt;
			push @variants, join('', @variant);
		}
	}
	return @variants;
}

1;

=back

=head1 BUGS AND LIMITATIONS

Uses a modified Illumina naming convention: base#barcode/pair

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut
