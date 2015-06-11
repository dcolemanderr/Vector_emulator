=pod

=head1 NAME

iTagger::Fastq - Simple object for Fastq sequence

=head1 SYNOPSIS

    my $rec=new Fastq( $hdr, $id, $base, $barcode, $pair, $seq, $qual, $qc_params, $barcodes, $barcode_variants);
    print $rec->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package iTagger::FastqDb::Fastq;

use warnings;
use strict;
use Carp;
use constant {
	CHARACTERS_PER_LINE    => 80,    # for formatting Fasta/Qual output only
	DEFAULT_MINLEN         => 20,
	DEFAULT_MEANQ          => 20,
	DEFAULT_WINSIZE        => 5,
	DEFAULT_LOW_COMPLEXITY => 0.8,
	DEFAULT_MAXN           => 3};
use iTagger::Stats qw(oddMedian round);

our $VERSION = 1.0;

=item new $hdr $seq $qual $qc_params $barcodes $variants

Initialize new sequence object. If the quality scores use Illumina scaling, the $qual_to_sanger flag *must* be set as the object
assumes and requires sanger-scaling.  Quality encoding method is determined by FastqDb class instead.

=cut

sub new
{
	my ($class, $hdr, $seq, $qual, $qc_params, $barcodes, $barcode_variants) = @_;
	confess("Missing hdr\n")           unless defined($hdr);
	confess("Missing seq for $hdr\n")  unless defined($seq);
	confess("Missing qual for $hdr\n") unless defined($qual);
    confess("Seq and qual don't match for: $hdr\n") unless length($seq) == length($qual);

	# INIT
	my $this = {
		seq      => $seq,     # complete sequence without newlines
		qual     => $qual,    # complete quality string without newlines
		filtered => undef,    # reason why filtered
        qc       => $qc_params
	};
	bless $this, $class;
	$this->_parse_header($hdr);      # populates id, base, pair, barcode values
	$this->convert_qual_to_sanger($this->{qc}->{qual_to_sanger});
	$this->trim_roche_mid($this->{qc}->{roche_mid_len}) if defined($this->{qc}->{roche_mid_len});
	$this->check_barcode($barcodes, $barcode_variants) if defined($barcodes);
	$this->qc;
	return $this;
}

=item convert_qual_to_sanger

Convert the quality string from Illumina to Sanger scaling.

=cut

sub convert_qual_to_sanger
{
	my ($this, $type) = @_;
    unless ( $type )
    {
        1;
    } elsif ( $type == 1 ) # illumina pre-1.8
    {
        $this->{qual} = join('', map { chr(ord($_) - 31) } split(//, $this->{qual}));
    } elsif ( $type == 2 ) # pac bio
    {
        $this->{qual} = join('', map { chr(int((ord($_) - 33)/2)+33) } split(//, $this->{qual}));
    } else
    {
        confess("Invalid convert qual code\n");
    }
}

=item trim_roche_mid

Trim the first x bases which are the 454 MID.

=cut

# TRIM THE ROCHE MOLECULAR ID FROM THE 5' END OF THE SEQUENCE
sub trim_roche_mid
{
	my ($this, $len) = @_;
	return unless defined($len);
	confess("Invalid MID length, $len\n") unless $len > 0;
	if (length($this->{seq}) <= $len)
	{
		$this->{seq} = $this->{qual} = undef;
		return;
	}
	$this->{barcode} = uc(substr($this->{seq}, 0, $len));
	$this->{seq}  = substr($this->{seq},  $len);
	$this->{qual} = substr($this->{qual}, $len);
}

=item check_barcode

Check if barcode is valid, perform 1-base error correction, or filter read.

=cut

sub check_barcode
{
	my ($this, $barcodes, $barcode_variants) = @_;
	return unless defined($barcodes) and defined($barcode_variants);
	my $barcode = $this->barcode;
	return unless $barcode;
	return if exists($barcodes->{$barcode});
	if (exists($barcode_variants->{$barcode}))
	{
		$this->barcode($barcode_variants->{$barcode});
	} else
	{
		$this->del('invalid_barcode');
	}
}

=item _parse_header ($hdr)

Private method used by constructor to parse a header and populate object's ID-related variables.  Doesn't return anything.

=cut

sub _parse_header
{
	my ($this, $hdr) = @_;
	$hdr = '@' . $hdr unless $hdr =~ /^@/;
    $this->{header} = $hdr;
    $this->{id}      = undef;
	$this->{base}    = undef;    # base ID only (e.g. "A"); always defined
	$this->{pair}    = undef;    # pair ID only (e.g. "1"); only defined if paired
	$this->{barcode} = undef;    # barcode sequence (always upper-case); only defined it barcoded
    $this->{control} = 0;

	if ($hdr =~ /^@(\w+_\d+_\d+_\w+_\w+_\w+\/\d+)\/\w+/)
	{
        # PacBio as of feb 2014
		$this->{base}    = $1;
	} elsif ($hdr =~ /^@(\S+) (\d?):([YN]):(\d+):([aAtTcCgGnN]*)/)
	{
        # current Illumina format (Casava 1.8.2)
		$this->{base}    = $1;
		$this->{pair}    = $2 if $2;
        $this->{filtered} = $3 eq 'Y' ? "casava" : undef;
        $this->{control} = $4;
		$this->{barcode} = uc($5);
	} elsif ($hdr =~ /^@(\S+)#([aAtTcCgGnN]+)\/([12])/)
	{
        # old Illumina format, JGI variant
		$this->{base}    = $1;
		$this->{barcode} = uc($2);
		$this->{pair}    = $3;
	} elsif ($hdr =~ /^@(\S+)\/([12])#([aAtTcCgGnN]+)/)
	{
        # old Illumina format
		$this->{base}    = $1;
		$this->{pair}    = $2;
		$this->{barcode} = uc($3);
	} elsif ($hdr =~ /^@(\S+)\/([12])/)
	{
        # old Illumina format
		$this->{base} = $1;
		$this->{pair} = $2;
	} elsif ($hdr =~ /^@(\S+)#([aAtTcCgGnN]+)/)
	{
        # old Illumina format
		$this->{base}    = $1;
		$this->{barcode} = uc($2);
	} elsif ($hdr =~ /^@(\S+)/)
	{
        # old Illumina format
		$this->{base} = $1;
	} elsif ($hdr =~ /^@(\S+) length=\d+ xy=\d+_\d+ region=\d run=\S/)
	{
		# Roche/454
		$this->{base} = $1;
	} elsif ($hdr =~ /^@(\S+)/)
	{
		# other/unknown (pairing, barcodes not parsable)
		$this->{base} = $1;
	} else
	{
		confess("Unable to parse sequence header: $hdr\n");
	}

    # DEFINE ID
    my $base = $this->{base};
    my $pair = defined($this->{pair}) ? $this->{pair} : '';
    my $filtered = defined($this->{filtered}) ? 'Y' : 'N';
    my $control = defined($this->{control}) ? $this->{control} : 0;
    my $barcode = defined($this->{barcode}) ? $this->{barcode} : '';
    $this->{id} = "$base $pair:$filtered:$control:$barcode";
}

=item header <$reformat>

Returns the object's header line.  Optional flag will output original ID rather than standard format (default).

=cut

sub header
{
    my ($this, $orig)=@_;
    return $orig ? $this->{header} : '@'.$this->{id};
}

=item id

Returns the object's ID, which is the sequence's unique identifier without comments which may be present in the header.
It cannot be changed directly; set header, base, barcode, or pair instead.

=cut

sub id { return shift->{id} }

=item base

Returns the base of the ID; both paired reads will have the same base.
Optionally set base of the ID.

=cut

sub base
{
    my ($this, $base) = @_;
    $this->{base} = $base if $base;
    return $this->{base};
}

=item barcode

If the read contains an Illumina molecular barcode ID, it will be returned; otherwise returns undef.
Supplying an optional argument will set the barcode; passing an empty string will clear the barcode.

=cut

sub barcode
{
	my ($this, $barcode) = @_;
	if (defined($barcode))
	{
		$this->{barcode} = $barcode ? $barcode : undef;
	}
	return $this->{barcode};
}

=item pair

Returns the read's ord in the pair; undef otherwise.

=cut

sub pair
{
    my ($this, $pair) = @_;
    if (defined($pair) and ($pair eq '' or $pair == 1 or $pair == 2) )
    {
        $this->{pair} = $pair;
    }
    return $this->{pair};
}

=item unpair

Method to clear pairing of a read (e.g. singleton).

=cut

sub unpair { shift->{pair} = undef }

=item seq

Returns the read's complete sequence, without newlines.

=cut

sub seq
{
	my $this = shift;
	return $this->{filtered} ? '' : $this->{seq};
}

=item len

Returns the length of the sequence (returns 0 if filtered).

=cut

sub len
{
	my $this = shift;
	return $this->{filtered} ? 0 : length($this->{seq});
}

=item revcomp

Reverse-complements a sequence and quality scores.

=cut

sub revcomp
{
	my $this = shift;
	return unless $this->{seq};
	$this->{seq} =~ tr/ATCGatcg/TAGCtagc/;
	my @seq = reverse split(//, $this->{seq});
	$this->{seq} = join('', @seq);
	my @qual = reverse split(//, $this->{qual});
	$this->{qual} = join('', @qual);
}

=item qual ($new_qual)

Returns the read's quality string, without newlines.

=cut

sub qual
{
	my $this = shift;
	return $this->{filtered} ? '' : $this->{qual};
}

=item qual_arrayref

Returns an arrayref of sanger-scaled quality scores.

=cut

sub qual_arrayref
{
	my $this = shift;
	return [] unless $this->{qual};
	return [] if $this->{filtered};
	my @qual = map { ord($_) - 33 } split(//, $this->{qual});
	return \@qual;
}

=item qual_stats

Returns statistics on quality scores.

=cut

sub qual_stats
{
    my $this = shift;
    return undef if $this->{filtered};
    return undef unless $this->{qual};
	my @qual = map { ord($_) - 33 } split(//, $this->{qual});
    my $n = scalar(@qual);
    # odd median
    my $oddMedian = oddMedian(\@qual);
    # min,max,sum
    my $min = 32767;
    my $max = -1;
    my $sum = 0;
    my %count = ();
    foreach my $q (@qual)
    {
        $min = $q if $q < $min;
        $max = $q if $q > $max;
        $sum += $q;
        $count{$q}++;
    }
    # mean
    my $mean = $sum/$n;
    # std dev
    my $d = 0;
    foreach my $q (@qual)
    {
        $d += ( $q - $mean ) **2;
    }
    my $stdev = $d/$n;
    # mode
    my @result = ();
    foreach (sort { $count{$b} <=> $count{$a} } keys %count)
    {
        last if @result && $count{$_} != $count{$result[0]};
        push(@result, $_);
    }
    my $mode = oddMedian(\@result);
    return ($n,$min,$max,$oddMedian,$mode,round($mean),round($stdev),\%count);

}

=item del

Mark this record as filtered.  Always returns 0.

=cut

sub del
{
	my ($this, $reason) = @_;
	confess("Missing filter reason\n") unless $reason;
    $this->{seq} = $this->{qual} = '';
	$this->{filtered} = $reason;
    return 0;
}

=item filtered

Returns the reason why the read was filtered; undef otherwise.

=cut

sub filtered { return shift->{filtered} }

=item pass

Returns 1 if read passes QC, 0 otherwise.

=cut

sub pass { return shift->{filtered} ? 0 : 1 }

=item output

Returns a multiline string of the sequence in Fastq format.  Filtered unpaired reads return an empty string, while filtered paired reads return an empty record, so pairing won't be broken.  Optional flag to use original read ID rather than standard format (default).

=cut

sub output
{
	my ($this, $orig) = @_;
    return '' if $this->filtered;
    return $this->header($orig) . "\n" . $this->{seq} . "\n+\n" . $this->{qual} . "\n";
}

=item output_fasta

Returns a multiline string of the sequence in Fasta format if read passed QC filters.
in the header.

=cut

sub output_fasta
{
	my $this = shift;
    return '' if $this->{filtered} or !$this->len;
    my $fasta = '>'.$this->id."\n";
    my $seq = $this->{seq};
    while (length($seq) > CHARACTERS_PER_LINE)
    {
        $fasta .= substr($seq, 0, CHARACTERS_PER_LINE) . "\n";
        $seq = substr($seq, CHARACTERS_PER_LINE);
    }
    $fasta .= $seq . "\n";
    return $fasta;
}

=item output_qual

Returns a multiline string of the sequence's quality scores in phred format only if read passed QC filters.

=cut

sub output_qual
{
	my $this = shift;
    return '' if $this->filtered or !$this->len;
    my @qual   = map { ord($_) - 33 } split(//, $this->{qual});
    my $output = ">" . $this->id . "\n";
    my $i      = CHARACTERS_PER_LINE - 1;
    while (scalar(@qual) > CHARACTERS_PER_LINE)
    {
        $output .= join(' ', @qual[ 0 .. $i ]) . "\n";
        @qual = @qual[ CHARACTERS_PER_LINE .. $#qual ];
    }
    $output .= join(' ', @qual) . "\n";
    return $output;
}

###############################################################################
## QC METHODS

=back

=head2 QC Methods

Each function returns 1 if read passes filter, 0 otherwise.

=over 5

=item qc

To perform QC steps defined in %qc_params.

=cut

sub qc
{
	my ($this) = @_;
    return 1 if $this->{qc}->{no_qc};
    return 0 unless $this->min_len($this->{qc}->{min_len});
    return 0 unless $this->trim3($this->{qc}->{trim3});
    return 0 unless $this->trim5_bp($this->{qc}->{trim5_bp});
    return 0 unless $this->trim3_bp($this->{qc}->{trim3_bp});
    return 0 unless $this->trim_to($this->{qc}->{trim_to});
    return 0 unless $this->trim3_exp_err($this->{qc}->{trim3_exp_err});
    return 0 unless $this->max_n($this->{qc}->{max_n});
    return 0 unless $this->low_complexity($this->{qc}->{low_complexity});
    return 0 unless $this->low_q($this->{qc}->{low_q});
    return 0 unless $this->exp_err($this->{qc}->{exp_err});
    return 0 unless $this->exp_err_per_kb($this->{qc}->{exp_err_per_kb});
    return 1;
}

=item trim3

Trim uncalled bases from 3' end of read (i.e. bp=N or qual=2).

=cut

sub trim3
{
	my ($this) = @_;
    my $seq = $this->seq;
    my $qual = $this->qual;
    my $len = length($seq);
    if ( $seq =~ s/N+$//g )
    {
        $qual = substr($qual, 0, length($seq));
    }
    if ( $qual =~ s/#+$//g )
    {
        $seq = substr($seq, 0, length($qual));
    }
    $this->{seq} = $seq;
    $this->{qual} = $qual;
    return $this->min_len($this->{qc}->{min_len});
}

=item trim5_bp

Remove the specified number of bases from the 5' end.

=cut

sub trim5_bp
{
    my ($this, $trim5_bp) = @_;
    return 1 unless defined($trim5_bp);
    my $len = $this->len - $trim5_bp;
    if ( $len > 0 )
    {
        $this->{seq} = substr($this->seq, $trim5_bp);
        $this->{qual} = substr($this->qual, $trim5_bp);
    } else
    {
        $this->{seq} = $this->{qual} = '';
    }
    return $this->min_len($this->{qc}->{min_len});
}

=item trim3_bp

Remove the specified number of bases from the 3' end.

=cut

sub trim3_bp
{
    my ($this, $trim3_bp) = @_;
    return 1 unless $trim3_bp;
    my $len = $this->len - $trim3_bp;
    if ( $len > 0 )
    {
        $this->{seq} = substr($this->seq, 0, $len);
        $this->{qual} = substr($this->qual, 0, $len);
    } else
    {
        $this->{seq} = $this->{qual} = '';
    }
    return $this->min_len($this->{qc}->{min_len});
}

=item min_len

If the sequence is shorter than the minimum length, the record is filtered.

=cut

sub min_len
{
    my ($this, $min_len) = @_;
    return 1 unless $min_len;
	return $this->del('too_short') if $this->len < $min_len;
    return 1;
}

=item max_n

If the sequence contains more than the allowed number of Ns, the record is filtered.

=cut

sub max_n
{
    my ($this, $max_n) = @_;
    return 1 unless defined($max_n);
	my $n = $this->{seq} =~ s/N/N/gi;
	return $this->del('too_many_N') if $n > $max_n;
    return $this->min_len($this->{qc}->{min_len});
}

=item low_complexity $pct_len

If the sequence is >= $pct_len mono- or di-nucleotide repeats, the record is filtered.

=cut

sub low_complexity
{
    my ($this, $pct_len) = @_;
    return 1 unless defined($pct_len);
	my $seq = $this->{seq};
	my $len = length($seq);
	return 0 unless $len;
	foreach my $nn (qw/AA TT CC GG CA GT CT GA AT CG/)
	{
		my $n = $seq =~ s/$nn/$nn/gi;
		return $this->del('low_complexity') if ($n*2/$len) >= $pct_len;
	}
    return 1;
}

=item low_q

Filter record if there are too many low-quality bases. Default minq=10.

=cut

sub low_q
{
    my ($this, $max ) = @_;
    return 1 unless $max;
	my $minq = $this->{qc}->{low_q};
    $minq = 10 unless $minq;
	my $n = 0;
	foreach my $q (@{$this->qual_arrayref}) { ++$n if $q < $minq }
	return $this->del('too_many_low_qual') if $n > $max;
    return 1;
}

=item trim3_exp_err

Trim 3' end of read once total expected error of 5' bases reaches max (default = 0.5).  Fail if read has less than minimum specified number of bp (default = 20).

=cut

sub trim3_exp_err
{
    my ($this, $max_exp_err) = @_;
    return 1 unless defined($max_exp_err);
    my $exp_err = 0;
    my $len = 0;
    foreach (split(//, $this->{qual}))
    {
        $exp_err += 10 ** (((ord($_) - 33) * -1)/10);
        last if $exp_err > $max_exp_err;
        ++$len;
    }
    $this->{seq} = substr($this->seq,0,$len);
    $this->{qual} = substr($this->qual,0,$len);
    return $this->min_len($this->{qc}->{min_len});
}

=item exp_err

Fail read if expected number of errors exceeds threshold.

=cut

sub exp_err
{
    my ($this, $max_exp_err) = @_;
    return 1 unless $max_exp_err;
    my $exp_err = 0;
    my @qs = split(//, $this->{qual});
    foreach my $s (@qs)
    {
        my $q = ord($s) - 33;
        my $ee = 10 ** (($q * -1)/10);
        $exp_err += $ee;
        return $this->del('expected_error') if $exp_err > $max_exp_err;
    }
    return 1;
}

=item exp_err_per_kb

Fail read if expected number of errors exceeds threshold.

=cut

sub exp_err_per_kb
{
    my ($this, $max_exp_err_per_kb) = @_;
    return 1 unless $max_exp_err_per_kb;
    my $max_exp_err = int($this->len * $max_exp_err_per_kb/1000 + 0.5);
    my $exp_err = 0;
    my @qs = split(//, $this->{qual});
    foreach my $s (@qs)
    {
        my $q = ord($s) - 33;
        my $ee = 10 ** (($q * -1)/10);
        $exp_err += $ee;
        return $this->del('expected_error_per_kb') if $exp_err > $max_exp_err;
    }
    return 1;
}

=item trim_to

Trim read to exactly this many bp by discarding 3' bases.

=cut

sub trim_to
{
    my ($this, $len) = @_;
    return 1 unless $len;
    return 0 unless $this->len >= $len;
    $this->{seq} = substr($this->seq,0,$len);
    $this->{qual} = substr($this->qual,0,$len);
    return 1;
}

=back

=head1 BUGS AND LIMITATIONS

Reads must be named in accordance with Illumina naming conventions.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut

1;
__END__
