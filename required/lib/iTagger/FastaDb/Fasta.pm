=pod

=head1 NAME

Fasta - Simple object for Fasta sequence

=head1 SYNOPSIS

    my $seq=new Fasta( $hdr, $seq );
    $seq->qc;
    print $seq->output;

=head1 DESCRIPTION

Object for a single read sequence, with methods for basic manipulation and quality control.

=head1 METHODS

=over 5

=cut

package iTagger::FastaDb::Fasta;

use warnings;
use strict;
use Carp;
use constant {
    CHARACTERS_PER_LINE => 80, # for formatting Fasta output
    CLOSE_ENOUGH_TO_END => 6,   # hits of adapters this close to end will result in trimming to be done to end
};

our $VERSION = 1.0;

=item new $hdr $seq

Initialize new sequence object.

=cut

sub new {
    my ($class,$hdr,$seq)=@_;
    chomp $hdr;
    chomp $seq;
    confess("Incomplete Fasta record: hdr=$hdr, seq=$seq\n") unless $hdr and $seq;
    # INIT
    my $this={
        seq=>$seq
    };
    bless $this,$class;
    $this->_parse_header($hdr);
    return $this;
}

=item _parse_header

Extract info from header, which could be formatted in several ways.

=cut

sub _parse_header {
    my ($this, $hdr) = @_;
	$hdr = '>' . $hdr unless $hdr =~ /^>/;
    $this->{header} = $hdr;
    $this->{id}=undef; # complete ID (e.g. "A/1#GATTACA"); always defined
    $this->{base}=undef; # base ID only (e.g. "A"); always defined
    $this->{pair}=undef; # pair ID only (e.g. "1"); only defined if paired
    $this->{barcode}=undef; # barcode sequence (upper-case); only defined it barcoded

	if ($hdr =~ /^>(\w+_\d+_\d+_\w+_\w+_\w+\/\d+)\/\w+/)
	{
        # PacBio as of feb 2014
		$this->{base}    = $1;
    } elsif ($hdr =~ /^>(\S+)#([aAtTcCgGnN]+)\/([12])/) { # Illumina barcoded, paired
        $this->{base}=$1;
        $this->{barcode}=uc($2);
        $this->{pair}=$3;
    } elsif ($hdr =~ /^>(\S+)\/([12])#([aAtTcCgGnN]+)/) { # Illumina barcoded, paired
        $this->{base}=$1;
        $this->{pair}=$2;
        $this->{barcode}=uc($3);
    } elsif ($hdr =~ /^>(\S+)\/([12])/) { # Illumina paired
        $this->{base}=$1;
        $this->{pair}=$2;
    } elsif ($hdr =~ /^>(\S+)#([aAtTcCgGnN]+)/) { # Illumina barcoded, unpaired
        $this->{base}=$1;
        $this->{barcode}=uc($2);
    } elsif ( $hdr =~ /^>(\S+)/ )
    {
        $this->{base}=$1; # general; parsing barcode/pairing not possible
    }

    $this->{id} = $this->{base};
    $this->{id} .= "/".$this->{pair} if $this->{pair};
    $this->{id} .= "#".$this->{barcode} if $this->{barcode};
}

=item header <$reformat>

Returns the object's header line.  Optional flag will cause original header to be used rather than standard format (default).

=cut

sub header
{
    my ($this, $orig)=@_;
    return $orig ? $this->{header} : '>'.$this->{id};
}

=item id

Returns the object's ID, which is the sequence's unique identifier without comments which may be present in the header.
It cannot be changed directly, but will be updated whenever the header, base, or pair is changed.

=cut

sub id { return shift->{id} }

=item base

If the read is paired, returns it's base ID (same for both members of pair); returns undefined if not paired.
Providing the optional argument will change the base (and the id and header).

=cut

sub base {
    my ($this,$base)=@_;
    if ( defined($base) ) {
        $this->{id}=$this->{base}=$base;
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$this->{barcode} if $this->{barcode}
    }
    return $this->{base};
}

=item barcode

If the read contains an Illumina molecular barcode ID, it will be returned; otherwise returns undef.
Supplying an optional argument will set the barcode.

=cut

sub barcode {
    my ($this,$barcode)=@_;
    if ($barcode) {
        $this->{barcode}=$barcode;
        $this->{id}=$this->{base};
        $this->{id} .= "/".$this->{pair} if $this->{pair};
        $this->{id} .= "#".$barcode;
    }
    return $this->{barcode};
}

=item pair

Returns the read's ord in the pair; undef otherwise.  It may be changed by supplying the optional extra argument.
To clear the pairing information, pass it an empty string, not NULL.

=cut

sub pair {
    my ($this,$pair)=@_;
    if (defined($pair)) {
        if ($pair) {
            $this->{pair}=$pair;
            $this->{id}=$this->{base};
            $this->{id} .= "/".$pair;
            $this->{id} .= "#".$this->{barcode} if $this->{barcode}
        } else {
            # delete pairing information (e.g. singleton)
            $this->{pair}=undef;
            $this->{id}=$this->{base};
            $this->{id} .= "#".$this->{barcode} if $this->{barcode}
        }
    }
    return $this->{pair};
}

=item seq ($new_seq)

Returns the read's complete sequence, without newlines.  Optional argument changes it.

=cut

sub seq {
    my ($this,$new_seq)=@_;
    if ($new_seq) {
        chomp $new_seq;
        $this->{seq}=$new_seq;
    }
    return $this->{seq};
}

=item revcomp

Reverse-complements a sequence.

=cut

sub revcomp {
    my $this=shift;
    return unless $this->{seq};
    $this->{seq} =~ tr/ATCGatcg/TAGCtagc/;
    my @seq=reverse split(//, $this->{seq});
    $this->{seq}=join('', @seq);
}

=item output

Returns a multiline string of the sequence in Fasta format.  Returns no output if sequence is empty.  Optional flag will cause original header to be output rather than standard format (default).

=cut

sub output {
    my ($this, $orig) = @_;
    return '' unless $this->{seq}; # will be empty if failed QC
    return $this->header($orig)."\n"._format($this->{seq});
}

# PRIVATE METHOD TO ADD NEWLINE CHARACTERS 
sub _format {
    my $old=shift;
    return '' unless $old;
    return $old unless CHARACTERS_PER_LINE;
    my $new='';
    while (length($old)> CHARACTERS_PER_LINE) {
        $new .= substr($old,0, CHARACTERS_PER_LINE)."\n";
        $old = substr($old, CHARACTERS_PER_LINE);
    }
    $new .= $old."\n";
    return $new;
}

=item len

Returns length of the seq

=cut

sub len {
    my $this = shift;
    return length($this->{seq});
}

=item qc $winsize $meanq $minlen $maxn

To perform minimum length filtering, and filtering reads with too many Ns.

=cut

sub qc {
    my ($this, $minlen, $maxn)=@_;
    $this->trim_terminal_Ns;
    $this->length_filter($minlen);
    $this->N_filter($maxn);
}

=item trim_terminal_Ns

Discard uninformative Ns from the ends of the sequence.

=cut

sub trim_terminal_Ns {
    my $this=shift;
    if ($this->seq =~ /^(N+)(.*)$/) {
        $this->seq($2);
    }
    if ($this->seq =~ /^(.*)N+$/) {
        $this->seq($1);
    }
}

=item length_filter

If the sequence is shorter than the minimum length, the sequence string is emptied so they will not be
returned by the output method.  Returns true if sequence was filtered.

=cut

sub length_filter {
    my ($this,$minlen)=@_;
    my $len=length($this->{seq});
    if ( !defined($minlen) or $len >= $minlen) { 
        return 0;
    } else {
        $this->{seq}='';
        return 1;
    }
}

=item N_filter

If the sequence contains more than the allowed number of Ns, the sequence string is emptied.
Returns true if sequence was filtered.

=cut

sub N_filter {
    my ($this,$maxn)=@_;
    return 1 unless defined($maxn);
    # COUNT NUMBER OF 'N' BASES IN TRIMMED SEQ
    my $tmpseq=$this->{seq};
    my $n= $tmpseq=~ s/N//gi;
    if ($n <= $maxn) {
        return 0;
    } else {
        # FAIL
        $this->{seq}='';
        return 1;
    }
}

=item low_complexity_filter $pct_len

If the sequence is >= $pct_len mono- or di-nucleotide repeats, clears the sequence string and returns true.

=cut

sub low_complexity_filter {
    my ($this,$pct_len)=@_;
    my $seq=$this->{seq};    
    $seq =~ s/\n//g;
    my $len = length($seq);
    my $filter=0;
    foreach my $nn (qw/AA TT CC GG CA GT CT GA AT CG/) {
        my $n = $seq =~ s/$nn/$nn/g;
        if ($n >= $pct_len/200*$len) {
            $filter = 1;
            last;
        }
    }
    if ($filter) {
        $this->{seq}='';
    }
    return $filter;
}

#=item trim
#
#Given start and end coordinates of adapter/primer sequence, trims sequence string and returns final length.
#
#=cut
#
#sub trim {
#    my ($this, $start, $end)=@_;
#    return unless defined($start) and defined($end);
#    ($start,$end)=($end,$start) if $end<$start; 
#    my $len=length($this->{seq});
#    #print STDERR "- trim read, ", $this->id, " ($len bp) from $start-$end\n"; # DEBUG
#
#    if ($start <= CLOSE_ENOUGH_TO_END and $len - $end <= CLOSE_ENOUGH_TO_END ) {
#        # fail read
#        $this->{seq}='';
#    } elsif ($start <= CLOSE_ENOUGH_TO_END ) {
#        # trim left
#        $this->{seq}=substr($this->{seq},$end);
#    } elsif ($len - $end <= CLOSE_ENOUGH_TO_END ) {
#        # trim right
#        $this->{seq}=substr($this->{seq},0,$start);
#    }
#    return length($this->{seq});
#    #print STDERR "\t+ length after trimming is ", length($this->{seq})," bp\n"; # DEBUG
#}

=back

=head1 BUGS AND LIMITATIONS

No support for paired reads.

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved. This program is free software; you can redistribute it
and/or modify it under the same terms as Perl itself.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut

1;
__END__
