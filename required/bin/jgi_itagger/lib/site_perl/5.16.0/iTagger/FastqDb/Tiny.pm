=pod

=head1 NAME

iTagger::FastqDb::Tiny

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

This module provides an iterator for Fastq files, plus basic QC functions.
This Tiny version of FastqDb does not buffer input, does not auto-detect
qual-encoding, and does not support barcodes.

=head1 METHODS

=over 5

=cut

package iTagger::FastqDb::Tiny;

use warnings;
use strict;
use Carp;
use IO::File;
use File::Spec;
use Env qw/TMPDIR/;
use iTagger::FastqDb::Fastq;

our $VERSION = 1.0;

=item new $file

Constructor accepts a path to an existing Fastq file (otherwise reads from standard input).

=cut

sub new
{
	my ( $class, $file, $params ) = @_;
	if ( defined($file) )
	{
		if ( $file eq '-' ) { $file = undef }
        elsif ( $file !~ /\|\s*$/ )
		{
			$file = readlink($file) if -l $file;
			$file = File::Spec->rel2abs($file);
			confess("Input file, $file, not found\n") unless -e $file;
		}
	}
	my $this = {
		file => $file,
		line => 0,
		paired => undef,
		counter => { pass  => 0, error => 0 },
		buffer => [] };
	bless $this, $class;
	$this->_init_qc_params($params);
	$this->_open_file;
	return $this;
}

=item DESTROY

Close file before going away

=cut

sub DESTROY
{
	my $this = shift;
	$this->close_file;
}

=item _open_file

Open file, possibly compressed, or stdin.

=cut

sub _open_file
{
	my $this = shift;
	$this->close_file if exists( $this->{fhr} );
	my $file = $this->{file};
	my $fh   = new IO::File;
	if ( !$file )
	{
		$fh = *STDIN;
	}
	elsif ( $file =~ /\.gz$/ )
	{
		open( $fh, "gunzip -c $file |" ) or confess("Unable to open gzipped infile, $file: $!\n");
	}
	elsif ( $file =~ /\.bz2$/ )
	{
		open( $fh, "bzcat $file |" ) or confess("Unable to open bzipped infile, $file: $!\n");
	}
	elsif ( -e $file )
	{
		open( $fh, '<', $file ) or confess("Unable to open infile, $file: $!\n");
	}
	else
	{
		open( $fh, $file ) or confess("Unable to open infile, $file: $!\n");
	}
	$this->{fhr}  = \$fh;
}

=item close_file

Close file

=cut

sub close_file
{
	my $this = shift;
	return unless exists( $this->{fhr} );
	if ( defined( $this->{fhr} ) )
	{
		my $fh = ${ $this->{fhr} };
		close($fh) if $this->{file};
	}
	delete $this->{fhr};
}

=item _init_qc_params

Initialize QC parameters hash with default values.

=cut

sub _init_qc_params
{
	my ( $this, $params ) = @_;
	$params = {} unless defined($params);
    # DEFAULT VALUES
	my %qc_params = (
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
        trim_to        => undef     # trim read to exactly this many bp
	);

	foreach my $key ( keys %$params )
	{
        if ( $key eq 'paired' ) { $this->{paired} = $params->{$key} ? 1:0 }
        else { $qc_params{$key} = $params->{$key} }
    }
	$this->{qc_params} = \%qc_params;
}

=item next_seq

An iterator for the Fastq database.  Returns the next valid read object.
If the reads are paired and one read is filtered, the entire pair is filtered.
The &next_pair method is also useful for returning both pairs at once.

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

Get next single read

=cut

sub _next_single
{
	my ($this) = @_;
	my $rec;
	while (1)
	{
		$rec = $this->_next_seq;
		last unless defined($rec);
		if ( $rec->filtered )
		{
			$this->_tally( $rec->filtered );
		} else
        {
            $this->{counter}->{pass} += 1;
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
        }
    } until ( $rec1 and $rec2 );
	return [ $rec1, $rec2  ];
}

=item _tally

Update counter for a filtered read

=cut

sub _tally
{
	my ( $this, $reason ) = @_;
	if ( exists( $this->{counter}->{$reason} ) )
	{
		$this->{counter}->{$reason} += 1;
	}
	else
	{
		$this->{counter}->{$reason} = 1;
	}
}

=item _next_seq

Get next seq from buffer or input

=cut

sub _next_seq
{
    my ($this) = @_;
	my $rec = scalar( @{ $this->{buffer} } ) ? shift @{ $this->{buffer} } : undef;
    return $rec if $rec;
	return undef unless exists( $this->{fhr} ) and defined( $this->{fhr} );
	my $fh = ${ $this->{fhr} };
    my $hdr = <$fh>;
	unless ( defined($hdr) )
	{
		$this->close_file;
		return undef;
	}
    chomp $hdr;
    my $seq = <$fh>;
    confess("Truncated Fastq file; incomplete record for: $hdr\n") unless defined($seq);
    chomp $seq;
    my $sep = <$fh>;
    confess("Truncated Fastq file; incomplete record for: $hdr\n") unless defined($sep);
    my $qual = <$fh>;
    confess("Truncated Fastq file; incomplete record for: $hdr\n") unless defined($qual);
    chomp $qual;
	return new iTagger::FastqDb::Fastq( $hdr, $seq, $qual, $this->{qc_params});
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

Return number of reads which pass, fail filters.

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

1;

=back

=head1 COPYRIGHT

Copyright (c) 2010 U.S. Department of Energy Joint Genome Institute

All right reserved.

=head1 AUTHOR

Edward Kirton <ESKirton@LBL.gov>

=cut
