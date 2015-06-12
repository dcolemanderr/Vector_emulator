=pod

=head1 NAME

FastaDb - Simple object for handling Fasta files

=head1 SYNOPSIS

    # simple script to QC reads (outputs to stdout)
    use FastaDb;
    my $file=shift or confess("Infile required\n");
    my $fastaFile = new FastaFile($file);
    while (my $fasta=$fastaFile->next_seq) {
        $fasta->qc; # see Fasta package for details
        print $fasta->output; # no output if read failed QC
    }

=head1 DESCRIPTION

This module provides an iterator for Fasta files, plus basic QC functions

=head1 METHODS

=over 5

=cut

package iTagger::FastaDb;

use warnings;
use strict;
use Carp;
use IO::File;
use File::Basename;
use Cwd;
use iTagger::FastaDb::Fasta;

our $VERSION = 1.0;

=item new

Construct new object from file path.

=cut

sub new {
    my ($class,$file)=@_;

    my $fh = new IO::File;
    if ( !$file )
    {
        $fh = *STDIN;
    } else {
        # USE FULL PATH
        my $cwd=getcwd.'/';
        my ($base,$dir,$suffix)=fileparse($file, qr/\.[^.]*/);
        $dir=$cwd if $dir eq './';
        $file="$dir$base$suffix";
        if ( $file =~ /\.gz$/ )
        {
            open($fh, "gunzip -c $file|") or confess($!);
        } elsif ( $file =~ /\.bz2$/ )
        {
            open($fh, "bzcat $file|") or confess($!);
        } else
        {
            open($fh, '<', $file) or confess("Unable to open fasta file, $file: $!\n");
        }
    }
 
    my $this = { file => $file, next_header => undef, fhr => \$fh };
    bless $this,$class;
    return $this;
}

=item DESTROY

Closes file before going away

=cut

sub DESTROY {
    my $this=shift;
    return unless exists($this->{fhr}) and defined($this->{fhr});
    my $fh=${$this->{fhr}};
    close $fh;
}

=item file ($new_file)

Returns the Fasta db complete path.  You may also use this method to open a new file for reading.

=cut

sub file {
    my ($this,$new_file)=@_;
    if ($new_file) {
        $this->close_file;
        my $fh=new IO::File('<'.$new_file) or confess("Unable to open Fasta file, $new_file: $!\n");
        $this->{fhr}=\$fh;
    }
    return $this->{file};
}

=item next_seq

An iterator for the fasta database; returns a single sequence as multiline text.

=cut

sub next_seq {
    my $this=shift;
    return undef unless exists($this->{fhr}) and defined($this->{fhr});

    # INIT VARS
    my $fh = ${$this->{fhr}};
    my $hdr=undef;
    my $seq='';
    if (defined($this->{next_header})) {
        $hdr=$this->{next_header};
        $this->{next_header}=undef;
    } else {
        # GET FIRST HEADER
        while (<$fh>) {
            chomp;
            if (/^>\S+/) {
                $hdr=$_;
                last;
            } elsif (/^#/ or ! $_) {
                next;
            } else {
                confess("Invalid Fasta input file; choked on line: $_\n");
            }
        }
    }
    unless (defined($hdr)) {
        $this->close_file;
        return undef;
    }

    # GET SEQ
    while (<$fh>) {
        chomp;
        if (/^>\S+/) {
            $this->{next_header}=$_;
            last;
        } elsif($_) {
            $seq .= $_;
        }
    }

    # CHECK IF END OF FILE OR PART
    if ( !defined($this->{next_header}) ) {
        $this->close_file;
    }

    # RETURN ONE FASTA SEQ
    return new iTagger::FastaDb::Fasta($hdr,$seq);
}

=item close_file

This will close the file and delete the file handle.  Normally this is called by the iterator upon reaching the end
of the file, so this method only needs to be called if you wish to close the file prematurely.

=cut

sub close_file {
    my $this=shift;
    return unless exists($this->{fhr}) and defined($this->{fhr});
    my $fh=${$this->{fhr}};
    close $fh;
    delete $this->{fhr};

}

1;

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
