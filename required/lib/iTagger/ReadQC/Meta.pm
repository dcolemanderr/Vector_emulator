=pod

=head1 NAME

iTagger::ReadQC::Meta

=head1 DESCRIPTION

Meta-analysis of many QC'ed libraries.  Produces summary and flags outliers.

=head1 FUNCTIONS

=over 5

=cut

package iTagger::ReadQC::Meta;

use strict;
use warnings;
use Carp;
use Env qw(TMPDIR ITAGGER_LOG_CONFIG);
use File::Copy;
use Config::Tiny;
use File::Spec::Functions qw(rel2abs);
use File::Basename qw(fileparse);
use constant { TARGET_STD_DEV => 2 };
require Exporter;
use iTagger::FastqDb;
use iTagger::FastqDb::Tiny;
use iTagger::Stats;
use Data::Dumper;
our $VERSION = 1.1;
our @ISA = qw(Exporter);
our @EXPORT = qw(new summary);

=item new

Constructor.  Config is either a hashref or the path to an .ini format file.

=cut

sub new
{
    my ($class, $config, $dir) = @_;
    confess("Missing args") unless $config and $dir;

    # INIT
    my $this = {};
    bless $this, $class;

    # OUTPUT DIR
    $this->{dir} = $dir;

    # CONFIG
    if ( ref($config) eq '' )
    {
        $config = Config::Tiny->read($config);

    } elsif ( ref($config) ne 'HASH' )

    {
        confess("Invalid config parameter; hashref or .ini file required");
    }

#	print Dumper $config;

    foreach my $section (qw/AMPLICON AMPQC/)
    {
        confess("Config file missing required section: $section") unless exists($config->{$section});
    }

    $this->{config} = $config;
    
    # LOGFILE
    my $logFile = "$dir/readQC.log";

    # PAIRED OR NOT?
    $this->{paired} = ( exists($config->{FLASH}) or exists($config->{PANDASEQ}) ) ? 1 : 0;

    # VERIFY OUTPUT DIR EXISTS
    confess("The specified folder does not exist: $dir") unless -e $dir;;

    return $this;
}

=item DESTROY

Close output file.

=cut

sub DESTROY
{
    my ($this) = @_;
    if ( exists($this->{fh}) and defined($this->{fh}) )
    {
        my $out = $this->{fh};
        close($out);
    }
}

=item summary

Create summary of libraries' read QC logs and identify outliers.

=cut

sub summary
{
    my ($this) = @_;
    my $config = $this->{config};
    my $dir = $this->{dir};

    # FIND LIBRARIES
    opendir(DIR, $this->{dir}) or confess($!);
    my @libs = grep { $_ !~ /^\./ } readdir(DIR);
    closedir(DIR);

    # READ LOGS
    my %summary = ();
    foreach my $name (@libs)
    {
        next if $name eq 'readQC.log';
        next unless -d "$dir/$name";
        my $file = "$dir/$name/readQC.log";
        unless ( -e $file )
        {
            confess("Lib $name missing readQC.log; there was probably an error processing this library -- rerun and try again.");
        }
        open(my $in, '<', $file) or confess("Unable to open file, $file: $!");
        my $hdr = <$in>;
        while (<$in>)
        {
            chomp;
            my ($key, $value) = split(/\t/);
            $value =~ s/,//g;
            $value = 0 unless $value;
            if ( exists($summary{$key}) )
            {
                $summary{$key} += $value;
            } else
            {
                $summary{$key} = $value;
            }
        }
        close($in);
    }
    $this->{summary} = \%summary;

    # OPEN FILE, WRITE HEADER
    open(my $out, '>', "$dir/readQC.log") or confess($!);
    print $out "#STEP\tCOUNT\tPCT OF PREVIOUS STEP\n";
    $this->{fh} = $out;

    # READS INPUT.  $num will always be current number of reads (i.e. used as denominator in following step).
    my $num = $this->log('Input', undef);
    return unless $num;

    # FILTER CONTAMINANTS.  THERE MAY BE ZERO OR MORE CONTAM DB FILES.
    if ( exists($config->{DUK}->{CONTAM_DB}) )
    {
        my @contamDbs = split(/,/, $config->{DUK}->{CONTAM_DB});
        my @dbNames = ();
        for ( my $i=0; $i<=$#contamDbs; $i++ )
        {
            my $contamDb = $contamDbs[$i];
            last unless $contamDb;
            my ($dbName, $dbDir, $dbSuffix) = fileparse($contamDb, qr/\.[^.]*/);
            my $aKey = "Filtered $dbName";
            return unless $num = $this->log($aKey, $num);
        }
    }

    # TRIM PRIMERS USING CUTADAPT
    if ( exists($config->{CUTADAPT}) )
    {
        return unless $num = $this->log('Primer-trimmed', $num);
    }

    # UNPAIRED READS
    unless ( $this->{paired} )
    {
        # OPTIONAL HARD-TRIM AND LENGTH FILTER
        return unless $num = $this->log('Length-filtered', $num);

        # EXPECTED ERROR FILTER
        return unless $num = $this->log('Quality-filtered', $num);

        # DEREPLICATION
        $num = $this->log('High quality, unique', $num);
        return;
    }

    ## PAIRED READS

    # MERGE PAIRED READS
    return unless $num = $this->log('Extended', $num);

    # HARD-TRIM
    return unless $num = $this->log('Length-filtered', $num);

    # EXPECTED ERROR FILTER
    return unless $num = $this->log("High quality", $num);

    # DEREPLICATE
    $num = $this->log("Unique", $num);
    return;
}

=item log

Method calculates percentage and writes single line of summary table to output file.

=cut

sub log
{
    my ($this, $key, $denom) = @_;
    my $out = $this->{fh};
    my $num = exists($this->{summary}->{$key}) ? $this->{summary}->{$key} : 0;
    $denom =~ s/,//g if defined($denom);
    my $pct = $denom ? (int($num/$denom*1000+0.5)/10).'%' : '';
    if ( $num )
    {
        my $tmp = reverse $num;
        $tmp =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
        $num = scalar reverse $tmp;
    }
    print $out $key, "\t", $num, "\t", $pct, "\n";
    return $num;
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

# TODO
# CHECK FOR OUTLIERS
#    for (my $i=0; $i<=$#summary; $i++)
#    {
#        my $label = $summary[$i];
#        my @stats = ();
#        for (my $j=0; $j<=$#allStats; $j++)
#        {
#            push @stats, $allStats[$j]->[$i];
#        }
#        my ($stdev, $mean) = iTagger::Stats::stdDev(\@stats);
#        for (my $j=0; $j<=$#allStats; $j++)
#        {
#            my $name = $libs[$j];
#            my $x = $allStats[$j]->[$i];
#            my $d = int( abs($mean-$x)/$stdev * 10 + 0.5)/10;
#            if ( $d > 2 )
#            {
#                if ( $x < $mean )
#                {
#                    $this->info("Lib $name, $label=$x is $d stdev lt mean");
#                } else
#                {
#                    $this->info("Lib $name, $label=$x is $d stdev gt mean");
#                }
#            }
#        }
#    }

1;

=back

=head1 AUTHORS

Edward Kirton, Julien Tremblay

=head1 COPYRIGHT

Copyright (c) 2013 US DOE Joint Genome Institute.  Use freely under the same license as Perl itself.  Refer to duk and flash documentation for their own copyright/license information.

=cut
