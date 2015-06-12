=pod

=head1 NAME

iTagger::Logger

=head1 DESCRIPTION

Logger class.  Currently wraps log4perl.  Persistent connections are not used because when many concurrent jobs are run, we may hit the db connection limit.  Consequently, connections are created for each logging event, so be brief.

=head1 METHODS

=cut

package iTagger::Logger;

use strict;
use warnings;
use Env qw(ITAGGER_LOG_CONFIG);
use Log::Log4perl;

our $VERSION = 1.3;
our @ISA = qw(Exporter);
our @EXPORT = qw(new logdie logwarn error info debug);

=head2 CONSTRUCTOR AND DESTRUCTOR

=over 5

=item new

Constructor.  Config is either a hashref or the path to an .ini format file.

=cut

sub new
{
    my ($class, $loggerName, $loggerConfigFile) = @_;
    logdie("Logger: logger name required\n") unless $loggerName;
    if ( $loggerConfigFile )
    {
        Log::Log4perl->init($loggerConfigFile);
    } elsif ( $ITAGGER_LOG_CONFIG )
    {
        Log::Log4perl->init($ITAGGER_LOG_CONFIG);
    } else
    {
        logdie("Logger: either log config file parameter or environment variable ITAGGER_LOG_CONFIG is required\n");
    }
    my $logger = Log::Log4perl->get_logger($loggerName);
    my $this =
    {
        loggerName => $loggerName,
        is => 
        {
            logdie => 1,
            trace => $logger->is_trace(),
            debug => $logger->is_debug(),
            info => $logger->is_info(),
            logwarn => $logger->is_warn(),
            error => $logger->is_error(),
            fatal => $logger->is_fatal()
        }
    };
    bless $this, $class;
    Log::Log4perl->remove_logger($logger);
    return $this;
}

=item _log

Private method handles all logging events.

=cut

sub _log
{
    my ($this, $type, $msg) = @_;
    return unless $this->{is}->{$type}; # avoid connecting to db if not going to insert entry
    if ( !defined($msg) ) { $msg = '' }
    elsif ( $msg =~ /^(.*)\n+$/ ) { $msg = $1 }
    my $logger = Log::Log4perl->get_logger($this->{loggerName});
    if ( $type eq 'logdie' ) { $logger->logdie($msg) }
    elsif ( $type eq 'logwarn' ) { $logger->logwarn($msg) }
    elsif ( $type eq 'error' ) { $logger->error_logdie($msg) }
    elsif ( $type eq 'info' ) { $logger->info($msg) }
    elsif ( $type eq 'debug' ) { $logger->debug($msg) }
    else { $logger->logdie("Invalid log type, '$type', for message: $msg\n") }
    Log::Log4perl->remove_logger($logger);
}

=item logdie

Log fatal message and logdie.

=cut

sub logdie { shift->_log('logdie', @_) }

=item logwarn

Log fatal message and logdie.

=cut

sub logwarn { shift->_log('logwarn', @_) }

=item error

Log error message.

=cut

sub error { shift->_log('error', @_) }

=item info

Log info message.

=cut

sub info { shift->_log('info', @_) }

=item debug

Log debug message.

=cut

sub debug { shift->_log('debug', @_) }

=item trace

Log trace message.

=cut

sub trace { shift->_log('trace', @_) }

1;

=back

=head1 AUTHORS

Edward Kirton, Julien Tremblay

=head1 COPYRIGHT

Copyright (c) 2013 US DOE Joint Genome Institute.  Use freely under the same license as Perl itself.  Refer to duk and flash documentation for their own copyright/license information.

=cut
