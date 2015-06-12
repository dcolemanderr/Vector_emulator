=pod

=head1 NAME

iTagger::Stats

=head1 DESCRIPTION

Flyweight stats module.

=head1 FUNCTIONS

=over 5

=cut

package iTagger::Stats;

use strict;
use warnings;
require Exporter;

our $VERSION = 1.0;
our @ISA = qw(Exporter);
our @EXPORT = qw(sum meanMedian mean oddMedian mode minMax stdDev var trimmedStdDev round);

=item sum

Calc sum of a list.

=cut

sub sum
{
    my $ar = shift;
    my $sum = 0;
    foreach my $x (@$ar) { $sum += $x }
    return $sum;
}

=item meanMedian

Calc mean-median of list.

=cut

sub meanMedian
{
	my $ar = shift;
	my @a = sort { $a <=> $b } @$ar;
	if (@a % 2)
	{
		return $a[ @a / 2 ];
	} else
	{
		return ($a[ @a / 2 - 1 ] + $a[ @a / 2 ]) / 2;
	}
}

=item mean

Calc mean of list

=cut

sub mean
{
	my ($ar) = @_;
	my $result;
	foreach (@$ar) { $result += $_ }
	return $result / scalar @$ar;
}

=item oddMedian

Calc odd-median of list

=cut

sub oddMedian
{
	my $ar = shift;
	my @a  = sort @$ar;
	return $a[ (@a - (0, 0, 1, 0)[ @a & 3 ]) / 2 ];
}

=item mode

Calc mode of list

=cut

sub mode
{
	my $ar = shift;
	my (%count, @result);

	foreach (@$ar) { $count{$_}++ }

	foreach (sort { $count{$b} <=> $count{$a} } keys %count)
	{
		last if @result && $count{$_} != $count{$result[0]};
		push(@result, $_);
	}

	# Uncomment the following line to return undef for nonunique modes.
	# return undef if @result > 1;

	# Return the odd median of the modes.
	return oddMedian \@result;    # odd_median() is defined earlier.
}

=item minMax

Calc min and max values in list.

=cut

sub minMax
{
	my $ar = shift;
	return (undef, undef) unless scalar(@$ar) > 0;
	my $min = $ar->[0];
	my $max = $ar->[0];
	foreach my $x (@$ar)
	{
		$min = $x if $x < $min;
		$max = $x if $x > $max;
	}
	return ($min, $max);
}

=item stdDev

Calc standard deviation and mean of list

=cut

sub stdDev
{
	my ($var, $mean) = var(@_);
    return (undef, undef) unless defined($mean);
    return (undef, $mean) unless defined($var);
	return (sqrt($var), $mean);
}

=item var

Calculate variance and mean of list.

=cut

sub var
{
	my $ar = shift;
	my $n  = scalar(@$ar);
	return undef unless $n > 1;
	my $mean = mean($ar);
	my $var  = 0;
	foreach my $x (@$ar) { $var += ($x - $mean)**2 }
	$var /= ($n - 1);
	return ($var, $mean);
}

=item trimmedStdDev

Calculate standard deviation and mean of trimmed list of values.
Trim parameter is fraction of low and high values discarded; for example, a
(default) trim parameter of 0.10 results in the smallest 10% and largest 10% of
values being discarded, for a total of 20% of the values discarded.
Note that a small list and high trim parameter may result in too few values to
calculate standard deviation, resulting in an undefined result.

=cut

sub trimmedStdDev
{
	my ($ar, $x) = @_;
	$x = 0.10 unless defined($x) and $x > 0;
	my $n = scalar(@$ar);
	my $t = int($x * $n + 0.5);
	return stdDev($ar) unless $t > 0;
	my @a = sort { $a <=> $b } @$ar;
	@a = @a[ $t .. ($n - $t - 1) ];
	return stdDev(\@a);
}

=item round

Round an item to desired number of decimal places (default=0)

=cut

sub round
{
	my ($x,$n)=@_;
	$n = 0 unless $n;    # OPTIONAL, NUMBER OF DECIMAL PLACES DESIRED
	return (int($x * 10**$n + .5) / 10**$n);
}

1;

=back

=cut
