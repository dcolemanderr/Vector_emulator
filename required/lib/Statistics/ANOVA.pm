package Statistics::ANOVA;

use 5.008008;
use strict;
use warnings;
use Carp qw(croak carp);
use Algorithm::Combinatorics qw(combinations partitions);
use Math::Cephes qw(:dists);
use Math::Cephes::Matrix qw(mat);
use Scalar::Util qw(looks_like_number);
use Statistics::Descriptive;
use Statistics::Lite qw(mean sum max min);
use vars qw($VERSION);
$VERSION = '0.08';

=head1 NAME

Statistics::ANOVA - Parametric and nonparametric 1-way analyses of variance for means-comparison and clustering per differences/trends over independent or repeated measures of groups or levels

=head1 VERSION

This is documentation for B<Version 0.08> of Statistics::ANOVA.

=head1 SYNOPSIS

 use Statistics::ANOVA 0.08;
 my $aov = Statistics::ANOVA->new();

 # Some data:
 my @gp1 = (qw/8 7 11 14 9/);
 my @gp2 = (qw/11 9 8 11 13/);

 # Load the data (names should be stringy or numerical, depending on level of measurement):
 $aov->load_data({1 => \@gp1, 2 => \@gp2}); # NB: hashref
 # or $aov->load_data([ [1, \@gp1], [2, \@gp2] ]);
 # or $aov->load_data([ [1, @gp1], [2, @gp2] ]);

 # Now here comes another one - too late for just comparing the two:
 my @gp3 = (qw/7 13 12 8 10/);
 $aov->add_data(3 => \@gp3);
 
 # All loaded as Statistics::Descriptive objects in hash named "data", so, e.g.,
 printf("Mean for Group $_ = %.2f\n", $aov->{'data'}->{$_}->mean) foreach sort {$a cmp $b} keys %{$aov->{'data'}};

 #  Test equality of variances before omnibus comparison:
 $aov->obrien()->dump(title => 'O\'Brien\'s test of equality of variances');
 $aov->levene()->dump(title => 'Levene\'s test of equality of variances');

 # 1.10 Independent nominal groups ANOVA - parametric testing:
 $aov->anova(independent => 1, parametric => 1)->dump(title => 'Indep. groups parametric ANOVA', eta_squared => 1, omega_squared => 1);
 # 1.11 Independent nominal groups ANOVA - NON-parametric:
 $aov->anova(independent => 1, parametric => 0)->dump(title => 'Kruskal-Wallis test');

 #  or if independent AND ordered groups/levels: test linear/non-linear trend:
 # 1.20 Independent ordinal groups ANOVA - parametric testing:
 $aov->anova(independent => 1, parametric => 1, ordinal => 1)->dump(title => 'Indep. groups parametric ANOVA: Linear trend');
 $aov->anova(independent => 1, parametric => 1, ordinal => -1)->dump(title => 'Indep. groups parametric ANOVA: Non-linear trend');
 # 1.21 Independent ordinal groups ANOVA - NONparametric testing:
 $aov->anova(independent => 1, parametric => 0, ordinal => 1)->dump(title => 'Jonckheere-Terpstra test');
 
 #  If they are repeated measures:
 # 2.10 Dependent nominal groups ANOVA - parametric testing:
 $aov->anova(independent => 0, parametric => 1)->dump(title => 'Dependent groups ANOVA');
 # 2.11 Dependent nominal groups ANOVA - NONparametric testing:
 $aov->anova(independent => 0, parametric => 0, f_equiv => 0)->dump(title => 'Friedman test');
 
 # or if repeated AND ordinal measures:
 # 2.20 Dependent ordinal groups ANOVA - parametric testing: NOT yet IMPLEMENTED
 #$aov->anova(independent => 0, parametric => 1)->dump(title => '');
 # 2.21 Dependent ordinal groups test - NONparametric testing:
 $aov->anova(independent => 0, parametric => 0, ordinal => 1, f_equiv => 0)->dump(title => 'Page test');

 # Get pairwise comparisons (nominality of the factor assumed):
 $aov->compare(independent => 1, parametric => 1, flag => 1, alpha => .05, dump => 1); # Indep. obs. F- (or t-)tests
 $aov->compare(independent => 0, parametric => 1, flag => 1, alpha => .05, dump => 1); # Paired obs. F (or t-)tests
 $aov->compare(independent => 1, parametric => 0, flag => 1, alpha => .05, dump => 1); # Wilcoxon (between-groups) sum-of-ranks (Dwass Procedure)
 $aov->compare(independent => 0, parametric => 0, flag => 1, alpha => .05, dump => 1); # Friedman-type (within-groups) sum-of-ranks
 
 # Cluster analysis: the least variant of all possible first-order binary-splits of the data (indep. assumed):
 $aov->cluster(parametric => 1, dump => 1); # Scott-Knott method (default)
 $aov->cluster(parametric => 0, dump => 1); # Worsley method

 print $aov->table(precision_p => 3, precision_s => 3);
 
 $aov->unload('g3'); # back to 2 datasets (g1 and g2) - each Statistics::Descriptive objects

=head1 DESCRIPTION

By setting the extended Boolean (0, 1 or -1) value of only three parameters (C<independent>, C<parametric> and C<ordinal>), this module gives you oneway parametric or non-parametric analyses-of-variance (ANOVAs) for either nominal groups or ordinal levels (trend analysis), and for either independent or dependent (repeated measures) observations within each group/level. For the Fisher-esque ANOVAs, you can also access estimates of proportion of variance accounted for (I<eta>-squared) and effect-size (I<omega>-squared), plus I<a priori> pairwise comparisons by the relevant independent or dependent I<t>-tests. Non-Fisher-esque, non-parametric tests comprise the Kruskal-Wallis, Friedman and Page tests, all with default accounting for ties in the calculation of ranks, and standardizing of the test-statistics. Simple parametric and non-parametric I<post hoc> clustering is also offered (Scott-Knott and Worsley methods). The module also provides for testing equality of variances (O'Brien and Levene tests) across independent groups.

A basic design principle has been to offer as few method calls as possible, and to steer queries into the proper underlying method by Boolean manipulation of a minimal parameter set. To get to any relevant test of your data, there are (only) three Boolean vars to set, and these are named adjectivally: C<independent>, C<parametric> and C<ordinal>. So the relevant tests are offered I<not> in the form of little packages specific to each test, but by specifying the attributes of the data: are they based on independent groups or repeated measures? do the groups/measures form independent categories or related levels? are relationships, likenesses or differences to be found? do the data support parametric testing? It seems more useful for programmed access to statistical algorithms to be initially sensitive to the answers to these and like questions, while keeping the methods of access and set of arguments the same across tests, and the actual tests relatively invisible, rather than bearing the statistical method up-front, each package offering a unique interface and algorithm for the same things, each dependent on a unique set of arguments and their correct norder. It could be argued that that's what an application should do, not a module. This module goes the other way, recognizing that the algorithms share arguments and internal methods, which is typically obscured by algorithm-specific packaging.

Reliability has been assayed by testing each method with at least two different published sources; and comparing the output with one or another open-source or commercial statistics package. The tests based on published examples are implemented during cpan-wise installation; see the "t" folder of L<cpan.org|www.cpan.org>'s installation-distribution of this module. News of unreliabilities are welcome; a fundamental one from Cathal Seoghie has already helped, namely, by pointing to the need to account for NaNs, empty, and invalid values.

=head1 METHODS

=cut

=head2 INTERFACE

Object-oriented. No subs are explicitly exported, no parameters are set for cross-method application. 

The class-object yet always holds the myriad of statistics produced by the last I<test> you called - i.e., it is fed with relevant values upon each method-call; a "p_value" at least. 

The methods generally only return the class object itself. 

Most tests set a common lot of statistical I<keys>(e.g., denoting between- and within-groups sums-of-squares), but a few are idiosyncratic per test. 

Meanwhile, just L<dump|dump> as much as you can, and bless your digestion.

=head3 new

 $aov = Statistics::ANOVA->new()

Create a new Statistics::ANOVA object for accessing the subs.

=cut

#-----------------------------------------------
sub new {
#-----------------------------------------------
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self= {};
    bless($self, $class);
    return $self;
}

=head2 HANDLING DATA

=head3 load

 $aov->load('aname', @data1)
 $aov->load('aname', \@data1)
 $aov->load(['aname', @data1], ['another_name', @data2])
 $aov->load(['aname', \@data1], ['another_name', \@data2])
 $aov->load({'aname' => \@data1, 'another_name' => \@data2})

I<Alias>: C<load_data>

Accepts data for analysis in any of the above-shown forms, but always with the requirement that:

=over

=item 1.

a single set of observations (the "group" or "level") is given a unique name, and 

=item 2.

you do not mix the methods, e.g., send a bit of a hashref, plus a bit of an aref. That's perverse. 

=back

The reason for all these options is that there are so many to be found in Perl's Statistics modules that it would be a pain for you to have to re-architecture your data just to use this module. So, namely, you can load:

=over

=item 1. 

a single C<name =E<gt> value> pair of a sample name, and a list (referenced or not) of data. This is a simple but useless load: Presumably, you will follow this up, when you're ready, with a call to L<add|add> so that there's actually something relevant to test by this module. Don't try another call to C<load> with new data because that will only clobber what you've just done.

=item 2.

a reference to an array of (referenced) arrays, where each of the latter arrays consists of a sample name occupying the first index, and then its sample data, as an array or yet another referenced array.

=item 3. 

a hash reference of named array references of data. This is the preferred method - the one that is first checked in the elongated C<if> clause that parses all this variety.

=back

The data are loaded into the class object by name, within a hash named C<data>, as L<Statistics::Descriptive::Full|Statistics::Descriptive> objects. So you can easily get at any descriptives for the groups you've loaded - e.g., $aov->{'data'}->{'aname'}->mean() - or you could get at the data again by going $aov->{'data'}->{'aname'}->get_data() - and so on; see L<Statistics::Descriptive> for the full Christmas tree of possibilities available to you.

The names of the data are up to you; whatever can be set as the key in a hash. But if you intend to do trend analysis, you should, as a rule, give only I<numerical> names to your groups/levels, at least defining their ordinality.

Each call L<unload|unload>s any previous loads.

Returns the Statistics::ANOVA object - nothing but its blessed self.

=cut

#-----------------------------------------------
sub load {
#-----------------------------------------------
    my $self = shift;
    $self->unload();
    $self->add(@_);
}
*load_data = \&load; # Alias

=head3 add

 $aov->add('another_name', @data2)
 $aov->add('another_name', \@data2)
 $aov->add(['another_name', @data2])
 $aov->add(['another_name', \@data2])
 $aov->add({'another_name' => \@data2})

I<Alias>: C<add_data>

Same as L<load|load> except that any previous loads are not L<unload|unload>ed. Again, the hash-referenced list is given preferential treatment.

=cut

#-----------------------------------------------
sub add {
#-----------------------------------------------
    my $self = shift;
    my ($name, $data) = ();

    if (ref $_[0] eq 'HASH') {
      while (($name, $data) = each %{$_[0]}) {
         if (ref $data) {
              $self->{'data'}->{$name} = Statistics::Descriptive::Full->new();
              $self->_mark_purge($name, $data);
              $self->{'data'}->{$name}->add_data($data);
         } 
      }
    }
    elsif (ref $_[0] eq 'ARRAY') {
        $self->add(_aref2href($_[0]));
    }
    else {
       $name = shift;
       $data = ref $_[0] eq 'ARRAY' ? $_[0] : scalar (@_) ? \@_ : croak 'No list of data for ANOVA';
       $self->{'data'}->{$name} = Statistics::Descriptive::Full->new();
       $self->_mark_purge($name, $data);
       $self->{'data'}->{$name}->add_data($data);
    }
}
*add_data = \&add; # Alias

=head3 unload

 $aov->unload()     # bye-bye everything
 $aov->unload('g1') # so long data named "g1"

I<Alias>: C<delete_data>

With nil or no known arguments, empties all cached data and calculations upon them, ensuring these will not be used for testing. This will be automatically called with each new load, but, to take care of any development, it could be good practice to call it yourself whenever switching from one dataset for testing to another.

Alternatively, supply one or more names of already loaded data to clobber just them out of existence; preserving any other loads.

=cut

#-----------------------------------------------
sub unload {
#-----------------------------------------------
    my ($self)  = shift;
    if (scalar(@_)) {
        foreach (@_) {
            delete $self->{'data'}->{$_};
            delete $self->{'_purge'}->{$_};
        }
    }
    else { # delete everything
        $self->{'data'} = {};
    }
    $self->{$_} = undef foreach qw/_eq_var purged/;
    #$self->{'_stat'}->{$_} = undef foreach keys %{$self->{'_stat'}};
    $self->{'_cleared'} = 1;
}
*delete_data = \&unload; # Alias

=head3 I<Missing/Invalid values>

Any data-points/observations sent to L<load|load> or L<add|add> that are undefined or not-a-number are marked for purging before being anova-tested or tested pairwise. The data arrays accessed as above, as L<Statistics::Descriptive::Full|Statistics::Descriptive>, will still show the original values. When, however, you call one of the anova or pairwise methods, the data must and will be purged of these invalid values before testing. 

When the C<independent> parameter equals 1 when sent to L<anova|anova>, L<compare|compare> or L<cluster|cluster>, each list is simply purged of any undefined or invalid values. This also occurs for the equality of variances tests.

When C<independent> parameter equals 0 when sent to L<anova|anova>, L<compare|compare> and L<cluster|cluster>, each list is purged of any value at all indices that, in any list, contain invalid values. So if two lists are (1, 4, 2) and (2, ' ', 3), the lists will have to become (1, 2) and (2, 3) to account for the bung value in the second list, and to keep all the observations appropriately paired.

The number of indices that were subject to purging is cached thus: $aov->{'purged'}. The L<dump|dump> method can also reveal this value. 

The C<looks_like_number> method in L<Scalar::Util|Scalar::Util/looks_like_number> is used for checking validity of values. (Although Params::Classify::is_number might be stricter, looks_like_number benchmarks at least a few thousand %s faster.) 

=cut

=head2 MEASURING SIGNIFICANCE

One generic method L<anova|anova> (a.k.a. aov, test) is used to access parametric or nonparametric tests for independent or dependent/related observations, and categorical prediction or trend analysis. Accessing the different statistical tests depends on setting I<three> parameters on a true/false basis: I<independent>, I<parametric> and I<ordinal>. The following describes the particular tests you get upon each possible combination of these alternatives.

=head3 1. INDEPENDENT groups/levels

=head4 1.10 PARAMETRIC test for NOMINAL groups

 $aov->anova(independent => 1, parametric => 1, ordinal => 0)

Offers the standard Fisher-esque ANOVA for a single factor applied to independent groups of data-sources (people, plots-of-land, cultures, etc.). The attribute C<independent> refers to whether or not each level of the independent variable was yielded by independent or related sources of data; e.g., If the same people provided you with ratings on the levels of "What pacifies me" and "What freaks me out," that would be a dependent, repeated measures situation; i.e., the same people, or plot-of-land, whatever, contributed an observation under each level of the predictive, factor (or the "independent variable", or "treatment"). In the latter case, C<independent> must equal zero. But if a unique set of people, or plot-of-land, or culture, yielded the measures over the predictive factor, then C<independent> = 1 - each data source per group/level is not derived from the same data-source.

=cut

#-----------------------------------------------
sub _aov_indep_param_cat {
#-----------------------------------------------
    my ($self, $data, %args) = @_;

    my ($ss_b, $df_w, $ss_w, $grand_mean, $df_b, $ms_b, $ms_w, $f, $f_prob, $desc) = ();
    
    ($ss_w, $df_w, $desc, $grand_mean) = _indep_error($data);# within sum-of-squares (ss_w) & error dfs (df_w) (et al.):

    $ss_b += $_->[0] * ($_->[1] - $grand_mean)**2 foreach @{$desc}; # between sum-of-squares; $desc is aref of arefs of counts [0] & means [1] per gp
    $df_b = scalar(keys(%{$data})) - 1;
    $ms_b = $ss_b / $df_b;
    $ms_w = $ss_w / $df_w;
    $f = $ms_b / $ms_w;
    $f_prob = fdtrc($df_b, $df_w, $f);

    # (Benchmarking shows it is 3%-12% faster to calculate these vars firstly into scalars and then feed them into the object, rather than calc them directly into the object, even for the simplest arithmetical operations)
    ($self->{'_stat'}->{'f_value'}, $self->{'_stat'}->{'p_value'}, $self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}, $self->{'_stat'}->{'ss_b'},	$self->{'_stat'}->{'ss_w'}, $self->{'_stat'}->{'ms_b'}, $self->{'_stat'}->{'ms_w'}) = ($f, $f_prob, $df_b, $df_w, $ss_b, $ss_w, $ms_b, $ms_w);
    $self->{'_dfree'} = 0;
}

=head4 1.11 PARAMETRIC test for ORDINAL levels

 $aov->anova(independent => 1, parametric => 1, ordinal => 1) # test linear trend
 $aov->anova(independent => 1, parametric => 1, ordinal => -1) # test non-linear trend

If the independent/treatment/between groups variable is actually measured on a continuous scale/is a quantitative factor, assess their B<linear trend>: Instead of asking "How sure can we be that the means-per-group are equal?", ask "How sure can we be that there is a departure from flatness of the means-per-level?". 

The essential difference is that in place of the the between (treatment) mean squares in the numerator is the linear sum of squares in which each "group" mean is weighted by the deviation of the level-value (the name of the "group") from the mean of the levels (and divided by the sum of the squares of these deviations).

If the number of observations per level is unequal, the module applies the simple I<unweighted> approach. This is recommended as a general rule by Maxwell and Delaney (1990), given that the I<weighted> approach might erroneously suggest a linear trend (unequal means) when, in fact, the trend is curvilinear (and by which the means balance out to equality); unless "there are strong theoretical reasons to believe that the only true population trend is linear" (p. 234). (But then you might be theoretically open to either. While remaining as the default, a future option might access the I<hierarchical, weighted> approach.)

To test if there is the possibility of a B<non-linear trend>, give the value of -1 to the C<ordinal> argument.

Note that the contrast coefficients are calculated directly from the values of the independent variable, rather than using a look-up table. This respects the actual distance between values, but requires that the names of the sample data, of the groups (or levels), have been given I<numerical> names when L<load|load>ed - i.e., such that the data-keys can be summed and averaged.

=cut

#-----------------------------------------------
sub _aov_indep_param_ord {
#-----------------------------------------------
    my ($self, $data, %args) = @_;

    my ($f, $levels, $ss_l, $ms_w, $ss_w, $df_b, $df_w, $desc, $grand_mean, $num) = ();
    
    $levels = scalar(keys(%{$data}));
    $ss_l = _indep_ss_ord($data);
    ($ss_w, $df_w, $desc, $grand_mean) = _indep_error($data);
    $ms_w = $ss_w / $df_w;
    
    if ($args{'ordinal'} == 1) { # don't need to calculate ss_b
        $num = $ss_l;
        $df_b = $levels - 1;
    }
    elsif ($args{'ordinal'} == -1)  { # calc treatment sum-of-squares (ss_b):
        my $ss_b;
        $ss_b += $_->[0] * ($_->[1] - $grand_mean)**2 foreach @{$desc}; # $desc is aref of arefs of counts [0] & means [1] per gp
        $num = ($ss_b - $ss_l) / ($levels - 2);
        $df_b = $levels - 2;
        $self->{'_stat'}->{'ss_b'} = $ss_b;
    }

    $f = $num / $ms_w;
    $self->{'_stat'}->{'p_value'} = fdtrc($df_b, $df_w, $f); # Math::Cephes function
    $self->{'_stat'}->{'ss_l'} = $ss_l;
    ($self->{'_stat'}->{'f_value'}, $self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}, $self->{'_stat'}->{'ss_w'}, $self->{'_stat'}->{'ms_w'}) = ($f, $df_b, $df_w, $ss_w, $ms_w);
    $self->{'_dfree'} = 0;
}

=head4 1.20 NONPARAMETRIC test for NOMINAL groups (Kruskal-Wallis test)

 $aov->anova(independent => 1, parametric => 0, ordinal => 0)

Performs a one-way independent groups ANOVA using the non-parametric B<Kruskal-Wallis> sum-of-ranks method for 3 or more groups of a single factor. Instead of an I<F>-value, there is now a I<H>-value. The I<p>-value is read off the chi-square distribution; note that this is unreliable when there are no more than 3 groups and all groups comprise 5 or fewer observations.

By default, this method accounts for and corrects for ties, but if C<correct_ties> = 0, I<H> is uncorrected. The correction involves giving each tied score the mean of the ranks for which it is tied (see Siegal, 1956, p. 188ff).

=cut

#-----------------------------------------------
sub _aov_indep_dfree_cat {
#-----------------------------------------------
    my ($self, $data, %args) = @_;

    my ($h, $df_b, $n) = _kw_stat($data, $args{'correct_ties'});
    if ($args{'f_equiv'}) {
        my $levels = scalar(keys(%{$data}));
        my $f = ( $h / ($levels - 1) ) / ( ($n - 1 - $h) / ($n - $levels) );
        my $df_w;
        $df_w += ($data->{$_}->count() - 1) foreach keys %{$data};
        $self->{'_stat'}->{'p_value'} = fdtrc($df_b, $df_w, $f);
        ($self->{'_stat'}->{'f_value'}, $self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}) = ($f, $df_b, $df_w);# Feeding-time
        $self->{'_dfree'} = 0;
    }
    else {
        $self->{'_stat'}->{'p_value'} = chdtrc($df_b, $h);
        ($self->{'_stat'}->{'h_value'}, $self->{'_stat'}->{'df_b'}) = ($h, $df_b);# Feeding-time
        $self->{'_dfree'} = 1;
    }
    return $self;
}

=head4 1.21 NONPARAMETRIC test for ORDINAL levels (Jonckheere-Terpstra test)

 $aov->anova(independent => 1, parametric => 0, ordinal => 1)

Performs the B<Jonckheere-Terpstra> test that respects the given, numerical order of the levels of the independent/treatment variable; unlike the Kruskal-Wallis test that will return the same result regardless of the order of levels. Its test statistic I<J> is the sum of I<k>(I<k> - 1)/2 Mann-Whitney I<U> counts (sum of 1/0 for min/max value of all possible pairs of observations, between groups). Rather than calculating the exact I<p>-value, the present implementation calculates an expected I<J> value and variance (sensitive to tied ranks), to provide a normalized I<J> for which the I<p>-value is read off the normal distribution. This is appropriate for "large" samples, e.g., greater-than 3 levels, with more than eight observations per level. Otherwise, read the value of C<$aov-E<gt>{'j_value'}> and look it up in a table of I<j>-values, such as in Hollander & Wolfe (1999), p. 649ff. The class object is fed:

 $aov->{'j_value'}   :  the observed value of J
 $aov->{'j_exp'}     :  the expected value of J
 $aov->{'j_var'}     :  the variance of J
 $aov->{'z_value'}   :  the normalized value of J
 $aov->{'p_value'}   :  the one-tailed probability of observing a value as great as or greater than z_value.

By default, the method accounts for and corrects for ties, but if C<correct_ties> = 0, I<j_var> is the usual "null" distribution variance, otherwise with an elaborate correction accounting for the number of tied "groups" and each of their sizes, as offered by Hollander & Wolfe (1999) Eq 6.19, p. 204.

=cut

#-----------------------------------------------
sub _aov_indep_dfree_ord {
#-----------------------------------------------
    my ($self, $data, %args) = @_;
    my ($n, $nj, $grand_n, $j_value, $exp, $term_v, $var, $g1, $g2, $g3, $pairs) = (0, 0);
    
    # Get between-group rankings for all possible pairwise splits of the data, accumulating J-statistic:
    my @all_pairs = combinations([keys %{$data}], 2); # Algorithm::Combinatorics function
    foreach $pairs (@all_pairs) {
        my ($ranks_href) = _ranks_between({$pairs->[0] => $data->{$pairs->[0]}, $pairs->[1] => $data->{$pairs->[1]}});
        # Use the ranks to calculate the pair of Mann-Whitney Us for each paired group:
        my ($n1, $n2) = ($data->{$pairs->[0]}->count, $data->{$pairs->[1]}->count); # get Ns per group in the pair
        my $nprod = $n1 * $n2;
        my ($n, @us) = ();
        foreach (keys %{$ranks_href}) {
            $n = $data->{$_}->count;
            push @us, $nprod + (($n * ($n + 1))/2) - sum(@{$ranks_href->{$_}});
        }
        $j_value += max(@us); # Accumulate J-statistic with the maximum of the two U-values
    }
    
    # Calc expected J and variance:
    my @ns = ();
    foreach (keys %{$data}) {
        $n = $data->{$_}->count;
        $grand_n += $n;
        $nj += $n**2;
        $term_v += $n**2 * ( 2 * $n + 3 );
        push @ns, $n;
    }
    $exp = ( $grand_n**2 - $nj ) / 4;

    if (defined $args{'correct_ties'} and $args{'correct_ties'} == 0) {
        $var = ( $grand_n**2 * ( 2 * $grand_n + 3) - $term_v ) / 72;
    }
    else {
        my $ng = scalar(keys(%{$data}));
        my ($uranks_href, $uties_var, $gn, $xtied) = _ranks_between($data);
        my ($term_a2, $term_a3, $term_b2, $term_b3, $term_c2, $term_c3) = ();
        $term_a2 += $_ * ($_ - 1) * (2 * $_ + 5) foreach @ns;
        $term_a3 += $_ * ($_ - 1) * (2 * $_ + 5) foreach @$xtied;
        $term_b2 += $_ * ($_ - 1) * ($_ -2) foreach @ns;
        $term_b3 += $_ * ($_ - 1) * ($_ - 2) foreach @$xtied;
        $term_c2 += $_ * ($_ - 1) foreach @ns;
        $term_c3 += $_ * ($_ - 1) foreach @$xtied;

        $var = ( 1/72 * ($gn * ($gn - 1) * (2 * $gn + 5) - $term_a2 - $term_a3) ) +
            ( 1/( 36 * $gn * ( $gn - 1 ) * ( $gn - 2 ) ) * $term_b2 * $term_b3 ) +
            ( 1/( 8 * $gn * ( $gn - 1 ) ) * $term_c2 * $term_c3);
    }
    
    my $z_value = ($j_value - $exp) / sqrt($var);
    my $p_value = (1 - ndtr(abs($z_value))); # Math::Cephes fn
    $p_value *= 2 unless $args{'tails'} and $args{'tails'} == 1;
    
    ($self->{'_stat'}->{'j_value'}, $self->{'_stat'}->{'j_exp'}, $self->{'_stat'}->{'j_var'}, $self->{'_stat'}->{'z_value'}, $self->{'_stat'}->{'p_value'}) = ($j_value, $exp, $var, $z_value, $p_value);

    return $self;    
}

=head3 2. DEPENDENT groups/levels (REPEATED MEASURES)

=head4 2.10 PARAMETRIC test for NOMINAL groups

 $aov->anova(independent => 0, parametric => 1, ordinal => 0, multivariate => 0|1)

The parametric repeated measures analysis of variance is performed. By default, this is performed using the traditional univariate, or "mixed-model," approach, with sphericity assumed (i.e., equal variances of all factor differences, within each factor and all possible pairs of factors). The assumption is met when there are only two levels of the repeated measures factor; but unequal variances might be a problem when there are more than two levels. 

[TO DO: In order to account for the possibility of violated sphericity, two strategies are typically used: either adjustment of the degrees-of-freedom by one or another method (e.g., Huynh-Feldt procedure) or a multivariate analysis of variance. Presently, the module permits (only) the latter option, given that it can also be used in the next stage for comparing individual means (whereas adjustments of the degrees-of-freedom in paired comparisons is not recommended) (Maxwell & Delaney, 1992, Ch. 11). In order to perform a multivariate analysis of variance, simply specify the parameter C<multivariate> and give it a value of 1. Alternatively, call C<anova_mv>.  In fact, the multivariate approach is recommended as the default procedure, which might be implemented in a future version.]

=cut

sub _aov_rmdep_cat_param {
    my ($self, $data, $samples, $levels, %args) = @_;
    my ($ms_b, $ms_w, $f, $f_prob, $ss_w, $df_w, $ss_b, $df_b) = ();

	# TO DO - multivariate option not yet available:
    if (! $args{'multivariate'}) {
 		($ss_w, $df_w, $ss_b, $df_b) = _rmdep_ss_uni($data, $samples, $levels);# error sum-of-squares (ss_w) & error df (et al.)
    	$ms_b = $ss_b / $df_b;
    	$ms_w = $ss_w / $df_w;
	    $f = $ms_b / $ms_w;    
	}
	else {
		($ss_w, $df_w, $ss_b, $df_b) = _rmdep_ss_multi($data, $samples, $levels);# error sum-of-squares (ss_w) & error df (et al.)
    	#$ms_b = $ss_b / $df_b;
    	#$ms_w = $ss_w / $df_w;
		$f = $ss_b / $ss_w;
	}
    $self->{'_stat'}->{'p_value'} = fdtrc($df_b, $df_w, $f); # Math::Cephes fn
    ($self->{'_stat'}->{'f_value'}, $self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}, $self->{'_stat'}->{'ss_b'},	$self->{'_stat'}->{'ss_w'}, $self->{'_stat'}->{'ms_b'}, $self->{'_stat'}->{'ms_w'}) = ($f, $df_b, $df_w, $ss_b, $ss_w, $ms_b, $ms_w);
    $self->{'_dfree'} = 0;
}

=head4 2.11 PARAMETRIC test for ORDINAL levels

[Not implemented.]

=cut

sub _aov_rmdep_ord_param {
    my ($self) = @_;
    $| = 1;
    carp ':-( Parametric trend analysis for dependent/repeated measures is not implemented';
    ($self->{'_stat'}->{'f_value'}, $self->{'_stat'}->{'p_value'}, $self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}, $self->{'_stat'}->{'ss_b'},	$self->{'_stat'}->{'ss_w'}, $self->{'_stat'}->{'ms_b'}, $self->{'_stat'}->{'ms_w'}) = (-9, -9, -9, -9, -9, -9, -9, -9);
    $self->{'_dfree'} = 0;
    return $self;
}

=head4 2.20 NONPARAMETRIC test for NOMINAL groups (Friedman test)

 $aov->anova(independent => 0, parametric => 0, ordinal => 0)

Performs the B<Friedman> nonparametric analysis of variance - for two or more dependent (matched, related) groups. A ranking procedure is used, but, unlike the case for independent groups, the ranks are taken at each common index of each group/level, i.e., within-groups, given that the observations at each index are given by the same data-source (person, plot, etc.). Some tie-handling code from L<Statistics::RankCorrelation|Statistics::RankCorrelation> is channelled here for this purpose.

The statistical attributes now within the class object (see L<anova|anova>) pertain to this test, e.g., $aov->{'chi_value'} gives the chi-square statistic from the Friedman test; and $aov->{'p_value'} gives the associated I<p>-value (area under the right-side, upper tail of the distribution). There is now no defined 'f_value'.

If I<f_equiv> => 1, then, instead of the I<chi>-value, and I<p>-value read off the I<chi>-square distribution, you get the I<F>-value equivalent, with the I<p>-value read off the I<F>-distribution.

By default, the method accounts for and corrects for ties, but if C<correct_ties> = 0, the test-statistic is uncorrected. The correction involves accounting for the number of tied groups at each index, as per Hollander & Wolfe (1999), Eq. 7.8, p. 274.

Please note that parametric trend analysis for dependent/repeated measures is not (as yet) implemented.

=cut

sub _aov_rmdep_cat_dfree {
    my ($self, $data, $samples, $levels, %args) = @_;
    my ($sum_squared_ranks, $chi, $df) = ();
    # Get column sum of ranks of row values:
    my ($ranks, $xtied) = _ranks_within($data, $samples);
    $sum_squared_ranks += $ranks->{$_}**2 foreach keys %{$ranks};
    
    if(defined $args{'correct_ties'} and $args{'correct_ties'} == 0) {
        $chi = ( 12 / ($samples * $levels * ($levels + 1) ) ) * $sum_squared_ranks - 3 * $samples * ($levels + 1);
    }
    else {
        my $num = 12 * $sum_squared_ranks - 3 * $samples**2 * $levels * ($levels + 1)**2;
        my $term_x;
        $term_x += _sumcubes($xtied->{$_}) - $levels foreach keys %{$xtied};
        my $den = $samples * $levels * ($levels + 1) - (1 / ($levels - 1)) * $term_x;
        $chi = $num / $den;
    }

    $df = $levels - 1;
    if ($args{'f_equiv'}) {
        my ($f, $df_b, $df_w, $f_prob) = ();
        $f = (($samples - 1) * $chi) / ($samples * ($levels - 1) - $chi);
        $df_b = $df;
        $df_w = ($samples - 1) * ($df_b);
        $f_prob = fdtrc($df_b, $df_w, $f); # Math::Cephes fn
        ($self->{'_stat'}->{'f_value'}, $self->{'_stat'}->{'p_value'}, $self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}) = ($f, $f_prob, $df_b, $df_w);# Feeding-time
        $self->{'_dfree'} = 0;
    }
    else {
        ($self->{'_stat'}->{'chi_value'}, $self->{'_stat'}->{'p_value'}, $self->{'_stat'}->{'df_b'}) = ($chi, chdtrc($df, $chi), $df); # Math::Cephes fn
        $self->{'_dfree'} = 1;
    }
    
    sub _sumcubes {
        my($tie_aref, $sum) = @_;
        $sum += $_**3 foreach @{$tie_aref};
        return $sum;
    }

    return $self;
}

=head4 2.21 NONPARAMETRIC test for ORDINAL levels (Page test)

 $aov->anova(independent => 0, parametric => 0, ordinal => 1, tails => 1|2)

This option implements the B<Page> (1963) test; see Hollander and Wolfe (1999, p. 284ff) for a review. Ranks are computed exactly as for the Friedman test, but the ranks are weighted according to the ordinal position of the group/level to which they pertain. Also, the test of significance is based on a standardized value, with the I<p>-value read off the normal distribution. Similarly to the relationship between the Kruskal-Wallis and Jonckheere-Terpstra tests for non-dependent observations, the Friedman test returns the same value regardless of the ordinality of the groups/levels, but the Page test respects - and requires - numerical labels of the groups/levels. The groups are weighted according to their Perl sort { $a <=> $b} order, so be sure to give sort-able names that reflect your pre-experimental, hypothetical ordering of the different treatments/groups.

With only two groups, the test statistic is equivalent to that provided by a B<sign test>.

The statistical attributes now within the class object (see L<anova|anova>) pertain to this test, and are chiefly:

 $aov->{'l_value'} : the observed test statistic (sum of ordered and weighted ranks)
 $aov->{'l_exp'}   : expected value of the test statistic
 $aov->{'l_var'}   : variance of the test statistic (given so many groups and observations)
 $aov->{'z_value'} : the standardized l_value
 $aov->{'p_value'} : the 2-tailed probability associated with the z_value (or 1-tailed if tails => 1).

Hollander and Wolfe (1999) describe how Page's I<L>-statistic is directly related to Spearman's rank-order correlation coefficient (see L<Statistics::RankCorrelation|Statistics::RankCorrelation>). They provide a simple "l2r" transformation, and this is also offered, so, just for this test, you can also read:

 $aov->{'r_value'} : estimate of the Spearman rank-order correlation coefficient
  based on the observed and predicted order of each associated group/level per observation.

=cut

#-----------------------------------------------
sub _aov_rmdep_ord_dfree {
#-----------------------------------------------
    my ($self, $data, $samples, $levels, %args) = @_;
    my ($c, $sum_ranks, $exp, $var, $r_value, $z_value, $p_value) = (0);
    
    # Get column/within-group sum of ranks of row values:
    my ($ranks) = _ranks_within($data, $samples);
    $sum_ranks += ++$c * $ranks->{$_} foreach sort {$a <=> $b} keys %{$ranks};
    
    # Standardize this test statistic:
    $exp = ( $samples * $levels * ($levels + 1)**2 ) / 4;
    $var = ( $samples * $levels**2 * ($levels + 1) * ($levels**2 - 1) ) / 144;
    $z_value = ( $sum_ranks - $exp ) / sqrt($var);
    $p_value = (1 - ndtr(abs($z_value))); # Math::Cephes fn
    $p_value *= 2 unless $args{'tails'} and $args{'tails'} == 1;
    $r_value = _l2r($sum_ranks, $samples, $levels);
    
    # Feed the class object:
    ($self->{'_stat'}->{'l_value'}, $self->{'_stat'}->{'l_exp'}, $self->{'_stat'}->{'l_var'}, $self->{'_stat'}->{'z_value'}, $self->{'_stat'}->{'p_value'}, $self->{'_stat'}->{'r_value'}) = ($sum_ranks, $exp, $var, $z_value, $p_value, $r_value);
    $self->{'_dfree'} = 1;
    ##my $s_value = $r_value / sqrt( (1 - $r_value**2) / ($samples - 2) );
    #my $rp_value = stdtr($samples-2, -1 * abs($s_value)); # Math::Cephes fn # returns the left(!) 1-tailed p_value
    #$rp_value *= 2 if defined $args{'tails'} and $args{'tails'} == 2;
    #$self->{'rp_value'} = $rp_value;
    return $self;
    
    sub _l2r {
        my ($l, $n, $levels) = @_;
        return ( ( ( 12 * $l) / ( $n * $levels * ($levels**2 - 1) ) ) - ( (3 * ($levels + 1)) / ($levels - 1)) ) ;
    }
}

=head3 anova

 $aov->anova(independent => 1|0, parametric => 1|0, ordinal => 0|1)

I<Aliases>: aov, test

Generic method to access all anova functions by specifying TRUE/FALSE values for C<independent>, C<parametric> and C<ordinal>. 

    Independent    Parametric  Ordinal    What you get
    1              1           0          Fisher-esque independent groups ANOVA
    1              1           1          Fisher-esque independent groups ANOVA with trend analysis
    1              0           0          Kruskal-Wallis independent groups ANOVA
    1              0           1          Jonckheere-Terpstra independent groups trend analysis    
    0              1           0          Fisher-esque dependent groups ANOVA (univariate or multivariate)
    0              1           1          (Fisher-esque dependent groups ANOVA with trend analysis; not implemented)
    0              0           0          Friedman's dependent groups ANOVA
    0              0           1          Page's dependent groups trend analysis

All methods return nothing but the class object after feeding it with the relevant statistics, which you can access by name, as follows:

 $aov->{'f_value'} (or $aov->{'chi_value'}, $aov->{'h_value'}, $aov->{'j_value'}, $aov->{'l_value'} and/or $aov->{'z_value'})
 $aov->{'p_value'} : associated with the test statistic
 $aov->{'df_b'} : the between-groups or treatment or numerator degree(s) of freedom
 $aov->{'df_w'} : the within-groups or error or denominator degree(s) of freedom (also given with F-equivalent Friedman test)
 $aov->{'ss_b'} : between-groups or treatment sum of squares
 $aov->{'ss_w'} : within-groups or error sum of squares
 $aov->{'ms_b'} : between-groups or treatment mean squares
 $aov->{'ms_w'} : within-groups or error mean squares

=cut

sub anova {
    my ($self, %args) = @_;
    
    foreach (qw/independent parametric nominal/) { # 'nominal' an option not yet implemented
        $args{$_} = 1 if !defined $args{$_} ;
    }
    
    if (!$self->{'_cleared'}) {
        $self->{$_} = undef foreach qw/df_b df_w f_value chi_value h_value j_value j_exp j_var l_value l_exp l_var z_value p_value ss_b ss_w ms_b ms_w _eq_var purged/;
        $self->{'_cleared'} = 1;
    }
    
    if ($args{'independent'}) {
        if ($args{'parametric'}) {
            if (!$args{'ordinal'}) {
                $self->_aov_indep_param_cat($self->_indep_data($args{'data'}), %args);
            }
            else {
                $self->_aov_indep_param_ord($self->_indep_data($args{'data'}), %args);
            }
        }
        else {
            if (!$args{'ordinal'}) {
                $self->_aov_indep_dfree_cat($self->_indep_data($args{'data'}), %args);
            }
            else {
                $self->_aov_indep_dfree_ord($self->_indep_data($args{'data'}), %args);
            }
        }
    }
    else {
        if ($args{'parametric'}) {
            if (!$args{'ordinal'}) {
                $self->_aov_rmdep_cat_param($self->_rmdep_data($args{'data'}), %args);
            }
            else {
                $self->_aov_rmdep_ord_param(%args);
            }
        }
        else {
            if (!$args{'ordinal'}) {
                $self->_aov_rmdep_cat_dfree($self->_rmdep_data($args{'data'}), %args);
            }
            else {
                $self->_aov_rmdep_ord_dfree($self->_rmdep_data($args{'data'}), %args);
            }
        }
    }
    #$self->{'_cleared'} = 0;

    return $self;
}
*aov = \&anova; # Alias
*test = \&anova; # Alias

=head3 Tests for equality of variances

=head4 obrien

 $aov->obrien()

I<Alias>: obrien_test

Performs B<O'Brien's> (1981) test for equality of variances within each group: based on transforming each observation in relation to its group variance and its deviation from its group mean; and performing an ANOVA on these transformed values (for which the group mean is equal to the variance of the original observations). The procedure is recognised to be robust against violations of normality (unlike I<F>-max) (Maxwell & Delaney, 1990).

The statistical attributes now within the class object (see L<anova|anova>) pertain to this test, e.g., $aov->{'f_value'} gives the I<F>-statistic for O'Brien's Test; and $aov->{'p_value'} gives the I<p>-value associated with the I<F>-statistic for O'Brien's Test.

=cut

#-----------------------------------------------
sub obrien {
#-----------------------------------------------
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to a hash of data for performing ANOVA';
    my $tdata = $self->_purge_in_list(\%data);# List-wise clean-up    

    my ($m, $v, $n, $sname, $sdata, @r, @data) = ();       
    $self->{'obrien'} = {};
    
    # Traverse each sample of data:
    while (($sname, $sdata) = each %{$tdata}){
        # For each group, compute the sample mean and the unbiased sample variance:
        ($m, $v, $n)  = ($sdata->mean, $sdata->variance, $sdata->count);
        # Transform each observation:
        @data = $sdata->get_data;
        foreach (@data) {
            push @r, (
                    # Numerator:
                    (
                        ( ($n - 1.5 ) * $n * (($_ - $m)**2) )
                        -
                        ( .5 * $v * ($n - 1) )
                    )
                    /
                    # Denominator:
                    (  ($n - 1) * ($n - 2)  )
                );
            }
            $self->{'obrien'}->{$sname} = Statistics::Descriptive::Full->new();
            $self->{'obrien'}->{$sname}->add_data(@r);
            @r = ();
            # Check that each group mean of the O'Briens are equal to the variance of the original data:
            if (sprintf('%.2f', $self->{'obrien'}->{$sname}->mean) != sprintf('%.2f', $v)) {
                croak "Mean for sample $sname does not equal variance";
            }
        }
    # Perform an ANOVA using the O'Brien values as the DV:
    $self->_aov_indep_param_cat($self->{'obrien'});
    $self->{'_eq_var'} = $self->{'_stat'}->{'p_value'} < .05 ? 0 : 1;
    return $self;
}
*obrien_test = \&obrien; # Alias

=head4 levene

 $aov->levene()

I<Alias>: levene_test

Performs B<Levene's> (1960) test for equality of variances within each group: an ANOVA of the absolute deviations, i.e., absolute value of each observation less its group mean.

The statistical attributes now within the class object (see L<anova|anova>) pertain to this test, e.g., $aov->{'f_value'} gives the I<F>-statistic for Levene's Test; and $aov->{'p_value'} gives the I<p>-value associated with the I<F>-statistic for Levene's Test.

=cut

#-----------------------------------------------
sub levene {
#-----------------------------------------------
    my ($self, %args) = @_;
    my %data = ref $args{'data'} eq 'HASH' ? %{$args{'data'}} : ref $self->{'data'} eq 'HASH' ? %{$self->{'data'}} : croak 'No reference to an associative array for performing ANOVA';
    my $tdata = $self->_purge_in_list(\%data);# List-wise clean-up    

    my ($m, $v, $n, @d) = ();       
    $self->{'levene'} = {};
    
    # Traverse each sample of data:
    while(my($sname, $sdata) = each %{$tdata}){
        # For each group, compute the sample mean and the unbiased sample variance:
        $m = $sdata->mean;
        $v = $sdata->variance;
        $n = $sdata->count;
        # For each observation, compute the absolute deviation:
        my @data = $sdata->get_data;
        my $m = $sdata->mean;
        push @d, abs($_ - $m) foreach @data;
        $self->{'levene'}->{$sname} = Statistics::Descriptive::Full->new();
        $self->{'levene'}->{$sname}->add_data(@d);
        @d = ();
    }

    # Perform an ANOVA using the abs. deviations as the DV:
    $self->_aov_indep_param_cat($self->{'levene'});
    $self->{'_eq_var'} = $self->{'_stat'}->{'p_value'} < .05 ? 0 : 1;
    return $self;
}
*levene_test = \&levene; # Alias

=head2 MEASURING EFFECT

Follow-up parametric ANOVAs.

=head3 eta_squared

Returns eta-squared if an ANOVA has been performed; otherwise C<croak>s. Also feeds $aov with the value, named 'eta_sq'. Values range from 0 to 1, 0 indicating no effect, 1 indicating difference between at least two DV means. Generally indicates the proportion of variance in the DV related to an effect.

=cut

#-----------------------------------------------
sub eta_squared {
#-----------------------------------------------
    my ($self) = @_;
    croak 'Need to run ANOVA to obtain requested statistic' if !defined $self->{'_stat'}->{'df_b'} || !defined $self->{'_stat'}->{'f_value'};
    $self->{'_stat'}->{'eta_sq'} = ( $self->{'_stat'}->{'df_b'} * $self->{'_stat'}->{'f_value'} ) / ($self->{'_stat'}->{'df_b'} * $self->{'_stat'}->{'f_value'} + $self->{'_stat'}->{'df_w'} );
    ##$ss_b / ( $ss_b / $ss_w );
    return $self->{'_stat'}->{'eta_sq'};
}

=head3 omega_squared

Returns the effect size statistic omega-squared if an ANOVA has been performed; otherwise C<croak>s. Also feeds $aov with the value, named 'omega_sq'. Generally, size is small where omega_sq = .01, medium if omega_sq = .059, and strong if omega_sq = .138. 

=cut

#-----------------------------------------------
sub omega_squared {
#-----------------------------------------------
    my ($self) = @_;
    my ($eta, $a, $n) = ();
    $eta = $self->eta_squared();
    $a = $eta * ($self->{'_stat'}->{'df_b'} + $self->{'_stat'}->{'df_w'}) - $self->{'_stat'}->{'df_b'};
    $n += $self->{'data'}->{$_}->count foreach keys %{$self->{'data'}};
    $self->{'_stat'}->{'omega_sq'} = $a / ( $a + $n * ( 1 - $eta ) ); 
    return $self->{'_stat'}->{'omega_sq'};
}

=head2 IDENTIFYING RELATIONSHIPS/DIFFERENCES

=head3 compare

 $aov->compare(independent => 1|0, parametric => 1|0, tails => 2|1, flag => 0|1, alpha => .05,
    adjust_p => 0|1, adjust_e => 1|0|2, use_t => 0|1, dump => 0|1, str => 0|1)

The method performs all possible pairwise comparisons, with the Bonferroni approach to control experiment-wise error-rate. The particular tests depend on whether or not you want parametric (default) or nonparametric tests, and if the observations of the factor have been made between groups (default) or by repeated measures.

[TO DO: Note: The following new procedures, as implemented from v0.07 onwards, are only relevant for independent designs; for repeated measures, the comparisons are, at this time, the same as provided in earlier versions, i.e., by multiple paired comparisons I<t>-tests.]

B<Parametric pairwise comparison>. If C<parametric> =E<gt> 1 (default), performs I<F>-tests on each possible pair of observations, with respect to the value of C<independent>. Unless you have just run an equality of variances test, the omnibus equality-of-variances is first tested by the L<O'Brien method|obrien>; otherwise, the result of the last such test is looked up.

=over

=item

I<If the variances are unequal> (I<p> E<lt> .05), the variance of each sample in the pair is used in the error-term of the I<F>-value, and the denominator degrees-of-freedom is adjusted accordingly. 

=item

I<If the variances are equal>, the mean-square error ($aov-E<gt>{'ms_w'}) is used in the denominator. 

=back

You get direct, unadjusted use of the mean-square error, however, with no prior test of equality-of-variances, if you specifically set the parameter C<adjust_e> =E<gt> 0. On the other hand, you can I<force> the procedure to use separate variances, and adjust the error term and degrees-of-freedom accordingly, even if the variances are equal at the .05 level, if C<adjust_e> =E<gt> 2.

B<Non-parametric pairwise comparison>. If C<parametric> =E<gt> 0: derives the I<z>-value and associated I<p>-value for the standardized (a) Wilcoxon (between-groups) sum-of-ranks if C<independent> =E<gt> 1 (B<Dwass> procedure), or (b) the Friedman-type (within-groups) sum-of-ranks if C<independent> =E<gt> 0. Naturally, the number of observations per group should be reasonably large for this I<p>-value to be appropriate; otherwise, look-up $aov->{'s_value'} in the appropriate tables.

Nominality is always assumed. Perhaps some Monte Carlo testing could be useful if the factor is, in fact, at least ordinal.

The I<p>-value is 2-tailed, by default, unless otherwise specified, as above.  If the value of the argument C<adjust_p> equals 1, then the probability values themselves will be adjusted according to the number of comparisons, alpha will remain as given or at .05. The correction is:

=for html <p>&nbsp;&nbsp;&nbsp; <i>p</i>' = 1 &ndash; (1 &ndash; <i>p</i>)<sup>N</sup></p>

where I<p> is the probability returned by the relevant comparison procedure, and I<N> is the number of pairwise comparisons performed.

By default, returns a hashref of hashrefs, the outer hash keyed by the pairs compared (as a comma-separated string), each with a hashref with keys named C<t_value>, C<p_value>, C<df>, C<sig> (= 1 or 0 depending on its being below or greater than/equal to C<alpha>).

Alternatively, if the value of C<str> =E<gt> 1, you just get back a referenced array of strings that describe the results, e.g., G1 - G2: t(39) = 2.378, 2p = 0.0224.

Give a value of 1 to C<dump> to automatically print these strings to STDOUT. (Distinct from earlier versions, there is no default dump to STDOUT of the results.)

The output strings are appended with an asterisk if the logical value of the optional attribute C<flag> equals 1 and the C<p_value> is less than the Bonferroni-adjusted alpha level. This alpha level, relative to the given or default alpha of .05, for the number of paired comparisons, is printed at the end of the list.

An alternative (actually, legacy from earlier version) is to use I<t>-tests, rather than I<F>-tests, and this you get if the argument C<use_t> =E<gt> 1. The module uses Perl's Statistics I<t>-test modules for this purpose, with no accounting for the variance issue.

=cut

sub compare {
    my ($self, %args) = @_;
    foreach (qw/independent parametric nominal adjust_e/) { # 'nominal' an option not yet implemented
        $args{$_} = 1 if !defined $args{$_} ;
    }
    my ($data, $s_value, $a3, $a4, $eq_var, $p_value, $cmp_fn, $p_str, $flag, $flag_str, $pairs, $alpha, @all_pairs, @strings, %res) = ();
    $args{'tails'} ||= 2;

    # Define the data and which routine to use based on values of independent and parametric:
    if ($args{'independent'}) {
        ($data) = $self->_indep_data($args{'data'});
        if ($args{'parametric'}) {
            if (!$args{'use_t'}) {
                $self->obrien() if !defined $self->{'_eq_var'};
                $eq_var = $self->{'_eq_var'};
                $self->anova(independent => 1, parametric => 1, ordinal => 0);
                $cmp_fn = \&_cmp_indep_param_cat;
            }
            else { # legacy offer:
                require Statistics::TTest;
                my $ttest = new Statistics::TTest;
                $cmp_fn = sub {
                    my $data_pairs = shift;
                    $ttest->load_data([$data_pairs->[0]->[1]->get_data()], [$data_pairs->[1]->[1]->get_data()]);
                    $p_value = $ttest->{'t_prob'}; # returns the 2-tailed p_value
                    $p_value /= 2 if $args{'tails'} == 1;    
                    return ($ttest->t_statistic, $p_value, $ttest->df);
                };
            }
        }
        else {
            #$cmp_fn = !$args{'ordinal'} ? \&_cmp_indep_dfree_cat : \&_cmp_indep_dfree_ord;
            $cmp_fn = \&_cmp_indep_dfree_cat;
        }
     }
     else {
        ($data) = $self->_rmdep_data($args{'data'});

        if ($args{'parametric'}) {
                $cmp_fn = \&_cmp_rmdep_param_cat;
        }
        else {
            carp 'Non-parametric multi-comparison procedure for dependent/repeated measures is not implemented';
        }
    }

    @all_pairs = combinations([keys(%{$data})], 2); # Algorithm::Combinatorics fn
    $alpha = $args{'alpha'} || .05;
    $alpha /= scalar(@all_pairs) if !$args{'adjust_p'}; # divide by number of comparisons

    # Compare each pair:
    foreach $pairs (@all_pairs) {
        $pairs = [sort {$a cmp $b} @{$pairs}];
        ($s_value, $p_value, $a3, $a4) = $cmp_fn->([ [$pairs->[0], $data->{$pairs->[0]}], [$pairs->[1], $data->{$pairs->[1]}] ], %args, eq_var => $eq_var, ms_w => $self->{'_stat'}->{'ms_w'});
         $p_value = _pcorrect($p_value, scalar(@all_pairs)) if $args{'adjust_p'};
         $a3 = _precisioned($args{'precision_s'}, $a3); # degrees-of-freedom
         $s_value = _precisioned($args{'precision_s'}, $s_value);
         $p_value = _precisioned($args{'precision_p'}, $p_value);
         $p_str = $args{'tails'} == 1 ? '1p' : '2p';
         $flag = $p_value < $alpha ? 1 : 0;
         $flag_str = $args{'flag'} ? $flag ? ' *' : '' : '';
         if ($args{'parametric'}) {
             $res{"$pairs->[0],$pairs->[1]"} = { t_value => $s_value, p_value => $p_value, df => $a3, flag => $flag};
             push @strings,  "($pairs->[0] - $pairs->[1]), t($a3) = $s_value, $p_str = $p_value" . $flag_str;
         }
         else {
             $res{"$pairs->[0],$pairs->[1]"} = { z_value => $s_value, p_value => $p_value, s_value => $a3, flag => $flag};
             push @strings, "($pairs->[0] - $pairs->[1]), Z(W) = $s_value, $p_str = $p_value" . $flag_str;
        }
     } # end loop
        
     if ($args{'dump'}) {
         print "$_\n" foreach @strings;
         print "Alpha = $alpha\n";
     }

    return $args{'str'} ? \@strings : \%res;
}

sub _cmp_indep_param_cat {
    my ($data_pairs, %args) = @_;
    my @nn = ($data_pairs->[0]->[1]->count(), $data_pairs->[1]->[1]->count());
    my @uu = ($data_pairs->[0]->[1]->mean(), $data_pairs->[1]->[1]->mean());
    my @cc = (1, -1); # contrast coefficients for pairwise comparison of means (ui - uj = 0)
    my ($f_value, $f_prob, $df_w, $lu, $ss_cc) = ();
    $lu += $cc[$_] * $uu[$_] foreach (0, 1); # estimate of linear combo of means; s/= zero if no diff.
    $ss_cc += $cc[$_]**2 / $nn[$_] foreach (0, 1); # weighted sum of squared contrast coefficients
    
    if ( ($args{'adjust_e'} && !$args{'eq_var'}) || $args{'adjust_e'} == 2) { # permits unequal variances; sep. error terms/contrast
        my $ss_u = $lu**2 / $ss_cc;
        my @vv = ($data_pairs->[0]->[1]->variance(), $data_pairs->[1]->[1]->variance());
        my ($denom, $den_a, $df_n, $df_d) = ();
        $den_a += ($cc[$_]**2 / $nn[$_])*$vv[$_] foreach (0, 1);
        $denom = $den_a / $ss_cc;
        $f_value = $ss_u / $denom; # Maxwell & Delaney Eq. 5.10 (p. 180)
        $df_n += ($cc[$_]**2 * $vv[$_] / $nn[$_]) foreach (0, 1);
        $df_d += ($cc[$_]**2 * $vv[$_] / $nn[$_])**2 / ($nn[$_] - 1) foreach (0, 1);
        $df_w = $df_n**2 / $df_d;
     }
     else { # use the pooled error-term (MSe):
         $f_value = $lu**2 / ( $ss_cc *  $args{'ms_w'} ); # Maxwell & Delaney Eq 4.37 (p. 176); assumes equal variances
         ##$f_value = ( $nn[0] * $nn[1] * ( $uu[0] - $uu[1])**2 ) / ( ( $nn[0] + $nn[1] ) * $self->{'_stat'}->{'ms_w'} );
         $df_w = ($nn[0] + $nn[1] - 2);
     }
     $f_prob = fdtrc(1, $df_w, $f_value); # s/be compared to alpha/nContrasts
     return ($f_value, $f_prob, 1, $df_w);
}

sub _cmp_rmdep_param_cat { # not completed - only use Dependant t-test for now:
    my ($data_pairs, %args) = @_;
    require Statistics::DependantTTest;
    my $ttest = new Statistics::DependantTTest; 
    $ttest->load_data($data_pairs->[0]->[0], $data_pairs->[0]->[1]->get_data());
    $ttest->load_data($data_pairs->[1]->[0], $data_pairs->[1]->[1]->get_data());
    my ($s_value, $df) = $ttest->perform_t_test($data_pairs->[0]->[0], $data_pairs->[1]->[0]);
    my $p_value = stdtr($df, -1 * abs($s_value)); # Math::Cephes fn # gets the left 1-tailed p_value
    $p_value *= 2 unless $args{'tails'} == 1;
    return ($s_value, $p_value, $df);
    #my ($i, $pairs, %diffs, @full_diffs, @null_diffs, @F_x_sum, @R_x_sum) = ();
    #my $plim = scalar(keys %{$data});##scalar(@{$all_pairs});
	#my $n;
}


sub _cmp_indep_dfree_cat { # Dwass-Steel procedure
    my ($data_pairs, %args) = @_;
    my ($n1, $n2, $nm, $sum, $exp, $var, $term_tie, $z_value, $p_value, $pstr) = ($data_pairs->[0]->[1]->count, $data_pairs->[1]->[1]->count); # Ns/group
    $nm = $data_pairs->[1]->[0]; # arbitrarily use second group as reference data
    my ($ranks_href, $ties_var, $gn, $xtied) = _ranks_between(_aref2href($data_pairs)); # get joint ranks
    $sum = sum(@{$ranks_href->{$nm}}); # calc. sum-of-ranks for the (arbitrarily) second member of the pair
    $exp = ( $n2 * ($gn + 1) ) / 2; # expected value
    $term_tie += ($_ - 1 ) * $_ * ($_ + 1) foreach @{$xtied};
    $var = ( ($n1 * $n2) / 24 ) * ( ($gn + 1 - ( $term_tie / ( ($gn) * ($gn - 1) ) ) ) );# variance
    $z_value =  ( $sum - $exp ) / sqrt($var) ; # standardized sum-of-ranks
    $p_value = ( 1 - ndtr( abs($z_value) ) ); # Math::Cephes fn
    $p_value *= 2 unless $args{'tails'} == 1;    
    return ($z_value, $p_value, $sum);
}

=head3 cluster

 $aov->cluster(parametric => 1|0, tails => 1|2, flag => 1|0, alpha => .05, adjust_p => 1|0, dump => 1|0, str => 1|0)

Identifies the first-order clusters arising by all possible binary splits of the groups/levels, selecting that which gives the lowest error sum-of-squares (or greatest between sum-of-squares).

B<Parametric binary clustering>. If C<parametric> =E<gt> 1 (default), the method is as described by Scott & Knott (1974). 

B<Non-parametric binary clustering>.  If C<parametric> =E<gt> 0, Worsley's (1977) nonparametric implementation of the Scott-Knott method is offered: the observed values are replaced by their rank, and the maximum Kruskal-Wallis statistic given by the binary partitions is selected; ties are corrected as described above for the Kruskal-Wallis test, rather than the random assignment method described by Worsley.

In both the parametric and nonparametric cases, the method returns a hashref like so:

 k_value : the test statistic pertaining to the selected partition
 p_value : probability associated with the test statistic
 clusters : referenced array of referenced arrays, each giving the names of groups/levels included in each cluster

(A future version might seek second-order partitions - i.e., doing this all over again for each of the first-order partitions.)

=cut

sub cluster {
    my ($self, %args) = @_;
    my ($ttest, $data) = ();
    my ($s_value, $p, $df, $p_value, $cmp_fn, $p_str, $flag, $pairs, $alpha, %means, @strings, %res) = ();
    $args{'tails'} ||= 2;
    $data = $self->_purge_in_list($self->{'data'});# List-wise clean-up
    $means{$_} = $data->{$_}->mean() foreach keys(%{$data});
    my @parts = partitions([keys %means], 2); # Algorithm::Combinatorics function # not reliably/substantially faster if undeclared
    my ($statistic, $pair) = $args{'parametric'} ? _cluster_param($data, @parts) : _cluster_dfree($data, @parts);
    foreach (@{$pair}) {
        print @{$_}, "\n";
    }
    $self->{'_stat'}->{'k_value'} = $statistic;
}

sub _cluster_param {
    my ($data, @parts) = @_;
    my ($d1, $d2, $p, $ss_w, $df_w, $desc, $kn, $levels, $gm, $gvar, $min_pr, @clusters, %ss_ws) = ();
    foreach $p(@parts) {  # Get within-groups sum-of-squares for each possible partition by 2:
        $d1 = Statistics::Descriptive::Full->new();
        $d2 = Statistics::Descriptive::Full->new();
        $d1->add_data($data->{$_}->get_data) foreach @{$p->[0]};
        $d2->add_data($data->{$_}->get_data) foreach @{$p->[1]};
        ($ss_w, $df_w, $desc, $gm, $gvar) = _indep_error({ 1 => $d1, 2 => $d2});# error SS & dfs (et al.):
        #$ss_b += $_->[0] * ($_->[1] - $grand_mean)**2 foreach @{$desc}; # treatment SS; $desc is aref of arefs of counts [0] & means [1] per gp
        $ss_ws{$ss_w} = [ [$p->[0], $d1], [$p->[1], $d2] ];# assume no ties; as aref not href ($pdata) for ease of keeping the names 
    }
    # Determine minimum within-groups sum-of-squares for which partition, & calc stat:
    $min_pr = $ss_ws{ min(keys %ss_ws) };
    $kn += $_->[1]->count * ($_->[1]->mean - $gm)**2 foreach @{$min_pr}; # calc max between SS
    $levels = $kn / $gvar;
    push @clusters, [sort {$a cmp $b} @{$_->[0]}] foreach @{$min_pr};
    return ($levels, \@clusters);
}

sub _cluster_dfree {
    my ($data, @parts) = @_;
    my ($d1, $d2, $p, $df_b, $n, $h, @clusters, %kws) = ();
    foreach $p(@parts) {
        $d1 = Statistics::Descriptive::Full->new();
        $d2 = Statistics::Descriptive::Full->new();
        $d1->add_data($data->{$_}->get_data) foreach @{$p->[0]};
        $d2->add_data($data->{$_}->get_data) foreach @{$p->[1]};
        ($h, $df_b, $n) = _kw_stat({ 1 => $d1, 2 => $d2});
        $kws{$h} = [$p, $df_b, $n];
    }
    $h = max(keys %kws);
    push @clusters, [sort {$a cmp $b} @{$_}] foreach @{$kws{$h}->[0]};
    return ($h, \@clusters);
}

# Could be 
# alpha => [.05, .005] (give as proportions 0-1 not percentages 1-100).
# tails => 1|2
# exp_tail = -1|1; (negative - less than zero - or positive)

sub _is_significant {
    my ($pvalue, $alpha, $exp_tail) = @_;
    $exp_tail ||= 2;
    if (ref $alpha) {
        # Assume aref:
        if ($pvalue > $alpha->[0]) {
            if ($exp_tail) {
                
            }
        }
    }
}

=head3 confidence

 $itv_str = $aov->(independent => 1|0, alpha => .05, name => 'aname', limits => 0) # get interval for single group as string
 $lim_aref = $aov->(independent => 1|0, alpha => .05, name => 'aname', limits => 1) # get upper & lower limits for single group as aref
 $itv_href = $aov->(independent => 1|0, alpha => .05, name => ['aname', 'bname'], limits => 0) # get interval for 2 groups as hashref keyed by group names
 $lim_href = $aov->(independent => 1|0, alpha => .05, name => ['aname','bname'], limits => 1) # get upper & lower limits for 2 groups as hashref of group-named arefs
 $itv_href = $aov->(independent => 1|0, alpha => .05, name => undef, limits => 0) # get intervals for all groups as hashref keyed by group names
 $lim_href = $aov->(independent => 1|0, alpha => .05, name => undef, limits => 1) # upper & lower limits for all groups as hashref 

Computes confidence intervals using (by default) the pooled estimate of variability over groups/levels, rather than the standard error within each group/level, as described by Masson and Loftus (2003). For a between groups design, the confidence interval (as usual) indicates that, at a certain level of probability, the true population mean is likely to be within the interval returned. For a within-subjects design, as any effect of the variability between subjects is eliminated, the confidence interval (alternatively) indicates the reliability of the how the sample means are distributed as an estimate of the how the population means are distributed.

In either case, there is an assumption that the variances within each condition are the same between the conditions (homogeneity of variances assumption).

Actual algorithm depends on whether the measures are obtained from seperate groups (independent => 1) or by repeated measures (independent => 0) (i.e., whether between-groups or within-groups design). Default is between-groups.

The option C<use_mse> can be set to equal 0 so that the (typical) standard error of the mean is used in place of the mean-square error. This is one option to use when the variances are unequal.

The option C<conditions> can, optionally, include a referenced array naming the particular conditions that should be included when calculating I<MSe>. By default, this is all the conditions, using I<MSe> from the omnibus ANOVA. This is one option to handle the case of unequal variances between conditions.

Default 

=cut

sub confidence {
#-----------------------------------------------
    my ($self, %args) = @_;
    #croak 'The <confidence> procedure is not yet implemented';
    croak 'Need to run ANOVA to obtain requested statistic' if !defined $self->{'_stat'}->{'df_w'} || !defined $self->{'_stat'}->{'ms_w'};

    my $data = $self->_purge_in_list($self->{'data'});# List-wise clean-up    

    # Init key params:
    my $indep = defined $args{'independent'} ? $args{'independent'} : 1;
    my $alpha = _init_alpha($args{'alpha'}); # default = .05
    my $tcrit = abs(stdtri($self->{'_stat'}->{'df_w'}, $alpha/2));
    my $limits = delete $args{'limits'} or 0;
    my $use_mse = defined $args{'use_mse'} ? $args{'use_mse'} : 1;
    my @names = defined $args{'name'} ? ref $args{'name'} ? @{$args{'name'}} : ($args{'name'}) : keys(%{$data});
    my @conditions = ref $args{'conditions'} ? @{$args{'conditions'}} : @names;
print "tcrit = $tcrit\n";
    my ($erv, $itv, %confints) = ();
    
    foreach (@names) {
        if ($use_mse) {
            my $mse;
            #if (ref $args{'conditions'}) {
            #    print "cond n = ",scalar(@{$args{'conditions'}}), "\n";
            #    print "cond g = ", scalar(keys(%{$data})), "\n";
            # and scalar(@{$args{'conditions'}}) < scalar(keys(%{$data})) ) {
                # get the mse for the particular comparison:
            #    my $aov_i = Statistics::ANOVA->new();
            #    $aov_i->add($_ => $data->{$_}->get_data) foreach @{$args{'conditions'}};
            #    $mse = $aov_i->{'_stat'}->{'ms_w'};
            #    print "mse by contrast = $mse\n";
            #}
            #else {
                $mse = $self->{'_stat'}->{'ms_w'};
            #    print "mse omni = $mse\n";
            #}
            $erv = sqrt( $mse / $data->{$_}->count );
        }
        else {
            $erv = $data->{$_}->standard_deviation / sqrt($data->{$_}->count - 1);
        }
        $itv = $erv * $tcrit;
        if ($limits) {
            $confints{$_} = [$data->{$_}->mean - $itv, $data->{$_}->mean + $itv];
        }
        else {
            $confints{$_} = $itv;
        }
    }
        
    return scalar(keys(%confints)) > 1 ? \%confints : $confints{$names[0]};
    
    my $levels = scalar(keys(%{$data}));
    my $n_comp = $levels * ($levels - 1) / 2;
    print "ncomp = $n_comp\n";
    
    $alpha /= $n_comp; # divide by number of comparisons
    my $f_val = fdtri(1, $self->{'_stat'}->{'df_w'}, (.05 / $n_comp ) ); # Math::Cephes fn
    my $w_var = sqrt($f_val);
    my $sum = 1;
    ##$sum +=  1 / $data->{$_}->count foreach keys(%{$data});
    my $cf = $w_var * sqrt( $self->{'_stat'}->{'ms_w'} * $sum);
    print "cf = $cf, f= $f_val, w = $w_var\n";
}

=head2 ACCESSING RESULTS

=head3 string

 $str = $aov->string(mse => 1, eta_squared => 1, omega_squared => 1, precision_p => integer, precision_s => integer)

Returns a statement of result, in the form of C<F(df_b, df_w) = f_value, p = p_value>; or, for Friedman test C<chi^2(df_b) = chi_value, p = p_value> (to the value of I<precision_p>, if any); and so on for other test statistics. Optionally also get MSe, eta_squared and omega_squared values appended to the string, where relevant. These and the test statistic are "sprintf"'d to the I<precision_s> specified (or, by default, not at all).

=cut

#-----------------------------------------------
sub string {
#-----------------------------------------------
    my ($self, %args) = @_;
    my $str;
    my $p_value = $args{'precision_p'} ? sprintf('%.' . $args{'precision_p'} . 'f', $self->{'_stat'}->{'p_value'}) : $self->{'_stat'}->{'p_value'};
    my $precision_s = $args{'precision_s'} || 0;
    if (defined $self->{'_stat'}->{'f_value'} && !$self->{'_dfree'}) {
        $str .= "F($self->{'_stat'}->{'df_b'}, $self->{'_stat'}->{'df_w'}) = ";
        $str .= _precisioned($precision_s, $self->{'_stat'}->{'f_value'} );
        $str .= ", p = $p_value,";
        $str .= ' MSe = ' . _precisioned($precision_s, $self->{'_stat'}->{'ms_w'}) . ',' if $args{'mse'};
        $str .= ' eta^2 = ' . _precisioned($precision_s, $self->eta_squared()). ',' if $args{'eta_squared'};
        $str .= ' omega^2 = ' . _precisioned($precision_s, $self->omega_squared()). ',' if $args{'omega_squared'};
        chop($str);
    }
    elsif (defined $self->{'_stat'}->{'h_value'}) { # Kruskal-Wallis statistic
        $str .= "H($self->{'_stat'}->{'df_b'}) = ";
        $str .= _precisioned($precision_s, $self->{'_stat'}->{'h_value'});
        $str .= ", p = $p_value";
    }
    elsif (defined $self->{'_stat'}->{'j_value'}) { # Jonckheere-Terpstra statistic
        $str .= "J = ";
        $str .= _precisioned($precision_s, $self->{'_stat'}->{'j_value'});
        $str .= ", p = $p_value";
    }
    elsif (defined $self->{'_stat'}->{'l_value'}) { # Page statistic
        $str .= "L = ";
        $str .= _precisioned($precision_s, $self->{'_stat'}->{'l_value'});
        $str .= ", p = $p_value";
    }
    elsif (defined $self->{'_stat'}->{'chi_value'}) { # Friedman statistic
        $str .= "chi^2($self->{'_stat'}->{'df_b'}) = ";
        $str .= _precisioned($precision_s, $self->{'_stat'}->{'chi_value'});
        $str .= ", p = $p_value";
    }
    else {
        $| = 1;
        croak 'Need to run omnibus test (anova) to obtain results string'
    } 
    return $str;
}

=head3 table

 $table = $aov->table(precision_p => integer, precision_s => integer);

Returns a table listing the degrees of freedom, sums of squares, and mean squares for the tested "factor" and "error" (between/within groups), and the I<F>- and I<p>-values. The test statistics are "sprintf"'d to the I<precision_s> specified (or, by default, not at all); the p value's precision can be specified by I<precision_p>. 

Up to this version, if calculating any of these values was not essential to calculation of the test statistic, the value will simply appear as a blank in the table. If the omnibus test last made was non-parametric, and no I<F>-value was calculated, then the table returned is entirely an empty string.

Formatting with right-justification where appropriate is left for user-joy.

=cut

#-----------------------------------------------
sub table {
#-----------------------------------------------
    my ($self, %args) = @_;
    my $tbl = '';
    my $precision_p = $args{'precision_p'} || 0;
    my $precision_s = $args{'precision_s'} || 0;
    # F-table:
    if (defined $self->{'_stat'}->{'f_value'} && !$self->{'_dfree'}) {
        $tbl .= "\t$_" foreach ('Df', 'Sum Sq', 'Mean Sq', 'F value', 'Pr(>F)');
        $tbl .=  "\n";
        $tbl .=  "$_\t" foreach ('Factor', $self->{'_stat'}->{'df_b'}); 
        $tbl .= _precisioned($precision_s, $_) . "\t" foreach ($self->{'_stat'}->{'ss_b'}, $self->{'_stat'}->{'ms_b'}, $self->{'_stat'}->{'f_value'});
        $tbl .= _precisioned($precision_p, $self->{'_stat'}->{'p_value'});
        $tbl .=  "\n";
        $tbl .=  "$_\t" foreach ('Error', $self->{'_stat'}->{'df_w'}); 
        $tbl .= _precisioned($precision_s, $_) . "\t" foreach ($self->{'_stat'}->{'ss_w'}, $self->{'_stat'}->{'ms_w'});
        $tbl .=  "\n";
    }
    return $tbl;
}

=head3 dump

 $aov->dump(title => 'ANOVA test', precision_p => integer, precision_s => integer, mse => 1, eta_squared => 1, omega_squared => 1, verbose => 1)

Prints the string returned by L<string|string>, or, if specified with the attribute I<table> => 1, the table returned by L<table|table>; and the string as well if I<string> => 1. A newline - "\n" - is appended at the end of the print of the string. Above this string or table, a title can also be printed, by giving a value to the optional C<title> attribute.

If I<verbose> => 1, then any curiosities arising in the calculations are noted at the end of other dumps. At the moment, this is only the number of observations that might have been purged were they identified as undefined or not-a-number upon loading/adding.

=cut

#-----------------------------------------------
sub dump {
#-----------------------------------------------
    my ($self, %args) = @_;
    print "$args{'title'}\n" if $args{'title'};
    if ($args{'table'}) {
        print $self->table(%args);
        print $self->string(%args), "\n" if $args{'string'};
    }
    else {
        print $self->string(%args), "\n";
    }
    print "Observations purged as undefined or not-a-number: ". $self->{'purged'} . "\n" if $self->{'purged'} && $args{'verbose'};
}

#-----------
# Private methods

sub _prepare { # get data and check group count for anova - called by _indep_data and _rmdep_data
    my ($self, $data) = @_;
    my $dataref = ref $data eq 'HASH' ? $data : ref $self->{'data'} eq 'HASH' ? $self->{'data'} : ref $data eq 'ARRAY' ? _aref2href($data) : ref $self->{'data'} eq 'ARRAY' ? _aref2href($self->{'data'}) : croak 'No reference to a hash of data for performing ANOVA';
    croak 'Not enough groups, if any, in the data for performing ANOVA' if scalar(keys(%{$dataref})) <= 1;
    return ($dataref);
}

sub _indep_data {
    my $self = shift;
    return $self->_purge_in_list( $self->_prepare(@_) ); # list-wise clean-up
}

sub _indep_ss_ord {
    my $data = shift;
    my @names = keys(%{$data});
    croak "Check names for groups: All need to be numerical for trend analysis" if grep { !looks_like_number($_)} @names;
    # (Benchmark: about 50% faster to "my" the vars apiece rather than as a list to "my" beforehand:)
    my $mean_t = mean( @names );
    my $sum_sample_contrasts = sum( map { $data->{$_}->mean * ($_ - $mean_t) } @names );
    my $sum_squared_coeffs = sum( map { ($_ - $mean_t)**2 / $data->{$_}->count }  @names ); # unweighted
    return $sum_sample_contrasts**2 / $sum_squared_coeffs;
}

sub _indep_error { # essentially to get error sum-squares and error degrees-of-freedom, but other params piggy-back upon it:
    my ($data, $mean, $samples, $ss_w, $df_w, @dat, @desc) = (shift);
    my $all = Statistics::Descriptive::Sparse->new();

    foreach (keys %{$data}) {
        @dat = $data->{$_}->get_data;
        $all->add_data(@dat);
        $mean = $data->{$_}->mean;
        $samples = $data->{$_}->count;
        $ss_w += ($_ - $mean)**2 foreach @dat;# accumulate within-groups SS, and df
        $df_w += ($samples - 1);
        push @desc, [$samples, $mean]; # for calculating ss_b
    }
    croak 'No within-groups data for performing ANOVA' if !$ss_w || !$df_w;
    return ($ss_w, $df_w, \@desc, $all->mean(), $all->variance());
}

sub _rmdep_data {
    my ($self) = @_;
    ##my $samples = _testdata_counts($self->_prepare(@_)); # prelim check of equals Ns required before purging
    my $data = $self->_purge_across_cases( $self->_prepare(@_) ); # assume equal Ns, test later
    return ($data, _testdata_counts($data), scalar(keys(%{$data}))); # data-hash, N-per-level, N-levels
}

sub _rmdep_ss_uni { # UNIVARIATE method: error and treatment sums-of-squares for rm anovas (et al.):
    my ($data, $samples, $levels, $ss_w, $df_w, $ss_b, $df_b, $grand_mean, $i, @i_means, %i_data, %j_means) = (shift, shift, shift);

    # Mean over each index:
    for ($i = 0; $i < $samples; $i++) {
        $i_data{$i} = Statistics::Descriptive::Full->new();
        $i_data{$i}->add_data( ($data->{$_}->get_data)[$i] ) foreach keys %{$data};
        $i_means[$i] = $i_data{$i}->mean();
    }

    $grand_mean += $_ foreach @i_means;
    $grand_mean /= scalar(@i_means);

    foreach (keys %{$data}) {
        $j_means{$_} = $data->{$_}->mean();
        $ss_b += ($j_means{$_} - $grand_mean)**2;
    }
    $ss_b *= $samples;

    foreach (keys %{$data}) {
        #my @sdata = $data->{$_}->get_data;
        for ($i = 0; $i < $samples; $i++) {
            $ss_w += (($data->{$_}->get_data)[$i] - $i_means[$i] - $j_means{$_} + $grand_mean)**2;
        }
    }

    $df_b = $levels - 1;
    $df_w = $df_b * ( $samples - 1);

    return ($ss_w, $df_w, $ss_b, $df_b, $grand_mean);
}    

# -------------------------------------------------------------------------------
# TO DO:
sub _rmdep_ss_multi { # MULTIVARIATE method
    croak "The multivariate method for performing a repeated-measures ANOVA is not yet implemented";
}

#-------------------------------------------------------------------------------------------
# RANKING methods for NONPARAMETRIC tests
#-------------------------------------------------------------------------------------------

sub _ranks_between { # used by Kruskal-Wallis, Jonckheere-Terpstra and Worsley-cluster tests
    my $data = _group(shift);
    my ($i, $rank, $ties_var, $lim, $n_groups, @xtied, @sorted, @groups, %ranks) = (0, 1, 0);
    @sorted = sort {$a <=> $b} keys %{$data};
    $lim = scalar @sorted;

    for ($i = 0; $i < $lim; $i++) {
        @groups = @{$data->{$sorted[$i]}};
        $n_groups = scalar @groups;
        if ($n_groups > 1) {
            $ties_var += ($n_groups**3 - $n_groups);
            push @{$ranks{$_}}, mean($rank .. $rank + $n_groups - 1) foreach @groups;
            $rank += $n_groups;
        }
        else {
            push @{$ranks{$groups[0]}}, $rank++;
        }
        push @xtied, $n_groups;
    }
    $rank--;
    return (\%ranks, $ties_var, $rank, \@xtied); # rank hash-of-arefs, tie-correction, N, ari of tied group Ns 

    sub _group {
        my $data = shift;
        my ($name, @dat, %grouped) = ();
        foreach $name(keys %{$data}) {
            @dat = $data->{$name}->get_data;
            foreach (@dat) {
                if (exists($grouped{$_})) {
        	        push @{$grouped{$_}}, $name;
        	    }
        	    else {
        		    $grouped{$_} = [$name];
            	}
     	    }
        }
        return \%grouped;
    }
}

sub _ranks_within { # used by Friedman and Page tests
    my ($data, $samples) = @_;
    my ($old, $cur, $i, $col, $rval, $ties, $av_rank, %ranks, %row_values) = (0, 0);
    my (%xtied) = ();
    for ($i = 0; $i < $samples; $i++) { # - set the averaged ranks, going down each index:
        # - list the values at this index in each data-array (being Statistics-Descriptive objects):
        # - a value might occur in more than one group at this index, so store an array of the groups:
        push @{ $row_values{ ($data->{$_}->get_data)[$i] } }, $_ foreach keys %{$data};
        # This loop adapted from Gene Boggs' "rank" function in Statistics-RankCorrelation:
        for $rval (sort { $a <=> $b } keys %row_values) {
            # Get the number of ties:
            $ties = scalar(@{ $row_values{$rval} });
            $cur += $ties;
            if ($ties > 1) {
                $av_rank = $old + ($ties + 1) / 2; # Average the tied data
                $ranks{$_} += $av_rank for @{ $row_values{$rval} };
                push @{$xtied{$i}}, $ties;
            }
            else {
                $ranks{ $row_values{$rval}[0] } += $cur; # Add the single rank to the list of ranks
                push @{$xtied{$i}}, $ties;
            }
            $old = $cur;
        }
        ($old, $cur, %row_values) = (0, 0);
    }
    return (\%ranks, \%xtied);
}

sub _kw_stat {
    my ($data, $do_corr) = @_;
    my ($ranks_href, $ties_var, $gn) = _ranks_between($data);
    my ($num, $h) = (0);
    $num += scalar(@{$ranks_href->{$_}}) * ( mean(@{$ranks_href->{$_}}) - ( ($gn + 1) / 2 ) )**2 foreach keys %{$ranks_href};        
    $h = 12 / ($gn *($gn + 1)) * $num;
    $h /= (1 - ($ties_var / ($gn**3 - $gn))) unless defined $do_corr and !$do_corr;
    return ($h, (scalar( keys %{$ranks_href}) - 1), $gn); # H, df, and grand N
}

sub _aref2href {
    my $aref = shift;
    my %hash = ();
    $aref = [$aref] if ! ref $aref->[0];
    foreach (@{$aref}) {
        if (ref $_->[1]) {
            $hash{$_->[0]} = $_->[1];
        }
        else {
            my $name = shift(@{$_});
            $hash{$name} = [@{$_}];
        }
    }
    return \%hash;
}

sub _testdata_counts { # ensure equal observations for dependent groups analyses
    my $data = shift;
    my $samples;
    foreach (keys %{$data}) {
        if (defined $samples) {
            my $tmp_count = $data->{$_}->count;
            if ($tmp_count != $samples) {
                croak 'Number of observations per group need to be equal for repeated measures ANOVA';
            }
            else {
                $samples = $tmp_count;
            }
        }
        else {
            $samples = $data->{$_}->count;
        }
    }
    croak 'One observation per group can not be tested in ANOVA' if $samples == 1;
    return $samples;
}

sub _mark_purge { # List keyed by sample-names of their indices where invalid values lie
    my ($self, $name, $dat, $n, $i) = @_;
    $n = scalar @{$dat};
    for ($i = 0; $i < $n; $i++) {
        $self->{'_purge'}->{$name}->{$i}++ if !looks_like_number($dat->[$i]);
    }
}

sub _purge_in_list {
    my ($self, $data) = @_;
    my ($i, $name, @clean, @data_single, %tdata) = ();
    $self->{'purged'} = 0;
    foreach $name(keys %{$data}) {
        @data_single = $data->{$name}->get_data();
        for ($i = 0; $i < scalar @data_single; $i++) {
            if ($self->{'_purge'}->{$name}->{$i}) {
                $self->{'purged'}++;
            }
            else {
                push @clean, $data_single[$i];
            }
        }
        croak "Empty data sent to ANOVA following purge of invalid value(s) in list < $name >" if !scalar(@clean);
        $tdata{$name} = Statistics::Descriptive::Full->new();
        $tdata{$name}->add_data(@clean);
        @clean = ();
    }
    return \%tdata;
}

sub _purge_across_cases {
    my ($self, $data) = @_;
    my ($i, $cull, $key, $val, $name, %invalid_ids, %clean, %tdata) = ();

    # List of all indices in all lists with invalid values; and copy of each group of data:
    foreach $name(keys %{$data}) {
        $clean{$name} = [$data->{$name}->get_data()];
        while(($key, $val) = each %{$self->{'_purge'}->{$name}}) {
            $invalid_ids{$key} += $val;
        }
    }
    $self->{'purged'} = scalar(keys(%invalid_ids)) || 0;

    # Purge by index (from highest to lowest):
    my @invalid_ids = reverse( sort {$a <=> $b} keys(%invalid_ids) );
    foreach $cull(@invalid_ids) {
       foreach $name(keys %clean) {
           splice(@{$clean{$name}}, $cull, 1);
       }
    }

    foreach (keys %clean) {
        croak "Empty data sent to ANOVA following purge of invalid value(s) in list < $_ >" if !scalar(@{$clean{$_}});
        $tdata{$_} = Statistics::Descriptive::Full->new();
        $tdata{$_}->add_data($clean{$_});
    }

    return \%tdata;
}

sub _pcorrect { # (1 - ( 1 - p)^N )
    return  1 - ( 1 - $_[0] )**$_[1] ;
}

sub _precisioned {
    return $_[0] ? sprintf('%.' . $_[0] . 'f', $_[1]) : (defined $_[1] ? $_[1] : ''); # don't lose any zero
}

sub _init_alpha {
    my $val = shift;
    if (defined $val) {
        if ($val > 0 && $val < 1) {
            return $val;
        }
        else {
            croak "Alpha value should be between 0 and 1, not '$val'."
        }
    }
    else {
        return .05;
    }
}

1;

__END__

=head1 REFERENCES

Hollander, M., & Wolfe, D. A. (1999). I<Nonparametric statistical methods>. New York, NY, US: Wiley.

Levene, H. (1960). Robust tests for equality of variances. In I. Olkins (Ed.), I<Contributions to probability and statistics>. Stanford, CA, US: Stanford University Press.

Masson, M. E. J., & Loftus, G. R. (2003). Using confidence intervals for graphically based data interpretation. I<Canadian Journal of Experimental Psychology>, I<57>, 203-220.

Maxwell, S. E., & Delaney, H. D. (1990). I<Designing experiments and analyzing data: A model comparison perspective.> Belmont, CA, US: Wadsworth.

O'Brien, R. G. (1981). A simple test for variance effects in experimental designs. I<Psychological Bulletin>, I<89>, 570-574.

Page, E. B. (1963). Ordered hypotheses for multiple treatments: A significance test for linear ranks. I<Journal of the American Statistical Association>, I<58>, 216-230.

Scott, A. J., & Knott, M. (1974). A cluster analysis method for grouping means in the analysis of variance. I<Biometrics>, I<30>, 507-512.

Siegal, S. (1956). I<Nonparametric statistics for the behavioral sciences>. New York, NY, US: McGraw-Hill

Worsley, K. J. (1977). A non-parametric extension of a cluster analysis method by Scott and Knott. I<Biometrics>, I<33>, 532-535.

=head1 SEE ALSO

L<Math::Cephes|lib::Math::Cephes> Probabilities for all tests are computed using this module's functions, rather than the "in-house" L<Statistics::Distributions|lib::Statistics::Distributions> module, as the former appears to be more accurate for larger values of F.

L<Statistics::Descriptive|Statistics::Descriptive> Fundamental calculations of means and variances are left up to this standard; any limitations/idiosyncrasies therein are naturally passed onto the present one; although the present one purges missing and non-numerical values, unlike Statistics::Descriptive, which gives them the value of "0" (and so returns erroneous descriptives).

L<Statistics::FisherPitman|lib::Statistics::FisherPitman> For an alternative to independent groups ANOVA when the variances are unequal.

L<Statistics::KruskalWallis|lib::Statistics::KruskalWallis> Offers Newman-Keuls for pairwise comparison by ranks. Also offers non-parametric independent groups ANOVA, but note it does not handle ties in rank occurring between two or more observations, nor correct for them; an erroneous H-value is calculated if ties exist in your data. Also does not handle missing/invalid values. Present module adapts its _grouped method.

L<Statistics::Sequences|lib::Statistics::Sequences> Offers methods for comparison of independent groups by split over central tendency.

L<Statistics::Table::F|Statistics::Table::F> Simply returns an F value. Does not handle missing values, treating them as zero and thus returning an erroneous F-value in these cases.

=head1 BUGS/LIMITATIONS/TO DO

Optimisation welcomed.

No adjustment is offered for violations of sphericity in repeated measures ANOVA, but the multivariate approach is a work in progress.

To do: Repeated measures planned/pairwise comparisons requires bringing to the level as for independent measures.

To do: Confidence rather than/in addition to comparison and cluster breakdowns.

To do: Parametric ANOVAs for ordinal level of measurement in repeated measures designs. 

To do: Verbose tracing to check-out what's actually happening internally, what defaults are being used, etc.

=head1 REVISION HISTORY

See CHANGES in installation distribution.

=head1 AUTHOR/LICENSE

=over 4

=item Copyright (c) 2006-2012 Roderick Garton

rgarton AT cpan DOT org

This program is free software. It may be used, redistributed and/or modified under the same terms as Perl-5.6.1 (or later) (see L<http://www.perl.com/perl/misc/Artistic.html>).

=item Disclaimer

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=back

=cut
