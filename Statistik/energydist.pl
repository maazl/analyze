#!/usr/bin/perl
use strict;

@ARGV >= 2 or die "usage: energydist <statistics_file> <filter_order> ft [ft [ft ...]]\n";

my $col = 1;

my $file = shift @ARGV;
my $order = 2*shift @ARGV;

my @fs;
my @fq;
@fq = sort { $a<=>$b }
map
{  unless (/^[0-9\.]*$/)
   {  # load data file
      print STDERR "L: $_\n";
      open D, $_ or die "Failed to load data file $_\n";
      my @d;
      map
      {  chomp;
         local @_;
         split;
         push @d, [$_[0], $_[5]*cos($_[6]*.017453292)];
      } <D>;
      close D;
      unshift @d, [0, $d[0][1]];
      push @d, [100000, $d[$#d][1]];
      push @fs, [@d];
      ()
   } else
   {  $_
   }
} @ARGV;
my @fqmin = @fq;
unshift @fqmin, 0;
my @fqmax = @fq;
push @fqmax, 1E5;


# vector operations
sub transpose(@)
{  # max dimension
   my $dim;
   foreach (@_)
   {  $dim = @$_ if @$_ > $dim;
   }
   my @r;
   for (my $i = 0; $i < $dim; ++$i)
   {  push @r, [map $$_[$i], @_];
   }
   return @r;
}

sub vadd(\@@)
{  my $aref = shift;
   $_ += shift foreach (@$aref); # existing columns
   push @$aref, @_;              # new columns
}

sub sum(@)
{  my $r;
   $r += shift while @_;
   return $r;
}

# filter
my @filter = transpose \@fqmin, \@fqmax, \@fs; # [band][high, low, weight]

# dofilter
# highpass: dofilter f0/f
# lowpass:  dofilter f/f0
sub dofilter($)
{  return $order
    ? 1/(1 + $_[0]**$order)
    : ((1 <=> $_[0]) + 1) /2; # ideal filter
}

sub doweight($$)
{  return 1 unless $_[1]; # no data => no weight
   my ($f,$d) = @_;
   # find position
   my $low = 0;
   my $high = @$d;
   #print "F: $low,$high\t$f\n";
   while ($low < $high)
   {  my $m = int(($low+$high)/2);
      #print "S: $m\t$$d[$m][0]\t$low,$high\n";
      my $r = $$d[$m][0] <=> $f;
      if ($r > 0)
      {  $high = $m;
      } elsif ($r < 0)
      {  $low = $m+1;
      } else
      {  $low = $m;
         $high = $m+1;
         last;
      }
   }
   --$high;
   # interpolate
   #my $r = ($f-$$d[$high][0]) / ($$d[$high+1][0]-$$d[$high][0]) * ($$d[$high+1][1]-$$d[$high][1]) + $$d[$high][1];
   #print "I: $low,$high\t$r\n";
   return ($f-$$d[$high][0]) / ($$d[$high+1][0]-$$d[$high][0]) * ($$d[$high+1][1]-$$d[$high][1]) + $$d[$high][1];
}

sub applyfilter($@)
{  return dofilter($_[0]/$_[2]) * dofilter($_[1]/($_[0]+1E-30)) / doweight($_[0], $_[3]);
}

# apply function & to each element of \@ prepended by the constant arguments @ and return the results.
# this is a generalized form of map.
sub operate21(&\@@)
{  my $rfn = shift;
   my $rarr = shift;
   return map { &$rfn(@_, @$_) } @$rarr;
}

# data
my @sums;

#sub storesum($)
#{  push @sums, $_[0];
#   $tsum += $_[0];
#}

# read data
{  open F, $file or die "Failed to open $file.\n";
   while (<F>)
   {  chomp;
      my @l = split;
      #$tsum += $l[$col];
      vadd @sums, map $l[$col] * $_, operate21 \&applyfilter, @filter, $l[0];
   }
   close F;
}

my $tsum = sum @sums;

# write result
{  print " fmin   fmax   %\n";
   printf "%5.0f %6.0f %6.2f\n", shift @fqmin, shift @fqmax, 100 * shift(@sums) / $tsum while (@sums);
}

