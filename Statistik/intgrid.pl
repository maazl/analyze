#!/usr/bin/perl
use strict;

sub vadd(\@@)
{  my $aref = shift;
   $_ += shift foreach (@$aref); # existing columns
   push @$aref, @_;              # new columns
}

sub vdiv(\@@)
{  my $aref = shift;
   $_ /= shift foreach (@$aref); # existing columns
   push @$aref, map 0, @_;       # new columns
}

sub vscale(\@$)
{  my $aref = shift;
   my $scale = shift;
   $_ *= $scale foreach (@$aref); # existing columns
}


my $steps = $ARGV[0] ? shift @ARGV : 10;
my $fbin = $ARGV[0] ? shift @ARGV : 0;
my $bin;

# read data
my @data = <>;
$_ = [split] foreach @data;
printf STDERR "%u lines read.\n", scalar @data;

# calculate integral
my @sum;
foreach (@data)
{  @$_ == 3 or die "Wrong number of columns\n";
   vadd @sum, @$_;
}
shift @sum;

# write data
my @step = @sum;
vscale @step, 1/$steps;
my @nxt = @step;
my @cur;
my $binc;
my @avg;
foreach my $d (@data)
{  vadd @cur, @$d[1..$#$d];
   # bin
   vadd @avg, @$d;
   next unless ++$binc >= int $$d[0]*$fbin;
   vscale @avg, 1/$binc;
   $binc = 0;
   # data
   my $f = shift @avg;
   vdiv @avg, @sum;

   my @cn = map
   {  #print STDERR "$f $cur[$_], $nxt[$_], $step[$_]\n";
      if ($cur[$_] > $nxt[$_])
      {  $nxt[$_] += $step[$_];
         $avg[$_]
      } else
      {  0
      }
   } (0..1);
   printf "%12g  %15.9g %15.9g %15.9g %15.9g\n", $f, @avg, @cn;

   undef @avg;
}
#print STDERR "@cur\n@sum\n@nxt";
