use strict;

sub prim($)
{  use integer;
   my $val = shift;
   my @result;

   sub checkfactor(\$$) # value, divisor
   {  my $rval = shift;
      my $div = shift;
      my $res = 0;
      while ($$rval % $div == 0)
      {  $$rval /= $div;
         ++$res;
      }
      return [$div, $res];
   }

   push @result, checkfactor $val, 2;

   my $div = 3;

   while ($val > 1)
   {  push @result, checkfactor $val, $div;
      $div += 2;
   }

   return @result;
}

sub showprim(@)
{  map { print "$$_[0]^$$_[1] "; } @_;
   print "\n";
}

sub kgV(@)
{  use integer;
   my @prims = map [prim $_], @_;
   my $res = 1;

   sub getnextmax(@)
   {  my $div;
      my $max = 0;
      map
      {  my $e = shift @$_;
         if ($e && $$e[1] > $max)
         {  $div = $$e[0];
            $max = $$e[1];
         }
      } @_;
      return ($div, $max);
   }

   while (1)
   {  my ($div, $p) = getnextmax @prims;
      last unless $div;
      $res *= $div ** $p;
   }

   return $res;
}

sub ggT(@)
{  use integer;
   my @prims = map [prim $_], @_;
   my $res = 1;

   sub getnextmin(@)
   {  my $div;
      my $min = 999999999;
      map
      {  my $e = shift @$_;
         if ($e)
         {  if ($$e[1] < $min)
            {  $div = $$e[0];
               $min = $$e[1];
            }
         } else
         {  $min = 0;
         }
      } @_;
      return ($div, $min);
   }

   while (1)
   {  my ($div, $p) = getnextmin @prims;
      last unless $div;
      $res *= $div ** $p;
   }

   return $res;
}


#map { print "$_ -> "; showprim prim $_; } @ARGV;

#print kgV(@ARGV), "\n";
#print ggT(@ARGV), "\n";


my $sfreq = 48000;

# ***** read config

my @freq;

while (<>)
{  chomp;
   next unless $_;
   next if /^#/;
   my @d = split;
   if ($d[0] =~ /:/)
   {  #range
      my @r = $d[0] =~ /^(\S+):(?:(\S+):)?(\S+)$/ or die "illegal range";
      $r[1] = 1 if $r[1] eq '';
      for (my $i = $r[0]; $i <= $r[2]; $i+=$r[1])
      {  $d[0] = $i;
         push @freq, [@d];
      }
   } else
   {  push @freq, \@d;
}  }


# ***** prepare
my $fqlow = ggT map $$_[0], @freq;
#my $fqlow = 1;
#my $fqhigh = kgV map $$_[0], @freq;

#use integer;
#my $totalsamp = ($sfreq + $fqlow/2) / $fqlow;
#no integer;
my $totalsamp = int $sfreq / $fqlow;

print STDERR "low frequency limit: $fqlow\n";
#print STDERR "high frequency limit: $fqhigh\n";
print STDERR "total samples: $totalsamp\n";
print STDERR "\n";


# ***** generate
my @samp;

my $PI = 3.14159265358979323846;

sub generatefq(@)
{  my ($fq, $vol, $ph) = @_;
   $ph = rand*2*$PI if !defined $ph;
   print STDERR "generating $fq\n";
   $fq /= $fqlow;
   for (my $i = 0; $i < $totalsamp; ++$i)
   {  @samp[$i] += $vol * sin 2*$PI * ($i * $fq / $totalsamp + $ph);
   }
}

generatefq @$_ foreach @freq;


# ***** normalize
my $sum = 0;
my $sum2 = 0;
my $max = 0;

foreach (@samp)
{  $sum += $_;
   $sum2 += $_ * $_;
   $max = abs $_ if abs $_ > $max;
}

$sum /= $totalsamp;
$sum2 /= $totalsamp;

print STDERR "average: $sum\n";
print STDERR "RMS: $sum2\n";
print STDERR "Maximum: $max\n";
print STDERR "\n";

sub scale(\@$)
{  my $rarr = shift;
   my $fact = shift;
   $_ *= $fact foreach @$rarr;
}

scale @samp, 32767/$max;


# ***** info

my $fscale = $sfreq / $totalsamp / $fqlow;

print STDERR "freq.\tintens.\tphase\n";
print STDERR $$_[0] * $fscale, "\t", $$_[1]*32767/$max, "\t$$_[2]\n" foreach @freq;


# ***** write

#print "$_\n" foreach @samp;

binmode STDOUT;
map { print pack "S", $_; } @samp;


