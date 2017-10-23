use strict;

@ARGV == 2 or die "Usage: $0 fbin ifile >ofile\n";
my $fbin = shift;

sub vadd(\@@)
{  my $r = shift;
   my $i;
   $$r[$i++] += $_ foreach @_;
}

sub vscale(\@$)
{  my ($r, $v) = @_;
   $_ *= $v foreach @$r;
}

my @d;
my @s;
my $n;
my $c;

while(<>)
{  chomp;
   s/^\s+//;
   @d = split /\s+/;
   if ($c == 0)
   {  # next bin
      if (@s)
      {  vscale @s, 1/$n;
         print join("\t", @s), "\n";
         undef @s;
      }
      $c = $n = int($d[0] * $fbin +1);
   }
   # add bin
   #print "L: $n $c @d\n";
   vadd @s, @d;
   --$c;
}
