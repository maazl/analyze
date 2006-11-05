use strict;

@ARGV or die "Usage: $0 size ifile >ofile\n";

my $avg = shift @ARGV;

sub vadd(\@@)
{  my $rarr = shift;
   foreach (@$rarr)
   {  $_ += shift;
   }
}

sub vsub(\@@)
{  my $rarr = shift;
   foreach (@$rarr)
   {  $_ -= shift;
   }
}

my @lastl;
my @current;

my $hist = 0;
while (++$hist < $avg && ($_ = <>))
{  chomp;
   my @d = split;
   $#current = $#d;
   vadd @current, @d;
   push @lastl, \@d;
}

while (<>)
{  chomp;
   my @d = split;
   vadd @current, @d;
   print join(" ", map sprintf("%12g", $_/$avg), @current), "\n";
   vsub @current, @{shift @lastl};
   push @lastl, \@d;
}


