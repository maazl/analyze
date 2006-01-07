#!/usr/bin/perl
use strict;

my @data;

sub vadd2(\@@)
{  my $aref = shift;
   my $x = shift;
   die "X axis incompatible found $x, expected $$aref[0]\n" if $$aref[0] && abs($x - $$aref[0])/$$aref[0] > 1.E-5;
   $_ += shift foreach (@$aref[1..$#$aref]); # existing columns
   push @$aref, @_;              # new columns
}

sub dataadd(@)
{  if (@data)
   {  die "Number of source lines incompatible.\n" if @data != @_;
      vadd2 @$_, @{shift()} foreach (@data); # existing lines
   } else
   {  # first time
      @data = @_;
   }
}

# read data
foreach (@ARGV)
{  open F, $_ or die " Failed to open $_\n";
   dataadd map { chomp; [split] } <F>;
   close F;
}

# write data
foreach (@data)
{  print join(' ', map { sprintf "%15.9g", $_ } @$_), "\n";
}

