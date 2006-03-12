#!/usr/bin/perl
use strict;

@ARGV >= 2 or die "usage: $0 source1 source2 [source3...] > destination\n";

my @data;
my $lastline;

foreach my $file (@ARGV)
{  warn "Reading $file...\n";
   open F, $file or die "Failed to open source file $file.\n";
   # check number of lines
   if ($lastline)
   {  while (<F>)
      {  chomp;
         split ' ';
         my $r = $data[$.-1];
         $r->[0] == $_[0] or die "X-axis is different at line $. in $file\n"
                               . "now: $_[0], previously: $r->[0]\n";
         $r->[1] += $_[1];
         $r->[2] = $_[2] if $_[2] > $r->[2];
      }
      $lastline == $. or die "The number of lines in the source files are inconsistent.\n"
                           . "$file: $., previously: $lastline.\n";
   } else
   {  while (<F>)
      {  chomp;
         push @data, [split ' '];
      }
      $lastline = $.;
   }
   close F;
}

foreach (@data)
{  $_->[1] /= @ARGV;
   print join("\t", @$_), "\n";
}

