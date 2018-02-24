#!/usr/bin/perl
use strict;
use 5.010;

@ARGV >= 3 or die "usage:\n"
	. "  datacheck.pl test_file column(s) check_file [tolerance abs, rel]\n"
	. "  datacheck.pl test_file column(s) { formula } [tolerance abs, rel]\n";

my $testfile = shift;
my @cols = split ",", shift;
my $check = shift;
my $tolerance_abs = shift; # might be 0
my $tolerance_rel = 1 + shift; # might be 1

sub openfile($)
{	open my $F, "<", $_[0] or die "Failed to open $_[0]: $!\n";
	return $F;
}
sub fetchrow($)
{	my $l;
	do
	{	defined($l = readline $_[0]) or return;
		#print "L: $l\n";
	} while $l =~ /^#/;
	chomp $l;
	$l =~ s/^s+//; $l =~ s/\s+$//;
	return split /\s*\t\s*/, $l;
}

sub is_nan($)
{	return $_[0] != $_[0];
}
sub cmpval($$$)
{	return if $_[0] == $_[1] || !defined $_[1];
	if (is_nan $_[0])
	{	return if is_nan $_[1];
	} elsif (!is_nan $_[1] && ($tolerance_abs || $tolerance_rel))
	{	return if abs($_[0] - $_[1]) <= $tolerance_abs + abs($_[1]) * $tolerance_rel;
	}
	die "Value does not compare equal: found $_[0], expected $_[1]. $testfile line $. column $_[2]\n";
}


my $FT = openfile $testfile;
if ($check =~ /^\{/)
{	#print "$check\n";
	eval "sub reffunc(\@) $check";
	die "Invalid formula: $@\n" if $@;
} else
{	my $FC = openfile $check;
	sub reffunc(@)
	{	return fetchrow $FC or die "Additional lines in test file $testfile.\n";
	}
}
my $i = 0;
while (my @data = fetchrow $FT)
{	unshift @data, $i++;
	my @ref = reffunc(@data);
	#print "C: [@data[@cols]] [@ref]\n";
	for (0..$#cols)
	{	cmpval $data[$cols[$_]], $ref[$_], $cols[$_];
	}
}
