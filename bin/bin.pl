use strict;

my $avg = shift @ARGV;

sub vadd(\@@) {
	my $rarr = shift;
	foreach (@$rarr) {
		$_ += shift;
	}
}

my @current;

open F, shift @ARGV;

my $line = 0;
while (<F>) {
	chomp;
	my @d = split;
	$#current = $#d;
	vadd @current, @d;
	if ( ++$line == $avg ) {
		print join( " ", map sprintf( "%12g", $_ / $avg ), @current ), "\n";
		undef @current;
		$line = 0;
	}
}

