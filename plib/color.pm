package color;

use strict;
use warnings;

sub rand_color {
	my $ids = shift;
	my %id2color;
	my %rs;
	my %gs;
	my %bs;
	my ($r, $g, $b);

	foreach my $id (@$ids) {	
		do {
			$r = int(rand 255);
			$g = int(rand 255);
			$b = int(rand 255);
		} while (exists $rs{$r} or exists $gs{$g} or exists $bs{$b});	
		$rs{$r} = 1;
		$gs{$g} = 1;
		$bs{$b} = 1;
		$id2color{$id} = "rgb($r, $g, $b)";
	}
	return \%id2color;
}

1;
