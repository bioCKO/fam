#!/usr/bin/env perl

use strict;
use warnings;
use Bio::TreeIO;

my $in = Bio::TreeIO->new(-format => "newick", -file => $ARGV[0]);
open OUT, ">$ARGV[1]" or die $!;

while (my $tree = $in->next_tree) {
	my @nodes = $tree->get_nodes();
	foreach my $nd (@nodes) {
		print OUT $nd->id, "\n" if $nd->is_Leaf ;
	}
}

