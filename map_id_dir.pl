#!/usr/bin/perl

use strict;
use warnings;

my $map_id_file = $ARGV[0];
my $in_dir = $ARGV[1];
my $out_dir = $ARGV[2];

if (! -d $out_dir) {
	mkdir ($out_dir, 0777) or die $!;
} else {
	die "$out_dir is exist!\n";
}

opendir DIR, $in_dir or die $!;
while (my $file = readdir DIR) {
	next if $file =~ /^\./;
	system("map_id.pl $map_id_file $in_dir/$file $out_dir/$file");
}


