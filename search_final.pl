#!/usr/bin/env perl

use strict;
use warnings;

my $hmm = $ARGV[0]
my $dir1 = $ARGV[1];
my $dir2 = $ARGV[2];
my $dir3 = $ARGV[3];

if (-d $dir3) {
    die "output dirtectory exist!";
} else {
    mkdir $dir3, 0777 or die $!;
}

opendir DIR, $dir2 or die $!;
while (my $file = readdir DIR) {
    next if $file =~ /^\./;
    if ($file =~ /^(\w+)\.scan$/) {
        my $spe = $1;
        my $genome_pep_file = "$dir1/${spe}.pep";
        my $scan_out_file = "$dir2/$file";
        system("parse_interpro.pl $hmm $spe $genome_pep_file $scan_out_file $dir3");
    }
}

