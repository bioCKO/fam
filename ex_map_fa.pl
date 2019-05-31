#!/usr/bin/env perl

use strict;
use warnings;
use util;
use fasta;

my $id2name = util::load_map_id($ARGV[0]);
my $fa = fasta->new($ARGV[1]);
open OUT, ">$ARGV[2]" or die $!;

my $id2Seq = $fa->{"id2seq"};

foreach my $id (keys %$id2name) {
    my $name = $id2name->{$id};
    if (exists $id2Seq->{$id}) {
        print OUT ">$name\n";
        print OUT "$id2Seq->{$id}\n";
    } else {
        print "$id: No seq found!\n";
    }
}
