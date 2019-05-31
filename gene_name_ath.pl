#!/usr/bin/env perl

use strict;
use warnings;

open AT2SYM, "<$ARGV[0]" or die $!;
open GENE, "<$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]" or die $!;

#T1G01010   ANAC001 NAC domain containing protein 1
#AT1G01010  NAC001  NAC domain containing protein 1

my %at2sym;
while (my $line = <AT2SYM>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @array = split /\t/, $line;
    if ($at2sym{$array[0]}) {
        my $syms = $at2sym{$array[0]};
        push @$syms, $array[1];
    } else {
        $at2sym{$array[0]} = [$array[1]];
    }
}
close AT2SYM;


while (my $line = <GENE>) {
    chomp $line;
    next if ($line =~ /^\s*$/);
    my $at = $line;
    my $sym = "";
    if ($at2sym{$at}) {
        my $array = $at2sym{$at};
        $sym = join("\t", @$array); 
        print OUT "$at\t$sym\n";
    } else {
        print OUT "$at\t$sym\n";
    }
    
}



