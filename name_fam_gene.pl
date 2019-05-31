#!/usr/bin/env perl

use strict;
use warnings;

my $spe = $ARGV[0];
open IN, "<$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]" or die $!;
open OUT1, ">$ARGV[3]" or die $!;

my @a = ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r',  's', 't',  'u', 'v', 'w', 'x',  'y', 'z');

my %target2querys; # {target_name => [[query_name, evalue], ...] , ...}

#MDP0000119457  AtCPK8  28.333  60  43  0   357 416 436 495 2.38e-04    34.3
while (my $line = <IN>) {
    chomp $line;
    my @items = split /\t/, $line;
    next if @items != 12;
    my $queryName = $items[0];
    my $targetName = $items[1];
    my $evalue = $items[10];
    if (exists $target2querys{$targetName}) {
        push @{$target2querys{$targetName}}, [$queryName, $evalue];
    } else {
        $target2querys{$targetName} = [[$queryName, $evalue]];
    }
}
        
foreach my $target (sort keys %target2querys) {
    my $name = substr($target, 1);
    my $querys = $target2querys{$target};
    my @sorted = sort {$a->[1] <=> $b->[1]} @$querys;
    if (@sorted == 1) {
   		my $q = $sorted[0];
        print OUT "$q->[0]\t$spe$name\n";
        print OUT1 "$q->[0]\t$target\t$q->[1]\n";
    } else {
        for (my $i = 0; $i <= $#sorted; $i++) {
            my $q = $sorted[$i];
            print OUT "$q->[0]\t$spe${name}$a[$i]\n";
            print OUT1 "$q->[0]\t${target}\t$q->[1]\n";
        }
    }
}

