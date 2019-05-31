#!/usr/bin/env perl


use strict;
use warnings;

open FILE, "<$ARGV[0]" || die $!;
open OUT, ">$ARGV[1]" or die $!;

my $tree;
my %id2name;
while (my $line = <FILE>) {
    chomp $line;
    if ($line =~ /^\ttranslate/) {
        while ($line = <FILE>) {
            chomp $line;
            last if $line =~ /;/;
            if ($line =~ /^\t+(\d+)\t([^,]+)/) {
                $id2name{$1} = $2;
            }
        }
    }
    
    if ($line =~ /^\s+tree[^\(]+(.+)$/) {
        $tree = $1;
    }
}
close FILE;

$tree =~ s/\s//g;

while ($tree =~ /\[[^\]]+\]/g) {
    if ($tree =~ /\)\[&prob=([^,]+),[^\]]+\]/) {
        my $d = $1;
        $d = $d * 1.0;
        $tree =~ s/\)\[[^\]]+\]/\)$d/;
    } else {
        $tree =~ s/\[[^\]]+\]//;
    }
    my $end = pos $tree;  
    pos $tree = $end;    
}

pos $tree = 0;

while ($tree =~ /[\(,]\d+:/g) {
    $tree =~ s/\((\d+):/\($id2name{$1}:/;
    $tree =~ s/,(\d+):/,$id2name{$1}:/;
    my $end = pos $tree;  
    pos $tree = $end;   
}

#$line =~ s/\):/\)/g;
print OUT "$tree\n";
close OUT;

