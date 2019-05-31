#!/usr/bin/perl

use strict;
use warnings;
use util;

my $id2name = util::load_map_id($ARGV[0]);
open IN, "<$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]" or die $!;

while (my $line = <IN>) {
    chomp $line;
    #next if $line =~ /^\s*$/;
    my $new_str = "";
    my $pre_end = 0;
    while ($line =~ /([\w\.-]+)/g) {
        my $id = $1;
        my $len = pos($line) - length($id) - $pre_end;
        my $frag = substr($line, $pre_end, $len);
    
        if (exists $id2name->{$id}) {
            $new_str = $new_str.$frag.$id2name->{$id};
        } else {
            $new_str = $new_str.$frag.$id;
        }
        $pre_end = pos($line);
    }
    
    if ($pre_end < length($line)) {
        my $frag = substr($line, $pre_end, length($line) - $pre_end);
        $new_str = $new_str.$frag;
    }
    print OUT "$new_str\n";
}


