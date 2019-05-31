#!/usr/bin/env perl

use strict;
use warnings;

open IN, "<$ARGV[0]" or die $!;
open OUT, ">$ARGV[1]" or die $!;

my @lines = <IN>;
my $line = join "", @lines;
$line =~ s/\s//g;

my $newStr = "";
my $preEnd = 0;
while ($line =~ /([\w\.]+):[\d\.]+/g) {
    my $id = $1;
    my $len = pos($line) - length($&) - $preEnd;
    my $frag = substr($line, $preEnd, $len);
    $newStr = $newStr.$frag.$id;
    $preEnd = pos($line);
}

if ($preEnd < length($line)) {
    my $frag = substr($line, $preEnd, length($line) - $preEnd);
    $newStr = $newStr.$frag;
}

$newStr =~ s/\):/\)/g;
print OUT "$newStr\n";
    

