#!/usr/bin/perl

use strict;
use warnings;

my @seqs;
my @ids;

my $cutoff = $ARGV[0];
open IN, "<$ARGV[1]" or die $!;
open OUT, ">$ARGV[2]" or die $!;

$/ = ">";
while (my $line = <IN>) {
    chomp $line;
    my ($header, $dna) = split(/\n/, $line, 2);
    next unless ( $header && $dna);
    if ($header =~ /^(\S+)/){
        my $id = $1;
        $dna =~ s/\s//g;
        my @seq = split("", $dna);
        push @ids, $id;
        push @seqs, \@seq;
   }
}
$/ = "\n";


my $row = scalar @ids;
my $col = scalar @{$seqs[0]};

#print "$row\t$col\n";

for (my $j = 0; $j < $col; $j++) {
    my $gap_num = 0;
    my $z = 0;
    for (my $i = 0; $i < $row; $i++) {
        my $char = $seqs[$i]->[$j];
        $gap_num++ if ($char eq "-");
        $z++ if ($char eq "Z" or $char eq "B");
    }
    if ($gap_num / $row > $cutoff or $z > 0) {
        for (my $i = 0; $i < $row; $i++) {
            $seqs[$i]->[$j] = "*";
        }
    }
}

for (my $i = 0; $i < $row; $i++) {
    my $id = $ids[$i];
    my $seq = join("", @{$seqs[$i]});
    $seq =~ s/\*//g;
    print OUT ">$id\n$seq\n";
}
