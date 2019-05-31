#!/usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;

my $ioIn = Bio::AlignIO->new(-format => "fasta", -file => $ARGV[0]);
my $ioOut = Bio::AlignIO->new(-format => "phylip", -file=> ">$ARGV[1]" , 
                                -idlength => 20, -interleaved => 1);

my $alnIn = $ioIn->next_aln();
$ioOut->write_aln($alnIn);
