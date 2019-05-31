#!/usr/bin/env perl

use strict;
use warnings;
use proscan;
use fasta;

my $domain_id = $ARGV[0];
my $spe = $ARGV[1];
my $pep_file = $ARGV[2];
my $scan_file = $ARGV[3];
my $out_dir = $ARGV[4];


my %id2blocks;

my $fa = fasta->new($pep_file);
my $id2seq = $fa->{"id2seq"};
my $scan = proscan->new($scan_file);
$scan->load_gene2blocks();
my $gene2blocks = $scan->{"gene2blocks"};

foreach my $gene_id (keys %$gene2blocks) {
    my $blocks = $gene2blocks->{$gene_id};
    foreach my $block (@$blocks) {
        my $motif_id = $block->{"motif_id"};
        if ($motif_id eq $domain_id) {
            if (exists $id2blocks{$gene_id}) {
                push @{$id2blocks{$gene_id}}, $block;
            } else {
                $id2blocks{$gene_id} = [$block];
            }
        }
    }
}

open ID, ">$out_dir/${spe}.id" or die $!;
open PEP, ">$out_dir/${spe}.pep" or die $!;
open DOM, ">$out_dir/${spe}.domain" or die $!;

foreach my $gene_id (sort keys %id2blocks) {
    my $pep = $id2seq->{$gene_id};
    print ID "$gene_id\n";
    print PEP ">$gene_id\n$pep\n";
    
    my $blocks = $id2blocks{$gene_id};
    my @sort_blocks = sort {$a->{'start'} <=> $b->{'start'}} @$blocks;
    # only print N-terminal domain 
    my $block = $sort_blocks[0];
    my $start = $block->{'start'};
    my $end = $block->{'end'};
    my $len = $end - $start + 1;
    my $frag = substr($pep, $start - 1, $len);
    print DOM ">$gene_id\n$frag\n";
}
  
close ID;
close PEP;
close DOM;




