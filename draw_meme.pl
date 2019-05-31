#!/usr/bin/env perl

use strict;
use warnings;
use util;
use meme;
use fasta;
use SVG;

my $gene_ids = util::load_ids_array($ARGV[0]);
my $fa = fasta->new($ARGV[1]);
my $meme = meme->new($ARGV[2]);

my $seq2len = $fa->get_id2seq_len($gene_ids);
my $gene2blocks = $meme->get_gene2blocks();
my $motif_ids = $meme->get_motif_ids();

my $thick = 10;
my $y_inter = $thick + $thick/2;
my $x_left = 200;
my $str_x_off = 10;
my $y_top = 50;
my $axis_y = $y_top - 20;
my $bar_width = 600;

my $max_len = 0;
foreach my $key (keys %$seq2len) {
    my $len = $seq2len->{$key};
    $max_len = $len if $len > $max_len;
}

my $scale = $max_len / $bar_width;
my $im_width = $bar_width + 2 * $x_left;
my $im_heigth = $y_inter * scalar(@$gene_ids) +  $y_inter * scalar(@$motif_ids) + 2 * $y_top;


my %motif2color;
{
    my %rs;
    my %gs;
    my %bs;
    my ($r, $g, $b);
    foreach my $id (@$motif_ids) {    
        do {
            $r = int(rand 255);
            $g = int(rand 255);
            $b = int(rand 255);
        } while (exists $rs{$r} or exists $gs{$g} or exists $bs{$b});   
        $rs{$r} = 1;
        $gs{$g} = 1;
        $bs{$b} = 1;
        $motif2color{$id} = "rgb($r, $g, $b)";
    }
}


my $svg = SVG->new("width" => $im_width, "height" => $im_heigth);
$svg->line("x1" => $x_left, "y1" => $axis_y, "x2" => $x_left + $max_len / $scale, "y2" => $axis_y, "stroke" => "black");

for (my $i = 0; $i < $max_len; $i = $i + 10) {
    $svg->line("x1" => $x_left + $i / $scale, "y1" => $axis_y - 5, "x2" => $x_left + $i / $scale, 
                    "y2" => $axis_y, "stroke" => "black", "stroke-width" => 0.5);
    if ($i%100 == 0) {
        $svg->line("x1" => $x_left + $i / $scale, "y1" => $axis_y - 10, "x2" => $x_left + $i / $scale, "y2" => $axis_y, "stroke" => "black", "stroke-width" => 0.5);
        $svg->text("x" => $x_left + $i / $scale, "y" => $axis_y - 10, "font" => "Arial", "font-size" => 8)->cdata("${i}");
    }
}


my $i = 0;
foreach my $gene_id (@$gene_ids) {
    my $y = $y_top + $i * $y_inter;
    my $y1 = $y - $thick / 2;
    my $y2 = $y + $thick / 2;
    $svg->text("x" => $str_x_off, "y" => $y2, "font_size" => $thick)->cdata($gene_id);
    $svg->line("x1" => $x_left, "y1" => $y, "x2" => $x_left + $seq2len->{$gene_id} / $scale,
                "y2" => $y, "stroke" => "black", "stroke-width" => 2);
    if (exists $gene2blocks->{$gene_id}) {
        my $blocks = $gene2blocks->{$gene_id};
        foreach my $block (@$blocks) {
            my $motif_id = $block->{"motif_id"};
            my $x = $x_left + $block->{"start"} / $scale;
            my $width = $block->get_len() / $scale;
            $svg->rectangle("x" => $x, "y" => $y1, "width" => $width, "height" => $thick,
                            "fill" => $motif2color{$motif_id}, "stroke" => "black",
                            "stroke-width" => 0.5);     
            #$svg->text("x" => $x + $width / 2 - 5, "y" => $y2 - 2, "-cdata" => $motif_id,
             #           "font" => "Arial", "font-size" => 8, "fill" => "white");
        }
    }
    $i++;
}


$y_top = $y_top + $i * $y_inter + 20;
my $j = 0;

my $id2motif = $meme->{"id2motif"};
    
my @sort_ids = sort {$a <=> $b} keys %{$id2motif};
foreach my $motif_id (sort {$a <=> $b} keys %{$id2motif}) {
    my $motif = $id2motif->{$motif_id};
    my $y = $y_top + $j * $y_inter;
    my $y1 = $y - $thick / 2;
    my $y2 = $y + $thick / 2;
    $svg->text("x" => $x_left - 20, "y" => $y2, "text" => $motif_id,
                "font-size" => $thick, "-cdata" => "Motif $motif_id");
    $svg->rectangle("x" => $x_left + 20, "y" => $y1, "width" => 10, "height" => $thick,
                    "fill" => $motif2color{$motif_id}, "stroke" => "black",
                    "stroke-width" => 0.5);
    $svg->text("x" => $x_left + 50, "y" => $y2, "-cdata" => $motif->{"regex"},
               "font-size" => $thick);
    $j++;
}

print $svg -> xmlify;


