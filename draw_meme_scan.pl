#!/usr/bin/env perl

use strict;
use warnings;
use meme;
use fasta;
use proscan;
use util;
use SVG;

my $id_file = $ARGV[0];
my $seq_file = $ARGV[1];
my $meme_file = $ARGV[2];
my $scan_file = $ARGV[3];


my $gene_ids = util::load_ids_array($id_file);
my $ids= util::load_ids_hash($id_file);
my $fa = fasta->new($seq_file);
my $meme = meme->new($meme_file);
my $seq2len = $fa->get_id2seq_len($gene_ids);
my $gene2motifs = $meme->get_gene2blocks();
my $motif_ids = $meme->get_motif_ids();
my $proscan = proscan->new($scan_file);
$proscan->load_gene2tools();
my $gene2tools = $proscan->{"gene2tools"};

my $scan_num = 0;

foreach my $gene_id (keys %$gene2tools) {
	if (exists $ids->{$gene_id}) {
		my $tool2blocks = $gene2tools->{$gene_id};
		my $tool_num = scalar keys %$tool2blocks;		
		$scan_num += $tool_num;
	}
}


my $x_str_left = 10;
my $thick = 5;
my $y_inter = $thick + $thick / 2;
my $x_left = 200;
my $y_top = 30;
my $bar_width = 500;

my $max_len = 0;
foreach my $key (keys %$seq2len) {
    my $len = $seq2len->{$key};
    $max_len = $len if $len > $max_len;   
}
my $scale = $max_len / $bar_width;
my $id_num = scalar(@$gene_ids) + $scan_num;

my $im_width = $bar_width + 2 * $x_left;
my $im_height = $y_inter * $id_num + $y_inter * scalar(@$motif_ids) + 2 * $y_top;
my $svg = SVG->new(width=>$im_width, height=>$im_height);

$svg->line("x1", $x_left, "y1", $y_top - 10, "x2", $x_left + $max_len / $scale, "y2", $y_top - 10, "stroke", "black", "stroke-width", 1);
for (my $i = 0; $i < $max_len; $i = $i + 10) {
    $svg->line("x1", $x_left + $i / $scale, "y1", $y_top - 15, "x2", $x_left + $i / $scale, "y2", $y_top - 10, "stroke", "black", "stroke-width", 0.5);
    if ($i % 100 == 0) {
        $svg->line("x1", $x_left + $i / $scale, "y1", $y_top - 20, "x2", $x_left + $i / $scale, "y2", $y_top - 10, "stroke", "black", "stroke-width", 0.5);
        $svg->text("x", $x_left + $i/ $scale, "y", $y_top - 20, "font", "Times", "font-size", 8)->cdata("${i}");
    }
}


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

my $i = 0;
foreach my $gene_id (@$gene_ids) {
    next if !exists $gene2motifs->{$gene_id};
    my $motifs = $gene2motifs->{$gene_id};
    my $y = int($y_top + $i * $y_inter);
    my $y1 = $y + $thick / 2;
    my $y2 = $y + $thick;
    $svg->text("x" => $x_str_left, "y" => $y2, "font-size" => $thick, "fill" => "blue", "-cdata" => $gene_id);
    $svg->line("x1" => $x_left, "y1" => $y1, 
                "x2" => $x_left + $seq2len->{$gene_id} / $scale, "y2" => $y1, 
                "stroke" => "black");
    foreach my $motif (@$motifs) {
        my $motif_id = $motif->{"motif_id"};
        my $x = $x_left + int($motif->{"start"} / $scale);      
        my $width = $motif->get_len() / $scale;
        $svg->rectangle("x" => $x, "y" => $y, "width" => $width, "height" => $thick, 
                        "fill" => $motif2color{$motif_id}, "stroke" => "black",
                        "stroke-width" => 0.5);      
        $svg->text("x" => $x, "y" => $y2, "-cdata" => $motif_id, 
                    "font-family" => "Times", "font-size" => $thick, 
                    "fill" => "white");
    }
    if (exists $gene2tools->{$gene_id}) {
        my $tool2blocks = $gene2tools->{$gene_id};
        foreach my $tool (sort keys %$tool2blocks) {
            my $blocks = $tool2blocks->{$tool};
            $y = int($y_top + ($i + 1) * $y_inter);
            $y1 = $y + $thick / 2;
            $y2 = $y + $thick;
            # $svg->text("x" => $x_str_left, "y" => $y2, "font-size" => $thick, 
            #             "-cdata" => $tool);
            $svg->line("x1" => $x_left, "y1" => $y1, 
                        "x2" => $x_left + $seq2len->{$gene_id} / $scale, 
                        "y2" => $y1, "stroke" => "black");
        
            foreach my $block (@$blocks) {
                my $domain_id = $block->{"motif_id"};
                my $x = $x_left + int($block->{"start"} / $scale);      
                my $width = $block->get_len() / $scale;
                $svg->rectangle("x" => $x, "y" => $y, "width" => $width, 
                                "height" => $thick, "fill" => "gray");       
                $svg->text("x" => $x,"y" => $y2, "-cdata" => $domain_id, 
                            "font-family" => "Times", "font-size", int($thick *8 / 10),
                            "fill" => "white");
            }
            $i++;
        }
    }
    $i=$i+1;
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

print $svg->xmlify;


