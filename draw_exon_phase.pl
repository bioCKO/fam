#!/usr/bin/env perl

use warnings;
use util;
use gff;
use proscan;
use SVG;

my $gene_ids = util::load_ids_array($ARGV[0]);
my $scan = proscan->new($ARGV[1]);
my $gff = gff->new($ARGV[2]);

$scan->load_gene2blocks();
my $gene2blocks = $scan->{gene2blocks};

foreach my $gene_id (@$gene_ids) {
    my $blocks = $gene2blocks->{$gene_id} if exists $gene2blocks->{$gene_id};
    foreach my $block (@$blocks) {
        $block->{start} = $block->{start} * 3 - 3;
        $block->{end} = $block->{end} * 3 - 1;
    }
}

my $gene2rna = $gff->{gene2long_rna};
my $rna2cdss = $gff->{rna2cdss};

my %hmm2color = ('PF02365' => 'blue');

my %gene2cdss;
my %color_block;
my $max_len = 0;

foreach my $gene_id (@$gene_ids) {
    if (! exists $gene2rna->{$gene_id}) {
        die "$gene_id has no rna";
    }
    my $rna = $gene2rna->{$gene_id};
    my $rna_id = $rna->[8];
    my $strand = $rna->[6];
    my $cdss = $rna2cdss->{$rna_id};

    if ($strand eq '+') {
        my $off_set = $cdss->[0]->[3];
        foreach my $cds (@$cdss) {
            $cds->[3] = $cds->[3] - $off_set + 1;
            $cds->[4] = $cds->[4] - $off_set + 1;
        }
    } else {
        my $off_set = $cdss->[scalar(@$cdss) - 1]->[4];
        foreach my $cds (@$cdss) {
            my $start = $cds->[3];
            my $end = $cds->[4];
            $cds->[3] = $off_set - $end + 1;
            $cds->[4] = $off_set - $start + 1;
        }
    }

    my @sorted_cdss = sort {$a->[3] <=> $b->[3]} @$cdss;

    $max_len = $sorted_cdss[$#sorted_cdss]->[4] if $sorted_cdss[$#sorted_cdss]->[4] > $max_len;

    my @cds_pos;
    foreach my $cds (@sorted_cdss) {
        push @cds_pos, ($cds->[3] .. $cds->[4]);
    }

    my $blocks = $gene2blocks->{$gene_id};
    foreach my $block (@$blocks) {
        my $motif_id = $block->{motif_id};
        if (exists $hmm2color{$motif_id}) {
            my $color = $hmm2color{$motif_id};
            my $pre = $cds_pos[$block->{"start"}];
            my $start = $pre;
            for (my $i = $block->{start} + 1; $i <= $block->{end}; $i++) {
                if ($cds_pos[$i] - $pre == 1) {
                    $pre = $cds_pos[$i];
                } else {
                    if (exists $color_block{$gene_id}) {
                        push @{$color_block{$gene_id}}, [ $start, $pre, $color ];
                    } else {
                        $color_block{$gene_id} = [ [ $start, $pre, $color ] ];
                    }
                    $pre = $cds_pos[$i];
                    $start = $pre;
                }
            }

            if (exists $color_block{$gene_id}) {
                push @{$color_block{$gene_id}}, [ $start, $pre, $color ];
            } else {
                $color_block{$gene_id} = [ [ $start, $pre, $color ] ];
            }
        }
    }

    my $exon_len = 0;
    for (my $i = 0; $i < scalar(@sorted_cdss) - 1; $i++) {
        my $cds = $sorted_cdss[$i];
        $exon_len += $cds->[4] - $cds->[3] + 1;
        $cds->[7] = $exon_len % 3;
    }

    $gene2cdss{$gene_id} = \@sorted_cdss;
}

my $str_x_off = 10;

my $thick = 10;
my $y_inter = $thick + $thick;
my $x_left = 200;
my $y_top = 100;
my $axis_y = $y_top - 20;
my $bar_width = 800;

my $scale = $max_len / $bar_width;
my $id_num = scalar(@$gene_ids);
my $im_width = $bar_width + $x_left * 2;
my $im_heigth = $y_inter * $id_num + 2 * $y_top;

my $svg = SVG->new(width => $im_width, height => $im_heigth);

$svg->line(x1      => $x_left, y1 => $axis_y, x2 => $x_left + $max_len / $scale, y2 => $axis_y, stroke => "black",
    'stroke-width' => 1);

for (my $i = 0; $i < $max_len; $i = $i + 100) {
    $svg->line(x1               => $x_left + $i / $scale, y1 => $axis_y - 5, x2 => $x_left + $i / $scale, y2 => $axis_y,
        stroke                  =>
        "black", 'stroke-width' => 0.5);
    if ($i % 2000 == 0) {
        my $num = $i / 1000;
        $svg->line(x1 => $x_left + $i / $scale, y1 => $axis_y - 10, x2 => $x_left + $i / $scale, y2 => $axis_y,
            stroke    => "black", 'stroke-width' => 0.5);
        $svg->text(x    => $x_left + $i / $scale, y => $axis_y - 20, font => "Arial",
            'font-size' => $thick, '-cdata' => "${num}k");
    }
}

my $i = 0;
foreach my $gene_id (@$gene_ids) {
    my $y = $y_top + $i * $y_inter;
    my $y1 = $y - $thick / 2;
    my $cdss = $gene2cdss{$gene_id};
    $svg->text(x => $str_x_off, y => $y + $thick / 2, 'font-size' => 8, '-cdata' => $gene_id);
    my $end = $cdss->[scalar(@$cdss) - 1]->[4];

    $svg->line(x1 => $x_left, y1 => $y, x2 => $x_left + $end / $scale, y2 => $y, stroke => "black", 'stroke-width' =>
        1);

    for (my $j = 0; $j < scalar(@$cdss); $j++) {
        my $cds = $cdss->[$j];
        my $x = $x_left + $cds->[3] / $scale;
        my $width = ($cds->[4] - $cds->[3] + 1) / $scale;

        $svg->rectangle(x => $x, y => $y1, width => int($width), height => $thick, fill => "gray");

        if ($j != scalar(@$cdss) - 1) {
            $svg->text(x => $x + $width + 1, y => $y - 1, font => "Arial", 'font-size' => $thick, '-cdata' =>
                $cds->[7]);
        }
    }

    my $domains = $color_block{$gene_id};
    foreach my $exon (@$domains) {
        my $x = $x_left + $exon->[0] / $scale;
        my $width = ($exon->[1] - $exon->[0] + 1) / $scale;
        $svg->rectangle(x => $x, y => $y1, width => $width, height => $thick, fill => $exon->[2]);
    }

    $i++;
}

print $svg->xmlify;
