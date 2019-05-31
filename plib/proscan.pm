package proscan;

use strict;
use warnings;
use block;

#AT2G14080  496d7c8e064019795e6d0f8edbde3915    1215    Pfam    PF01582 TIR domain  55  233 2.0E-55 T   24-02-2016
sub new {
    my ($class, $file) = @_;
    my $self = bless { }, $class;
    $self->{"scan_file"} = $file;
    $self->{"gene2blocks"} = { };
    $self->{"gene2tools"} = { };
    return $self;
}

sub load_gene2tools {
    my $self = shift;
    open FILE, "<$self->{scan_file}" or die $!;
    while (my $line = <FILE>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;
        my @items = split /\t/, $line;
        my $gene_id = $items[0];
        my $tool = $items[3];
        my $motif_id = $items[4];
        my $start = $items[6];
        my $end = $items[7];
        my $domain_len = $end - $start + 1;
        my $gene2tools = $self->{"gene2tools"};
        $motif_id = "TMhelix" if $tool eq "TMHMM";
        my $block = block->new(motif_id => $motif_id, gene_id => $gene_id, start => $start, end => $end);
        if (exists $gene2tools->{$gene_id}) {
            my $tool2blocks = $gene2tools->{$gene_id};
            if (exists $tool2blocks->{$tool}) {
                push @{$tool2blocks->{$tool}}, $block;
            } else {
                $tool2blocks->{$tool} = [ $block ];
            }
        } else {
            $gene2tools->{$gene_id} = { $tool => [ $block ] }
        }
    }
    close FILE;
    return
}

sub load_gene2blocks {
    my $self = shift;
    open FILE, "<$self->{scan_file}" or die $!;
    while (my $line = <FILE>) {
        chomp $line;
        next if $line =~ /^\s*$/;
        next if $line =~ /^#/;
        my @items = split /\t/, $line;
        my $gene_id = $items[0];
        my $motif_id = $items[4];
        my $start = $items[6];
        my $end = $items[7];
        my $e_value = $items[8];
        my $block = block->new(motif_id => $motif_id, gene_id => $gene_id, start => $start, end => $end, e_value =>
            $e_value);
        my $gene2blocks = $self->{gene2blocks};
        if (exists $gene2blocks->{$gene_id}) {
            push @{$gene2blocks->{$gene_id}}, $block;
        } else {
            $gene2blocks->{$gene_id} = [ $block ];
        }
    }
    close FILE;
}

1;
