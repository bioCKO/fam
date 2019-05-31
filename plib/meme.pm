package meme;

use strict;
use warnings;
use block;
use motif;

sub new {
	my ($class, $meme_file) = @_;
	my $self = bless {}, $class;
    $self->{"id2motif"} = {};
	open FILE, "<$meme_file" or die $!;
	my $i = 0;
    while (my $line = <FILE>) {
		chomp $line;
		next if $line =~ /^\s*$/;
		#Motif AKTGEKEWYFFSPRDRKYPNGSRTNRATAAGYWKATGKDKAI MEME-2 in BLOCKS format
		#--------------------------------------------------------------------------------
		#BL   MOTIF AKTGEKEWYFFSPRDRKYPNGSRTNRATAAGYWKATGKDKAI width=42 seqs=50
		#Pbr030137                (   66) AKMGEKEWYFFCVRDRKYPTGLRTNRATEAGYWKATGKDKEI  1 

		if ($line =~ /Motif (\S+) MEME-(\d+) in BLOCKS format/) {
			my $seq = $1;
	        my $motif_id = $2;
	        
            my $motif = motif->new("motif_id" => $motif_id, "regex"=> $seq);
            $self->{"id2motif"}->{$motif_id} = $motif;
			while ($line = <FILE>) {
				chomp $line;
				next if $line =~ /^\s*$/;
				if ($line =~ /^(\S+) +\( +(\d+)\) +(\S+)/) {
                	my $gene_id = $1;
                    my $start = $2;
		            my $seq = $3;
                    my $end = $start + length($seq) - 1;
                    my $block = block->new("motif_id" => $motif_id, "gene_id" => $gene_id,
											"start" => $start, "end" => $end, "seq" => $seq);
                    $motif->add_block($block);
				}
				
				last if $line =~ /^\/\//;
			}
		}	
	}
	close FILE;
	return $self;
}
               

sub get_gene2blocks {
	my $self = shift;
	my %gene2blocks_all;
	foreach my $motif_id (keys %{$self->{"id2motif"}}) {
		my $motif = $self->{"id2motif"}->{$motif_id};
		my $gene2blocks = $motif->{"gene2blocks"};
		foreach my $gene_id (keys %{$gene2blocks}) {
			my $blocks = $gene2blocks->{$gene_id};
			foreach my $block (@$blocks) {
				if (exists $gene2blocks_all{$gene_id}) {
                	push @{$gene2blocks_all{$gene_id}}, $block;
                } else {
                    $gene2blocks_all{$gene_id} = [$block];
                }
            }
        }
    }
	return \%gene2blocks_all;
}

sub get_motif_ids {
	my $self = shift;
    my @ids = keys %{$self->{"id2motif"}};
	return \@ids;
}
  
1; 
