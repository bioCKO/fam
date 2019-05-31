package motif;

use strict;
use warnings;


sub new {
	my ($class, %arg) = @_;
	my $self = bless {}, $class;
    $self->{"motif_id"} = $arg{"motif_id"};
    $self->{"regex"} = $arg{"regex"};
    $self->{"gene2blocks"} = {};
	return $self;
}

sub add_block {
	my ($self, $block) = @_;
	my $gene2blocks = $self->{"gene2blocks"};
	my $gene_id = $block->{"gene_id"};
    if (exists $gene2blocks->{$gene_id}) {
        push @{$gene2blocks->{$gene_id}}, $block;
	} else {
		$gene2blocks->{$gene_id} = [$block];
	}
}

1;