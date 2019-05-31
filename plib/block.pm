package block;

use strict;
use warnings;

sub new {
    my ($class, %arg) = @_;
    my $self = bless {}, $class;
    $self->{"motif_id"} = $arg{"motif_id"};
    $self->{"gene_id"} = $arg{"gene_id"};
    $self->{"start"} = $arg{"start"};
    $self->{"end"} = $arg{"end"};
    $self->{"seq"} = $arg{"seq"} if exists $arg{"seq"};
    $self->{"e_value"} = $arg{"e_value"} if exists $arg{"e_value"};
    return $self;
}

sub get_len {
    my $self = shift;
    return $self->{"end"} - $self->{"start"} + 1;
}

1;
