package fasta;

use strict;
use warnings;
use Bio::SeqIO;

sub new {
    my ($class, $file) = @_;
    my $self = bless {}, $class;
    $self->{"id2seq"}  = {};
    $self->{"ids"}     = [];
    $self->{"read_num"} = 0;
    $self->{"nt_num"}   = 0;
    $self->read_fasta_file($file);
    return $self;
}

sub read_fasta_file {
    my ($self, $file) = @_;
    my $seq_io = Bio::SeqIO->new(-file => $file, -format => "fasta");
    while (my $seq_obj = $seq_io->next_seq) {
        my $id  = $seq_obj->id;
        my $seq = $seq_obj->seq;
        $seq =~ s/\s//g;
        $self->{"nt_num"} += length($seq);
        $self->{"read_num"}++;
        push @{$self->{"ids"}}, $id;
        $self->{"id2seq"}->{$id} = $seq;
    }
}

sub get_seqs_by_ids {
    my ( $self, $ids ) = @_;
    my @seqs;
    foreach my $id (@$ids) {
        my $seq = $self->{"id2seq"}->{$id};
        push @seqs, \$seq;
    }
    return \@seqs;
}

sub get_ids {
    my $self = shift;
    return $self->{"ids"};
}

sub get_seqs {
    my $self = shift;
    my @v    = values %{ $self->{"id2seq"} };
    return \@v;
}

sub get_id2seq_len {
    my ( $self, $ids ) = @_;
    my %id2seq_len;
    my $id2seq = $self->{"id2seq"};
    if ($ids) {
        foreach my $id (@$ids) {
            $id2seq_len{$id} = length( $id2seq->{$id} ) if exists $id2seq->{$id};
        }
    }
    else {
        foreach my $id ( sort keys %$id2seq ) {
            $id2seq_len{$id} = length( $id2seq->{$id} );
        }
    }
    return \%id2seq_len;
}

1;
