package util;

use strict;
use warnings;
use Text::Trim;

sub load_ids_array {
    my $file = shift;
    my @ids;
    open FILE, "<$file" or die $!;
    while (my $line = <FILE>) {
		trim $line;
        next if $line =~ /^\s*$/;
        push @ids, $line;
    }
    close FILE;
    return \@ids;
}


sub load_ids_hash {
    my $file = shift;
    my %ids;
    open FILE, "<$file" or die $!;
    while (my $line = <FILE>) {
        trim $line;
        next if $line =~ /^\s*$/;
        $ids{$line} = 1;
    }
    close FILE;
    return \%ids;
}

sub load_map_id {
    my $file = shift;
    my %id2name;
    open FILE, "<$file" or die $!;
    while (my $line = <FILE>) {
        trim $line;
        my @items = split /\t/, $line;
        $id2name{$items[0]} = $items[1] if @items > 1;
    }
    close FILE;
    return \%id2name;
}     
            
sub rev_comp {
    my $seq = shift;
    $seq = reverse $seq;
    $seq =~ tr/ATGCatgc/TACGtacg/;
    return $seq;
}

sub intersect_all {
	my $array = shift;
	my $size = scalar @$array;	
	my $ids = $array->[0];
	for (my $i = 1; $i < $size; $i++) {
		$ids = intersect($ids, $array->[$i]);
	}
	return $ids;
}

sub intersect {
	my ($ids1, $ids2) = @_;
	my %ids;
	foreach my $id1 (keys %$ids1) {
		if (exists $ids2->{$id1}) {			
			$ids{$id1} = 1;
		}
	}
	return \%ids;
}

1;
