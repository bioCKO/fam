package gff;

use strict;
use warnings;
use Text::Trim;

sub new {
    my ($class, $file) = @_;
    my $self = bless {}, $class;
    $self->{chr2genes} = {};
    $self->{gene2rnas} = {};
    $self->{gene2long_rna} = {};
    $self->{rna2cdss} = {};
    $self->{rna2utrs} = {};
    $self->{rna2exons} = {};
    $self->{gene_id2gene} = {};
    $self->read_gff($file);
    return $self;
}

sub read_gff {
    my ($self, $file) = @_;
    open FILE, "<$file" or die $!;
    while (my $line = <FILE>) {
        trim $line;
        next if $line =~ /^#/;
        my @items = split /\t/, $line;
        next if @items != 9;
        my $chr_id = $items[0];
        my $source = $items[1];
        my $type = $items[2];
        my $start = $items[3];
        my $end = $items[4];
        my $score = $items[5];
        my $strand = $items[6];
        my $phase = $items[7];
        my $attr = $items[8];
        print STDERR "Error position: $line\n" if $start > $end;
        if ($type eq "gene") {
            if ($attr =~ /^ID=([^;]+)/) {
                my $gene_id = $1;
                my @gene = @items[0..7];
                push @gene, $gene_id; 
                $self->{gene_id2gene}->{$gene_id} = \@gene;
                if (exists $self->{chr2genes}->{$chr_id}) {
                    push @{$self->{chr2genes}->{$chr_id}}, \@gene;
                } else {
                    $self->{chr2genes}->{$chr_id} = [\@gene]
                }
             } else {
                print STDERR "Error gene: $line\n";
             }
        }

        if ($type eq "mRNA") {
			if ($attr =~ /^ID=([^;]+);Name=[^;]+;Parent=([^;]+)$/) {
				my $rna_id = $1;
				my $gene_id = $2;
				my @rna = @items[0..7];
				push @rna, $rna_id;
				if (exists $self->{gene2rnas}->{$gene_id}) {
					push @{$self->{gene2rnas}->{$gene_id}}, \@rna;
				} else {
					$self->{gene2rnas}->{$gene_id} = [\@rna];
				}
			} else {
                print STDERR "Error mRNA: $line\n";
            }
		}

        if ($type eq "CDS") {
            if ($attr =~ /^Parent=([^;]+)$/) {
                my $rna_id = $1;
                my @cds = @items[0..7];
                push @cds, $rna_id;
                if (exists $self->{rna2cdss}->{$rna_id}) {
                    push @{$self->{rna2cdss}->{$rna_id}}, \@cds;
                } else {
                    $self->{rna2cdss}->{$rna_id} = [\@cds];
                }
                if (exists $self->{rna2exons}->{$rna_id}) {
                    push @{$self->{rna2exons}->{$rna_id}}, \@cds;
                } else {
                    $self->{rna2exons}->{$rna_id} = [\@cds];
                }
            } else {
                print STDERR "Error CDS: $line\n";
            }
        }
        
		if ($type eq "UTR") {
            if ($attr =~ /^Parent=([^;]+)$/) {
                my $rna_id = $1;
                my @exon = @items[0..7];
                push @exon, $rna_id;
                if (exists $self->{rna2exons}->{$rna_id}) {
                    push @{$self->{rna2exons}->{$rna_id}}, \@exon;
                } else {
                    $self->{rna2exons}->{$rna_id} = [\@exon];
                }
            } else {
                print STDERR "Error exon: $line\n";
            }
        }
	} #end while
	close FILE;
     
    foreach my $chr_id (keys %{$self->{chr2genes}}) {
        my $genes = $self->{chr2genes}->{$chr_id};
        my @sorted = sort {$a->[3]<=>$b->[3]} (@$genes);
        $self->{chr2genes}->{$chr_id} = \@sorted;
    }
        
    foreach my $gene_id (keys %{$self->{gene2rnas}}) {
        my $rnas = $self->{gene2rnas}->{$gene_id};
        my @sorted = sort {$a->[3]<=>$b->[3]} (@$rnas); 
        $self->{gene2rnas}->{$gene_id} = \@sorted;
    }
      
    foreach my $rna_id (keys %{$self->{rna2cdss}}) {
        my $cdss = $self->{rna2cdss}->{$rna_id};
        my @sorted = sort {$a->[3]<=>$b->[3]} (@$cdss); 
        $self->{rna2cdss}->{$rna_id} = \@sorted;
    }
    
    foreach my $rna_id (keys %{$self->{rna2exons}}) {
        my $exons = $self->{rna2exons}->{$rna_id};
        my @sorted = sort {$a->[3]<=>$b->[3]} (@$exons);
        $self->{rna2exons}->{$rna_id} = \@sorted;
    }
    
    $self->ex_longest_rna();
}

sub ex_longest_rna {
    my $self = shift;
    foreach my $gene_id (keys %{$self->{gene2rnas}}) {
        my $rnas = $self->{gene2rnas}->{$gene_id};
        my $max_len = 0;
        my $long_rna;
        foreach my $rna (@$rnas) {
            my $rna_id = $rna->[8];
            my $cdss = $self->{rna2cdss}->{$rna_id};
            my $gene_len = 0;
            foreach my $cds (@$cdss) {
                my $cds_len = $cds->[4] - $cds->[3] + 1;
                $gene_len += $cds_len;
            }
         
            if ($gene_len > $max_len) {
                $max_len = $gene_len;
                $long_rna = $rna;
            }
        }
        $self->{gene2long_rna}->{$gene_id} = $long_rna;
    }
}

1;
