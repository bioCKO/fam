#!/usr/bin/perl

use strict;
use warnings;

use Bio::AlignIO;

my $in=Bio::AlignIO->new(-format=>"fasta",-file=>$ARGV[0]);
my $out=Bio::AlignIO->new(-format=>"nexus",-file=>">$ARGV[1]", -show_symbols=>0,-show_endblock=>0);
#my $io_in = Bio::AlignIO->new(-format=>"phylip", -file=>$ARGV[0], -idlength=>5, -interleaved=>1);
#my $io_out = Bio::AlignIO->new(-format=>'fasta', -file=>">$ARGV[1]");

if(my $aln=$in->next_aln){
    $out->write_aln($aln);
    $out->_print("BEGIN MRBAYES;
lset rates=gamma;
prset aamodelpr=fixed(wag);
mcmc ngen=1000000 samplefreq=100 nchains=4 nruns=2 stopval=0.01 stoprule=yes;
sumt burnin=100;
sump burnin=100;
END;");
}else{
    warn("Cannot find any alignments in file.\n");
}
