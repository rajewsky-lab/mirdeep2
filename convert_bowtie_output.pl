#!/usr/bin/perl

use strict;

my @line;
my @gseq;
my @changes;
my $pos;
my $ont;
my $mm=0;
my @edit="";
my $reverse = 0;


open IN,"<$ARGV[0]" or die "usage: $0 reads_mapped.bwt\n";

while(<IN>){
    @line = split(/\t/);
    if($line[1] eq "-" ){
        $line[4] = reverse($line[4]);
        $line[4] =~ tr/ACGTN/TGCAN/;
	}
    @gseq = split(//,lc $line[4]);
    @edit= split(//,("m" x length($line[4])));
    $mm=0;



    if($line[7]){
		
        @changes = split(/,/,$line[7]);
        #$mm=scalar @changes;
        foreach(@changes){
            if(/(\d+):(\w+)\>\w+/){
				$mm++;
				$gseq[$1] = lc $2;
				## this was not outcommented before
                #$gseq[$1] =~ tr/acgt/tgca/;
                @edit[$1] = "M";   
            }
        }
    }
 

    my @id =split(/\s/,$line[0]);
    my @db = split(/\s/,$line[2]);

    print "$id[0]\t",length($line[4]),"\t1\t",length($line[4]),"\t",lc $line[4],"\t$db[0]\t",length($line[4]),"\t",$line[3]+1,"\t",$line[3]+length($line[4]),"\t",@gseq,"\t",$line[1],"\t",$mm,"\t",@edit,"\n";
}
close IN;
