#!/usr/bin/perl

# miRDeep2 convert-bowtie-output perl script
# Copyright (C) 2009 - 2011  Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
