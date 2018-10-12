#!/usr/bin/perl

# miRDeep2 SAM-reads-collapse perl script
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


if(not $ARGV[1]){
	die "
Usage:  $0  sam_file  reads_output_file

This script parses a sam file and collapses
reads. The reads_output_file can later on be
used as input for miRDeep2.

";
}

open IN,"<$ARGV[0]" or die "Sam file not found\n";



while(<IN>){
    next if(/^@/);
    @line=split();

    $foo=$line[1];
    $rev= 0x10;
    $strand = '+';
    if(($foo & $rev) == $rev){$strand = '-';}
        if($strand eq "-"){
            $line[9] = reverse($line[9]);
            $line[9] =~ tr/ACGTUacgtu/TGCAATGCAA/;
        }
    if(not $seq{$line[9]}){
    $seq{$line[9]}="$line[0]";
}else{
    $seq{$line[9]}.=",$line[0]";

}
}
close IN;

open OUT,">$ARGV[1]" or die "File $output could not be created\n";
open OUT2,">reads_N_to_1.txt" or die "File reads_N_to_1.txt could not be created\n";
open OUT3,">read_1_to_1.txt" or die "file not created\n";
my $c=0;

my $pref='seq';
my $css=0;

for my $k( keys %seq){
    @ids=split(",",$seq{$k});
    #$c=($#ids)
    $finalid='';
    if($ids[0] =~ /(\S+)_x(\d+)$/){
        $c = $2;
        $finalid=$1;
    }else{
        $finalid=$ids[0];
        $c=1;
    }
    foreach my $i(@ids[1..$#ids]){
    if($ids[$i] =~ /_x(\d+)$/){
        $c += $1;
    }else{
        $c++;
    }
        print OUT2 "$i,$ids[0]\n";
    }
    $css++;
    print OUT ">${pref}_${css}_x$c\n$k\n";
    print OUT3 "$ids[0],${pref}_${css}_x$c\n";
}
close OUT;
close OUT2;
close OUT3;
