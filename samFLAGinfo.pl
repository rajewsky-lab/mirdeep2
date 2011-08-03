#!/usr/bin/perl


## convert decimal number to binary and then to hexadecimal to get read FLAG identification 

use strict;

if(not $ARGV[0]){
    die "\nPlease specify the FLAG number you want to have translated\n
usage:\n\tperl samFLAGinfo.pl [int]\n\n";
} 


my %bwa_codes;
while(<DATA>){
    chomp;
    if(/^(\d+)\s+(.+)$/){
        $bwa_codes{$1} = $2;
        #print "$1\t$2\n";
    }
}



my $in=$ARGV[0];

my $rest;

$rest = $in;

my @arr;

while($rest ne 0){
    push(@arr,$rest%2);
    $rest = int($rest/2);
}



my $hex;
my $bin;

print "\n";

my $pre="0x0";

for(my $i=0; $i < scalar @arr; $i++){
    $bin = $arr[$i] * 2**$i;
    $hex = sprintf("%x", $bin);
    $pre="0x0";  
    if($arr[$i] ne 0){
        $pre.='0' x (3-length($hex));
        $pre.=$hex;
        print "$arr[$i] * 2^$i = ",$bin,"\thex = $pre \t $bwa_codes{$hex}\n";
    }
}

print "\n";

__DATA__
0   .
1	the read is paired in sequencing
2	the read is mapped in a proper pair
4	the read sequence is unmapped
8	the mate is unmapped
10	read is mapped to minus strand (given seq in col 10 is therefore the reverse complement of the plus strand)
20	strand of the mate
40	the read is the first read in a pair
80	the read is the second read in a pair
100	the alignment is not primary
200 the read fails plattform/vendor quality checks
400 the read is either a PCR duplicate or an optical duplicate
