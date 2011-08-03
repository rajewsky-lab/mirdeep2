#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $usage =
"$0 file_fasta

finds the total read count if '_x\\d' format is used.
";

my $file_fasta=shift or die $usage;

my %options=();
getopts("",\%options);


my $count_reads=0;

open (FILE, "<$file_fasta") or die "can not open $file_fasta\n";

while (my $line=<FILE>){

    if($line=~/_x(\d+)/){
	$count_reads+=$1;
    }
}

print "$count_reads\n";







