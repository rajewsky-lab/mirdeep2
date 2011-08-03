#!/usr/bin/perl

use strict;

my $file = $ARGV[0] or die "Usage: \n\n $0 fastq file > outputfile\n";;

open IN,"<$file" or die "File $file not found error here\n";
my $c = 0;
my @line = ();
my $id ='';
my $seq = '';

my $processed = 0;

while(<IN>){
    chomp;
    $c++;
    if($c == 1){
        $processed++;
#        print STDERR "$processed reads processed\r";
        @line = split();
        $id = $line[0];
        $id =~ s/\@//;
        
    }elsif($c == 2){
        $seq = $_;
    }elsif($c == 4){
        $c = 0;
        print ">seq_$processed\n$seq\n";
        $id ='';
        @line =();
        $seq ='';
    }else{}
}
close IN;
exit;
#print STDERR "$processed reads processed\n\n";
