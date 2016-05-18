#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;

my $usage =
"$0 fastafile

This script parses a fastafile.

-a int    only output entries where the sequence is minimum int nts long
-b        remove all entries that have a sequence that contains letters
          other than a,c,g,t,u,n,A,C,G,T,U,N.
-s        output progress
";

my $file_fasta=shift or die $usage;

my %options=();
getopts("a:bs",\%options);


my $running=0;


parse_fasta($file_fasta);


sub parse_fasta{

    my ($file_fasta) = @_;
    my ($id,$seq) = ();

    open (FASTA, "<$file_fasta") or die "can not open $file_fasta\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)/)
        {
            $id  = $1;
	    $seq = "";
            while (<FASTA>){
                chomp;
                if (/^>(\S+)/){
		    resolve($id,$seq);
                    $id    = $1;
                    $seq   = "";
                    next;
                }
                $seq .= $_;
            }
        }
    }
    resolve($id,$seq);
    close FASTA;
    return;
}



sub resolve{

    my($id,$seq)=@_;

    $running++;

    if($options{s}){print STDERR "$running\r";}

    my $lng=length $seq;

    if($options{a} and $lng<$options{a}){
        print STDERR ">$id\n$seq\n";
	return;
    }

    if($options{b} and not($seq=~/^(a|c|g|t|u|n)+$/i)){
	
        print STDERR ">$id\n$seq\n";
	return;
    }

    print ">$id\n";
	
    print "$seq\n";

    return;
}

















