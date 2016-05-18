#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $usage =
"$0 _seq.txt

This script parses _seq.txt to fasta format. Options:
-a    format is qseq.txt

Example of use:
$0 mouse_seq.txt > mouse.fa
";

my $file_solexa=shift or die $usage;

my %options=();
getopts("a",\%options);

my $running=0;



parse_file_solexa($file_solexa);

exit;




sub parse_file_solexa{

    my ($file) = @_;

    open (FILE, "<$file") or die "can not open $file\n";
    while (<FILE>){
		next if(/^\s*$/);
	my @fields=split(/\s+/);

	my $seq;

	if($options{a}){

	    $seq=$fields[8];

	}else{	

	    $seq=$fields[4];
	}

	chomp $seq;

	$seq=~s/\./N/g;

	print ">seq\_$running\n";
	
	print "$seq\n";

	$running++;
    }
    return;
}

