#!/usr/bin/perl

use strict;

my $id;
my $seq;
my @species = split(",",$ARGV[1]);
die "Error: Specify species or comma separated list of species\n\nUsage:\n\tperl $0 miRBase_fasta_file species [mature/star]\n\te.g. $0 mature.fa cel mature\n\n" if(not $ARGV[1]);
my $in = 0;

open IN,"<$ARGV[0]" or die "$ARGV[0] not found $!";

my $first =1;

my $str;

my $sp;

my $type='mature';
if($ARGV[2]){
	if(lc($ARGV[2]) eq 'mature' or lc($ARGV[2]) eq 'star'){
		$type = $ARGV[2];
	}else{
		die 'argv 2 must either be mature or star\n';
	}
}

foreach my $sp(@species){

$in = 0;

open IN,"<$ARGV[0]" or die "$ARGV[0] not found $!";

$first =1;

$str;

my $id;

#die $type;

while(<IN>){
    chomp;
    if(/(>$sp\S+)/){
		$id=$1;
		if($type eq 'mature' and $id =~/\*$/){ $in= 0; next;}
		if($type eq 'star' and $id !~ /\*/){$in=0; next;}
	#print "=== $id\n";
       if($first){
            print "$id\n";
            $first = 0;
        }else{
            print "\n$id\n";
        }
        $in =1;
        next;
     }elsif(/>/){
         $in = 0;
     }elsif($in == 1){
	#print "\n==== $_ ===\n";
       $str = uc($_);
       $str =~ tr/U/T/;
       print $str;
     }else{
	$id='';
	$str='';
	$in=0;
     }
}

if(not $first){
    print "\n";
}
close IN;
}
