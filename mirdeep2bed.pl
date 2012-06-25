#!/usr/bin/perl 

use strict;

open IN,"$ARGV[0]" or die "No csv file given for bed conversion\n";

my ($known,$novel, $not, $line,$thres,$score,$line,@l,$strand,$label);

while(<IN>){
	if(/novel miRNAs predicted by miRDeep2/){
		$novel=1;
		$known=0;
		$not=0;
	}elsif(/mature miRBase miRNAs detected/){
		$novel=0;
		$known=1;
		$not=0;
	}elsif(/miRBase miRNAs not detected/){
		last;
	}else{
		@l=split();
		if($l[$#l] =~ /(\S+):(\d+)\.\.(\d+):(\S)/){
			$strand='255,0,0' if($4 eq '+');
			$strand='0,0,255' if($4 eq '-');
			$label="known" if($known);
			$label="novel" if($novel);
			print "$1\t$2\t$3\t${label}:$l[0]\t$l[1]\t$4\t$2\t$3\t$strand\n";
		}
	}
}
