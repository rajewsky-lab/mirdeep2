#!/usr/bin/env perl

# miRDeep2 get-precursors perl script
# Copyright (C) 2012 - 2014  Sebastian Mackowiak
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

#######################################################################
#
#
# This script extracts the mirdeep2 precursors from the result.csv file
# and writes a fasta file and a bed file for these
# See 'Usage' for more information.
#
# S.M. Jan. 2019
#
######################################


use strict;
use Getopt::Std;


my %options=();
getopts("r:s:dpht:mkT:o:",\%options);
if(not -f $options{'r'} or $options{'h'}){
	die "
Usage: get_mirdeep2_precursors.pl -r result.csv  [options]

[options]
	-s [int]\toutput only precursors with min-score >= [int]
	-t [int]\toutput only precursors with score     <  [int]
	-d      \toutput dna instead of rna
	-p      \tmake simple id with the name only
	-h      \tprint this help message
	-m      \tget_mature instead of precursor
	-k      \tget_star instead of precursor
	-T      \tTrackname for bedfiles
	-o      \toutdir
";
}

my $od="result_mirnas";
$od=$options{'o'} if($options{'o'});

if(not -d $od){mkdir $od;}

if($options{'T'} !~ /\S+/){
$options{'T'} = 'notTrackname';
}

my @l1=split(/result/,$options{'r'});
my ($timestamp,$dmp)=split(".csv",$l1[1]);


my ($known,$novel, $not, $line,$thres,$score,$strandcol,$bedh1,$pcoord);

my $bedh="browser position BPOS
browser hide all
track name=\"TNAME\" description=\"TDESC\" visibility=2
itemRgb=\"On\";
";
my $first=1;


my @line;

$thres = -50;

if(defined $options{'s'}){ $thres = $options{'s'};}
$score=$thres;
my $max='na';
my $maxs=999999999999999999999999999;
if(defined $options{'t'}){
	$max=$options{'t'};
	$maxs=$max;
}
open IN,"$options{'r'}" or die "No results.csv file given\n";

my $seqcol=15;
$seqcol=13 if($options{'m'});
$seqcol=14 if($options{'k'});

my @names=qw(0 1 2 3 4 5 6 7 8 9 10 11 12 mature star pres);

while(<IN>){
	if(/novel miRNAs predicted by miRDeep2/){
		$novel=1;
		$known=0;
		$not=0;
		$line=<IN>;
        $first=1;
		open OUT,">$od/novel_$names[$seqcol]${timestamp}_score${score}_to_$max.fa";
		open BED,">$od/novel_$names[$seqcol]${timestamp}_score${score}_to_$max.bed";
		next;
	}elsif(/miRBase miRNAs detected/){
        close OUT;
        close BED;
		$novel=0;
		$known=1;
		$not=0;
		$line=<IN>;
        $first=1;
		close OUT;
        close BED;
		open OUT,">$od/known_$names[$seqcol]${timestamp}_score${score}_to_$max.fa";
        open BED,">$od/known_$names[$seqcol]${timestamp}_score${score}_to_$max.bed";
		next;
	}elsif(/miRBase miRNAs not detected by miRDeep2/){
		$novel=0;
		$known=0;
		$not=1;
        $first=1;
		close OUT;
        close BED;
		open OUT,">$od/not_$names[$seqcol]${timestamp}_score${score}_to_$max.fa";
		next;
	}else{}
	next if(/^\s*$/);
			

	if($novel or $known){
		chomp;
		@line=split(/\t/);
		my $coord='na';
		$coord=$line[16] if($line[16]);
		my $c5p = index($line[15],$line[13]);
		my $c3p = index($line[15],$line[14]);
		my $l5p = length($line[13]);
		my $l3p = length($line[14]);

		if($known){
			if($line[1] >= $thres and $line[1] < $maxs){
				if($options{'d'}){
					$line[$seqcol] =~ tr/uU/tT/;
				}
				if($options{'p'}){
					if($line[0] =~ /\|([a-zA-Z0-9_-]*)$/){
						$line[0] = $1;
					}
					print OUT ">$line[0]\n",uc($line[$seqcol]),"\n";

				}else{
					print OUT ">${line[0]}_${line[9]}_x${line[5]}_coord:",$coord,"_score:$line[1]\n",uc($line[$seqcol]),"\n";
				}
                chomp $coord;

                $strandcol="0,0,255";
                if($coord =~ /^(\S+):(\d+)\.\.(\d+):(\S)$/){
					my ($c,$b,$e,$s) = ($1,$2,$3,$4);
					if($options{'m'}){
						if($s eq '+'){
							$b+=$c5p;
							$e=$b+$l5p;
						}else{
							$b=$e-$c5p-$l5p;
							$e=$e-$c5p;
						}
					}
					if($options{'k'}){
						if($s eq '+'){
							$b+=$c3p;
							$e=$b+$l3p;
						}else{
							$b=$e-$c3p-$l3p;
							$e=$e-$c3p;
						}
					}
                  $pcoord="$1:$b-$e";
                  $strandcol="255,0,0" if($4 eq "+");
                  if($first){
                     $first=0;
                     $bedh1=$bedh;
                     $bedh1 =~ s/TNAME/$options{'T'}.known_miRNAs/;
                     $bedh1 =~ s/TDESC/known miRNAs detected by miRDeep2 for $options{'T'}/;
                     $bedh1 =~ s/BPOS/$pcoord/;
                     print BED "$bedh1";
                  }
                  print BED "$c\t$b\t$e\t$line[0]\t$line[1]\t$s\t$b\t$e\t$strandcol\n";
                }
			}
		}else{
			if($line[1] >= $thres and $line[1] < $maxs){
				if($options{'d'}){
					$line[$seqcol] =~ tr/uU/tT/;
				}
				if($options{'p'}){
					if($line[0] =~ /\|([a-zA-Z0-9_-]*)$/){
						$line[0] = $1;
					}
					print OUT ">$line[0]\n",uc($line[$seqcol]),"\n";

				}else{
					print OUT ">${line[0]}_x${line[5]}_coord:",$coord,"_score:$line[1]\n",uc($line[$seqcol]),"\n";
				}
                chomp $coord;
                $strandcol="0,0,255";
                if($coord =~ /^(\S+):(\d+)\.\.(\d+):(\S)$/){
					my ($c,$b,$e,$s) = ($1,$2,$3,$4);
					if($options{'m'}){
						if($s eq '+'){
							$b+=$c5p;
							$e=$b+$l5p;
						}else{
							$b=$e-$c5p-$l5p;
							$e=$e-$c5p;
						}
					}
					if($options{'k'}){
						if($s eq '+'){
							$b+=$c3p;
							$e=$b+$l3p;
						}else{
							$b=$e-$c3p-$l3p;
							$e=$e-$c3p;
						}
					}
					#$pcoord="$1:$b-$e";
                  $strandcol="255,0,0" if($4 eq "+");
                  if($first){
                     $first=0;
                     $bedh1=$bedh;
                     $bedh1 =~ s/TNAME/$options{'T'}.novel_miRNAs/;
                     $bedh1 =~ s/TDESC/novel miRNAs detected by miRDeep2 for $options{'T'}/;
                     $bedh1 =~ s/BPOS/$pcoord/;
                     print BED "$bedh1";
                  }
                  print BED "$c\t$b\t$e\t$line[0]\t$line[1]\t$s\t$b\t$e\t$strandcol\n";
                }

			}
			next;
		}
	}
	if($not){
		chomp;
		@line=split(/\t/);
		if($options{'d'}){
			$line[$seqcol] =~ tr/uU/tT/;
		}
		if($options{'p'}){

		}else{
			print OUT ">${line[0]}_x${line[4]}\n",uc($line[$seqcol]),"\n";
		}
		next;
	}
}
close OUT;
