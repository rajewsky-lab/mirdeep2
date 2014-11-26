#!/usr/bin/perl

#######################################################################
#
#
# This script extracts the mirdeep2 precursors from the result.csv file 
# and writes a fasta file and a bed file for these
# See 'Usage' for more information.
# 
# Bed files are only created for the precursor sequences. 
#
# S.M. Jan. 2012
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
        open BED,">$od/not_$names[$seqcol]${timestamp}_score${score}_to_$max.bed";
		next;
	}else{}
	next if(/^\s*$/);
	if($novel or $known){
		chomp;
		@line=split(/\t/);
		my $coord='na';
		$coord=$line[16] if($line[16]);
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
                  $pcoord="$1:$2-$3";
                  $strandcol="255,0,0" if($4 eq "+");
                  if($first){ 
                     $first=0;
                     $bedh1=$bedh;
                     $bedh1 =~ s/TNAME/$options{'T'}.known_miRNAs/;
                     $bedh1 =~ s/TDESC/known miRNAs detected by miRDeep2 for $options{'T'}/;
                     $bedh1 =~ s/BPOS/$pcoord/;
                     print BED "$bedh1";
                  }
                  print BED "$1\t$2\t$3\t$line[0]\t$line[1]\t$4\t$2\t$3\t$strandcol\n";
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
                  $strandcol="255,0,0" if($4 eq "+");
                  if($first){ 
                     $first=0;
                     $bedh1=$bedh;
                     $bedh1 =~ s/TNAME/$options{'T'}.novel_miRNAs/;
                     $bedh1 =~ s/TDESC/novel miRNAs detected by miRDeep2 for $options{'T'}/;
                     $bedh1 =~ s/BPOS/$pcoord/;
                     print BED "$bedh1";
                  }
                  print BED "$1\t$2\t$3\t$line[0]\t$line[1]\t$4\t$2\t$3\t$strandcol\n";
                }

			}
			next;
		}
	}
	if($not){
		chomp;
		@line=split(/\t/);
		my $coord='na';
		$coord=$line[16] if($line[16]);
		if($options{'d'}){
			$line[$seqcol] =~ tr/uU/tT/;
		}
		if($options{'p'}){
			
		}else{
			print OUT ">${line[0]}_x${line[4]}\n",uc($line[$seqcol]),"\n";
		}
        chomp $coord;
        $strandcol="0,0,255";
        if($coord =~ /^(\S+):(\d+)\.\.(\d+):(\S)$/){
          $strandcol="255,0,0" if($4 eq "+");
          if($first){ 
                     $first=0;
                     $bedh1=$bedh;
                     $bedh1 =~ s/TNAME/$options{'T'}.not_detected_miRNAs/;
                     $bedh1 =~ s/TDESC/miRNAs not detected by miRDeep2 for $options{'T'}/;
                     $bedh1 =~ s/BPOS/$pcoord/;
                     print BED "$bedh1";
                  }
          print BED "$1\t$2\t$3\t$line[0]\t$line[1]\t$4\t$2\t$3\t$strandcol\n";
        }



		next;
	}
}
close OUT;
