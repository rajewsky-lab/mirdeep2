#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;
use File::Copy;
use File::Path;


my $usage =
"$0 file_reads file_precursors read_align_edit_distance

This script prepares the signature file for miRDeep. Options:

-a file  Fasta file with the sequences of known mature miRNAs for the species. 
         These sequences will not influence the miRDeep scoring, but will 
         subsequently make it easy to estimate sensitivity of the run.
-b       Output progress to screen

";


#input
my $file_reads=shift or die $usage;
my $file_precursors=shift or die $usage;
my $read_align_edit_distance = shift or die $usage;


#options
my %options=();
getopts("a:bo:",\%options);

#global variables
my $ltime=time();
my $dir="dir_prepare_signature$ltime";
my $outfile=$options{'o'} or die "no outfile specified with options{'o'}\n";
#MAIN

if($options{b}){print STDERR "preparing signature file\n";}
mkdir $dir;
copy($file_precursors,$dir);

if($options{b}){print STDERR "constructing index of precursors\n";}
system("bowtie-build $file_precursors $dir/precursors.ebwt > /dev/null");

if($options{b}){print STDERR "mapping reads to precursors\n";}
system("bowtie -f -v $read_align_edit_distance -a --best --strata --norc $dir/precursors.ebwt $file_reads $dir/reads_vs_precursors.bwt 2> /dev/null");

system("convert_bowtie_output.pl $dir/reads_vs_precursors.bwt > $dir/reads_vs_precursors.arf");

if($options{a}){

    my $file_mature=$options{a};

    if($options{b}){print STDERR "mapping reference mature miRNAs to precursors\n";}
    system("bowtie -f -v 0 -a --best --strata --norc $dir/precursors.ebwt $options{a} $dir/mature_vs_precursors.bwt 2> /dev/null");

    system("convert_bowtie_output.pl $dir/mature_vs_precursors.bwt > $dir/mature_vs_precursors.arf");

    if($options{b}){print STDERR "sorting rows\n";}
    cat_files("$dir/reads_vs_precursors.arf","$dir/mature_vs_precursors.arf","$dir/signature_unsorted.arf");
    
#    system("cat $dir/reads_vs_precursors.arf $dir/mature_vs_precursors.arf > $dir/signature_unsorted.arf");
    #Sortarf("$dir/signature_unsorted.arf");
	presort("$dir/signature_unsorted.arf");
	system("sort -nk1 $dir/signature_unsorted.arf.tmp > $dir/signature_unsorted.arf.tmp2");
	system("cut -f2-14 $dir/signature_unsorted.arf.tmp2 > $outfile");
#	system("sort -V -k6 $dir/signature_unsorted.arf > $outfile");
	

}else{

    if($options{b}){print STDERR "sorting rows\n";}
	#Sortarf("$dir/reads_vs_precursors.arf");
	presort("$dir/reads_vs_precursors.arf");
	system("sort -nk1 $dir/reads_vs_precursors.arf.tmp > $dir/reads_vs_precursors.arf.tmp2");
	system("cut -f2-14 $dir/reads_vs_precursors.arf.tmp2 > $outfile");
#	system("sort -V -k6 $dir/reads_vs_precursors.arf > $outfile");
}


## remove temporary directory
#rmtree($dir);

if($options{b}){print STDERR "signature file prepared\n\n";}



############################################################
############################################################
##                                                        ##
## Subroutines                                            ##
##                                                        ##
############################################################
############################################################

sub presort{
	my $file=shift;
	open IK,"$file" or die "no arf file given\n";
	open IKT,">$file.tmp" or die "tmp file could not be opened\n";
	
	my %index=();
	my $count=0;
	my @l;

	while(<IK>){
		@l=split();
		if(not $index{$l[5]}){
			$count++;
			$index{$l[5]}=$count;
		}
		print IKT "$index{$l[5]}\t$_";
	}
	close IK;
	close IKT;
}


sub Sortarf {
	my $file=shift;
	open IN,"<$file" or die "FILE $file not found $!\n";

	my %hash;
	my %hash2;
	my $counter=0;
	my $line;
	while($line = <IN>){
		chomp $line;
#	                 readid  len     start   end   readseq   genid   glen  gstart    gend   gseq    strand  mm      nt-info
		if($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)$/){
			$hash{$counter} = $6;
			$hash2{$counter} = $line;
			$counter++;
			
		}else{
			print "wrong line format\n $line\n";
		}
	}
	close IN;

	## sorting rows by 6th column 
	for(sort {$hash{$a} cmp $hash{$b} } keys %hash){
		print "$hash2{$_}\n";                         ## print hash
	}
}

exit;



sub cat_files{
    
    my($file_1,$file_2,$file_out)=@_;


    open OUT, ">$file_out" or die "cannot print to $file_out\n";

    open IN_1, "<$file_1" or die "cannot read from $file_1\n";
    
    while(my $line = <IN_1>){

	print OUT "$line";
    }

    close IN_1;

    open IN_2, "<$file_2" or die "cannot read from $file_2\n";
    
    while(my $line = <IN_2>){

	print OUT "$line";
    }

    close IN_2;

    close OUT;

    return;
}


