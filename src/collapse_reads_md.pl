#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

my $usage =
"$0 file_fasta prefix

Collapses reads in the fasta file to make each sequence entry unique. Each collapsed
entry will have an identifier that follows an underscore '_' separated format. Example:

>mmu_1189273_x10

The first field 'mmu' shows which sample the sequence is from. This prefix is given on
the command line, and must consist of exactly three alphabetic letters. The second field
'118273' is a running number that secures that each identifier is unique. The third
field '_x10' shows how many times the given sequence was present in the dataset.

-a    outputs progress

example use:
collapse_reads.pl reads.fa mmu
";

my $file_fasta=shift or die $usage;
my $prefix=shift or die $usage;


my %options=();
getopts("a",\%options);

my %hash;

test_prefix($prefix);

parse_file_fasta_seqkey(\$file_fasta,\%hash);

print_hash_seqkey(\%hash);



sub parse_file_fasta_seqkey{
    
    my($file,$hash)=@_;
    my($id,$seq)=();
 
    if($options{a}){print STDERR "reading file into hash\n";}
    my $running_1=0;
   
    open (FASTA, "<$$file") or die "can not open $$file\n";
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
				
				my $cnt=find_cnt($id);
				$seq=~tr/[acgtun\.]/[ACGTTNN]/;
				$$hash{$seq}+=$cnt;
				$running_1++;
				if($options{a}){print STDERR "$running_1\r";}
				
				$id   = $1;
				$seq  = "";
				next;
			}
			$seq.= $_;
	    }
	}
    }
    
    my $cnt=find_cnt($id);
    $seq=~tr/[acgtun\.]/[ACGTTNN]/;
    $$hash{$seq}+=$cnt;
    $running_1++;
    if($options{a}){print STDERR "$running_1\r";}
    
    close FASTA;
}

sub print_hash_seqkey{
    
    my ($hash)=@_;
    if($options{a}){print STDERR "sorting hash\n";}
    my $running_2=0;
    if($options{a}){print STDERR "printing hash\n";}
    foreach my $key(sort {$$hash{$b} <=> $$hash{$a}} keys %$hash){
	
		my $cnt=$$hash{$key};
		
		print ">$prefix\_$running_2\_x$cnt\n$key\n";
		$running_2+=$cnt;
		if($options{a}){print STDERR "$running_2\r";}
    }
}




sub find_cnt{
    
    #finds the frequency of a given read query from its id.
    
    my($query)=@_;
    
    if($query=~/_x(\d+)/){
		
		my $cnt=$1;
		
		return $cnt;
		
    }else{
	
		return 1;
    }
}


sub test_prefix{

    my $prefix=shift;
	
    unless($prefix=~/^\w\w\w$/ and !($prefix=~/_/)){
		
		die "prefix $prefix does not contain exactly three alphabet letters\n";
    }
    return;
}
