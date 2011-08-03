#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;

my $usage =
"$0 file_arf

Parses an arf file, and discards select mappings. Any files inputted with
the options should be single-column format, with a single id on each line.

-a int     Discard mappings of edit distance higher than this
-b int     Discard mappings of read queries shorter than this
-c int     Discard mappings of read queries longer than this
-d file    Discard read queries not in this file
-e file    Discard read queries in this file
-f file    Discard reference dbs not in this file
-g file    Discard reference dbs in this file
-h         Discard remaining suboptimal mappings
-i int     Discard remaining suboptimal mappings and discard any
           reads that have more remaining mappings than this
-j         Remove any unmatched nts in the very 3' end
-k         Output progress to standard output
";

#input
my $file_arf=shift or die $usage;

#options
my %options=();
getopts("a:b:c:d:e:f:g:hi:jk",\%options);

#global hashes
my %hash_fasta;
my %hash_edits;
my %hash_edit_best;
my %hash_nr_mappings;

my %hash_queries_incl;
my %hash_queries_excl;
my %hash_dbs_incl;
my %hash_dbs_excl;

#running index		
my $running=0;

#scan mode - a flag variable
my $scan=0;


#MAIN

#if($options{k}){print STDERR "parsing mappings\n";}

if($options{d}){parse_file_ids(\$options{d},\%hash_queries_incl);}
if($options{e}){parse_file_ids(\$options{e},\%hash_queries_excl);}
if($options{f}){parse_file_ids(\$options{f},\%hash_dbs_incl);}
if($options{g}){parse_file_ids(\$options{g},\%hash_dbs_excl);}

if($options{h} or $options{i}){scan($file_arf);}

parse($file_arf);

exit;



sub parse_file_arf{
    
    my ($file_arf)=@_;
  
    open (FILE_ARF, "<$file_arf") or die "can not open $file_arf\n";
    while (my $line=<FILE_ARF>){
	
	if($line=~/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
	    
	    my $query=$1;
	    my $query_map_lng=$2;
	    my $query_beg=$3;
	    my $query_end=$4;
	    my $query_seq=$5;
	    my $db=$6;
	    my $db_map_lng=$7;
	    my $db_beg=$8;
	    my $db_end=$9;
	    my $db_seq=$10;
	    my $strand=$11;
	    my $edits=$12;
	    my $edit_string=$13;

	    $running++;

#	    if($options{k}){print STDERR "$running\r";}

	    if($options{j}){($query_map_lng,$query_end,$query_seq,$db_map_lng,$db_end,$db_seq,$edits,$edit_string)=remove_trailing_nts($query_map_lng,$query_end,$query_seq,$db_map_lng,$db_end,$db_seq,$edits,$edit_string);}

	    #discard mappings based on criteria a-g
	    if(defined $options{a} and $options{a}<$edits){next;}
	    if($options{b} and $query_map_lng<$options{b}){next;}
	    if($options{c} and $options{c}<$query_map_lng){next;}
	    if($options{d} and not(defined $hash_queries_incl{$query})){next;}
	    if($options{e} and (defined $hash_queries_excl{$query})){next;}
	    if($options{f} and not(defined $hash_dbs_incl{$db})){next;}
	    if($options{g} and (defined $hash_dbs_excl{$db})){next;}

	    #unless suboptimal mappings should be disregarded, there is no need to read into hash
	    unless($options{h} or $options{i}){
		print "$query\t$query_map_lng\t$query_beg\t$query_end\t$query_seq\t$db\t$db_map_lng\t$db_beg\t$db_end\t$db_seq\t$strand\t$edits\t$edit_string\n";
		next;
	    }

	    #in scanning mode, fill the edit hash
	    if($scan){
		$hash_edits{$query}{$edits}++;
	    #in parsing mode, evaluate if the line should be printed	
	    }else{
		my $evaluation=evaluate_query($query,$edits);
		if($evaluation){print "$query\t$query_map_lng\t$query_beg\t$query_end\t$query_seq\t$db\t$db_map_lng\t$db_beg\t$db_end\t$db_seq\t$strand\t$edits\t$edit_string\n";}
	    }
	}
    }
}


sub remove_trailing_nts{

    my($query_map_lng,$query_end,$query_seq,$db_map_lng,$db_end,$db_seq,$edits,$edit_string)=@_;

    while($edit_string=~/M$/){

	$query_map_lng--;
	$query_end--;
	chop $query_seq;
	$db_map_lng--;
	$db_end--;
	chop $db_seq;
	$edits--;
	chop $edit_string;
    }

    return ($query_map_lng,$query_end,$query_seq,$db_map_lng,$db_end,$db_seq,$edits,$edit_string);
}



sub fill_hash{

    my @queries=sort keys %hash_edits;
    foreach my $query(@queries){

	#find the best mapping (lowest edit distance) for the query read
	my @edits=sort {$a<=>$b} keys %{$hash_edits{$query}};
	my $edit_best=$edits[0];

	#find the number of mappings of this edit distance
	$hash_edit_best{$query}=$edit_best;
	$hash_nr_mappings{$query}=$hash_edits{$query}{$edit_best};
    }
    return;
}    




sub evaluate_query{

    my ($query,$edits)=@_;

    #mapping should not be suboptimal
    if($hash_edit_best{$query}<$edits){return 0;}

    #read should not map more times than designated
    if($options{i} and $options{i}<$hash_nr_mappings{$query}){return 0;}

    return 1;
}


sub parse_file_ids{

    my ($file,$hash) = @_;

    #read id file into hash
    if($options{k}){print STDERR "reading id file into memory\n";}
    open (FILE, "<$$file") or die "can not open $$file\n";
    while (my $line=<FILE>){

	if($line=~/^(\S+)/){
	    my $id=$1;
	    $$hash{$id}=1;
	}
    }
}


sub scan{

    my $file=shift;
	if($options{k}){
		my $lines=`cat $file | wc -l`;
		chomp $lines;
		print STDERR "scanning mappings, total=$lines\n";
	}

    
	$scan=1;
    parse_file_arf($file);
 
    $scan=0;
    if($options{k}){print STDERR "resolving best mappings for each read\n";}
    fill_hash();

    $running=0;

    return;
}


sub parse{

    if($options{k}){print STDERR "parsing and printing mappings\n";}
    parse_file_arf($file_arf);
    if($options{k}){print STDERR "mappings printed\n\n";}

}
