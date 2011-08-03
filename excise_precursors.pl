#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;



my $usage =
"$0 file_fasta file_arf precursor.coords

This script excises potential miRNA precursor sequences from a genome.
The fasta file should be the relevant genome, and the arf file should
contain the read mappings.
The precursors.coords designates where to write the precursor genome coordinates to.

-a integer   Only excise if the potential mature microRNA is represented
             by a number of reads equal to or higher than the integer
             (default 2).
-b           Output progress to screen
";


my $file_fasta=shift or die $usage;
my $file_arf=shift or die $usage;

my $coord_file=shift or die $usage;

open PF,">$coord_file";

my %options=();
getopts("a:b",\%options);

my %hash_db_lng;
my %hash_pos;

my $count_lines=0;
my $count_excisions=0;

my $freq_min=2;
if($options{a}){$freq_min=$options{a}};

if($options{b}){print STDERR "finding lengths of genome contigs\n";}
#parse_genome_lengths(\$file_fasta,\%hash_db_lng);

parse_file_arf($file_arf);

if($options{b}){print STDERR "reading the genome into memory and excising potential precursors\n";}
parse_genome_and_excise($file_fasta);

if($options{b}){print STDERR "potential precursors excised\n";}

close PF;

exit;




sub parse_file_arf{

    my($file)=@_;

    my $lines=`cat $file | wc -l`;
    chomp $lines;

    if($options{b}){print STDERR "reading the mapping file into memory, total lines=$lines\n";}

    open(FILENAME, $file) or die "Could not open file $file";
    
    while (my $line = <FILENAME>){
	if($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){

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

# 	    #exit if the db is not contained in fasta file
# 	    unless($hash_db_lng{$db}){
# 		die "$db not present in the genome fasta file\n";
# 	    }

	    my $freq=find_freq($query);

	    #reverse positions on the minus strand - allows for faster string ops later
#	    if($strand eq "-"){
#		my $db_lng=$hash_db_lng{$db};
#		($db_beg,$db_end)=reverse_positions($db_lng,$db_beg,$db_end);
#	    }

	    #read into position hash
	    insertfeature($db,$strand,$db_beg,$db_end,$freq);

	    $count_lines++;

	    if($options{b}){
#		print STDERR "$count_lines\r";
	    }
	}
    }
    close FILENAME;
}




sub insertfeature{

    my($db,$strand,$db_beg,$db_end,$freq)=@_;
 
    $hash_pos{$db}{$strand}{$db_beg}{$db_end}+=$freq;
 
}




sub parse_genome_lengths{
    my ($file,$hash_db_lng) = @_;
    my ($id, $desc, $sequence) = ();

    open (FASTA, "<$$file") or die "can not open $$file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $$hash_db_lng{$id} = length $sequence;
		    $id         = $1;
		    $desc       = $2;
		    $sequence   = "";
		    next;
                }
                $sequence .= $_;
            }
        }
    }
    $$hash_db_lng{$id} = length $sequence;
    close FASTA;
}




sub parse_genome_and_excise{
    my ($file) = @_;
    my ($id, $desc, $sequence) = ();

    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
        {
            $id       = $1;
            $desc     = $2;
            $sequence = "";
            while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
                    excise(\$id,\$sequence);
                    $id         = $1;
                    $desc       = $2;
                    $sequence   = "";
                    next;
                }
                $sequence .= $_;
            }
        }
    }
    excise(\$id,\$sequence);
    close FASTA;
    return;
}





sub excise{

    my ($db,$db_seq)=@_;
    
    my @strands=sort keys %{$hash_pos{$$db}};
    foreach my $strand(@strands){

#	my $db_seq=();
	#use the genome contig or the reverse complement?
#	if($strand eq "+"){
	    #$db_seq=$sequence;
#	}else{
	    #$db_seq=revcom($db_seq);
#	}
	
	#total length of chromosome or genome contig
	my $db_lng=length($$db_seq);
	
	#the most 3' position that has yet been excised
	my $db_limit=0; 
	
	#scans from 5' to 3' end of the chromosome
	my @db_begs=sort {$a<=>$b} keys %{$hash_pos{$$db}{$strand}};
	foreach my $db_beg(@db_begs){
	    
	    #reads with the same 5' position can have different 3' positions
	    my @db_ends=sort {$a<=>$b} keys %{$hash_pos{$$db}{$strand}{$db_beg}};
	    foreach my $db_end(@db_ends){
		
		#height of read stack with those exact 5' and 3' positions
		my $freq=$hash_pos{$$db}{$strand}{$db_beg}{$db_end};
		
		#what is the highest read stack downstream (up to 70 nt downstream)?
		my $freq_max_ds=find_freq_max_downstream($$db,$strand,$db_beg,$db_end);
		
		#if read stack to low, if higher read stack downstream or if this locus has already
		#been excised, then continue to next stack
		if($freq<$freq_min or $freq<$freq_max_ds or $db_beg<$db_limit){next;}
		
		#else excise to sequences, corresponding to the read stack being the mature sequence
		#in the 5' arm or the 3' arm of the precursor hairpin
		my $excise_beg=$db_beg-70;
		my $excise_end=$db_end+20;
		
		#print out in fasta format
		print_positions($db,\$strand,$db_seq,\$db_lng,\$excise_beg,\$excise_end);
		
		$excise_beg=$db_beg-20;
		$excise_end=$db_end+70;
		
		print_positions($db,\$strand,$db_seq,\$db_lng,\$excise_beg,\$excise_end);
		
		#the most 3' position that has yet been excised
		$db_limit=$excise_end;
		
	    }
	}
    }
    
    return;
}



sub print_positions{

    my($db,$strand,$db_seq,$db_lng,$excise_beg,$excise_end)=@_;

    my $excise_seq=excise_position($db_seq,$db_lng,$excise_beg,$excise_end);
	if($$strand eq "-"){$excise_seq = revcom($excise_seq);}

    print ">$$db\_$count_excisions\n$excise_seq\n";
	print PF ">$$db\_$count_excisions\t$$strand\t$$excise_beg\t$$excise_end\n";

    $count_excisions++;

    return;
}	      





sub find_freq_max_downstream{
    
    my ($db,$strand,$db_beg,$db_end)=@_;

    my $freq_max=0;
    
    for(my $pos_beg=$db_beg+1; $pos_beg<=$db_end+70; $pos_beg++){

	if(defined $hash_pos{$db}{$strand}{$pos_beg}){
	
	    my @pos_ends=sort {$a<=>$b} keys %{$hash_pos{$db}{$strand}{$pos_beg}};
	    
	    foreach my $pos_end(@pos_ends){
		
		my $freq=$hash_pos{$db}{$strand}{$pos_beg}{$pos_end};

		if($freq>$freq_max){
		    $freq_max=$freq;
		}
	    }
	}
    }
    return $freq_max;
}




sub excise_position{

    my($db_seq,$db_lng,$excise_beg,$excise_end)=@_;

    my $excise_beg_limit=max2(1,$$excise_beg);
    my $excise_end_limit=min2($$db_lng,$$excise_end);

    my $excise_lng=$excise_end_limit-$excise_beg_limit+1;
    my $excise_seq=substr($$db_seq,$excise_beg_limit-1,$excise_lng);
 
    return $excise_seq;
}



sub find_strand{

    #A subroutine to find the strand, parsing different blast formats
    my($other)=@_;

    my $strand="+";

    if($other=~/-/){
	$strand="-";
    }

    if($other=~/minus/i){
	$strand="-";
    }

    return($strand);
}


sub reverse_positions{

    #A subroutine to find positions relative to the minus strand
    my($length,$begin,$end)=@_;

    my $new_end=$length-$begin+1;
    my $new_beg=$length-$end+1;

    return($new_beg,$new_end);
}



sub rev{

    my($sequence)=@_;

    my $rev=reverse $sequence;   

    return $rev;
}

sub com{

    my($sequence)=@_;

    $sequence=~tr/acgtuACGTU/TGCAATGCAA/;   
 
    return $sequence;
}

sub revcom{

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}

    
sub find_freq{
    #finds the frequency of a given read query from its id.

    my($query)=@_;

    if($query=~/_x(\d+)/){
	my $freq=$1;
	return $freq;
    }else{
	print STDERR "Problem with read format\n";
	return 1;
    }
}


sub max2 {
        my($a, $b) = @_;
        return ($a>$b ? $a : $b);
}

sub min2  {
        my($a, $b) = @_;
        return ($a<$b ? $a : $b);
}
