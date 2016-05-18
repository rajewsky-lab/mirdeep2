#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;




my $usage=
"$0 file_signature file_structure

This script identifies potential precursors whose structure is basically consistent
with Dicer recognition. The relevant ids are outputted one id per line.
";




#signature file in arf format
my $file_arf=shift or die $usage;

#structure file outputted from RNAfold
my $file_struct=shift or die $usage;

#options
my %options=();
getopts("",\%options);


my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;
my %hash_query;
my %hash_comp;
my %hash_bp;

my $db_old;



#parse structure file outputted from RNAfold
parse_file_struct($file_struct);


#parse signature file in arf format and resolve each potential precursor
parse_file_arf($file_arf);

exit;






sub parse_file_arf{

#    read through the signature blastparsed file, fills up a hash with information on queries
#    (deep sequences) mapping to the current db (potential precursor) and resolve each
#    potential precursor in turn

    my($file)=@_;

    open(FILENAME, $file) or die "Could not open file $file";
    
    while (my $line = <FILENAME>){
	if($line=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){

	    my $query=$1;
	    my $query_lng=$2;
	    my $query_beg=$3;
	    my $query_end=$4;
	    my $query_seq=$5;
	    my $db=$6;
	    my $db_lng=$7;
	    my $db_beg=$8;
	    my $db_end=$9;
	    my $db_seq=$10;
	    my $strand=$11;
	    my $edits=$12;
	    my $edit_string=$13;

	    #only reads that map sense to the potential precursor are considered
	    if($strand eq "-"){next;}
	    
	    #if the new line concerns a new db (potential precursor) then the old db must be resolved
	    if($db_old and $db_old ne $db){
		resolve_potential_precursor();
	    }

	    #resolve the number of reads that the deep sequence represents
	    my $freq=find_freq($query);

	    #read information of the query (deep sequence) into hash
	    $hash_query{$query}{"db_beg"}=$db_beg;
	    $hash_query{$query}{"db_end"}=$db_end;
	    $hash_query{$query}{"strand"}=$strand;
	    $hash_query{$query}{"freq"}=$freq;

	    $db_old=$db;
	}
    }
    resolve_potential_precursor();
}

sub resolve_potential_precursor{
    
#    dissects the potential precursor in parts by filling hashes, and tests if it passes the
#    initial structure filter

    fill_structure();
    
    fill_pri();

    fill_mature();
   
    fill_star();

    fill_loop();

    if(pass_filtering_initial()){

	my $seq=$hash_seq{$db_old};
	my $struct=$hash_struct{$db_old};
	
	print "$db_old\n";
    }
    
    reset_variables();
    
    return;
}



sub pass_filtering_initial{

    #test if the structure forms a plausible hairpin
    unless(pass_filtering_structure()){return 0;}

    return 1;

}



sub pass_filtering_structure{

    #The potential precursor must form a hairpin with miRNA precursor-like characteristics

    #return value
    my $ret=1;

    #potential mature, star, loop and lower flank parts must be identifiable
    unless(test_components()){return 0;}

    #no bifurcations
    unless(no_bifurcations_precursor()){$ret=0;}

    #minimum 14 base pairings in duplex
    unless(bp_duplex()>=0.6){$ret=0;}

    #not more than 6 nt difference between mature and star length
    unless(-6<diff_lng() and diff_lng()<6){$ret=0;}

    return $ret;
}





sub test_components{

    #tests whether potential mature, star, loop and lower flank parts are identifiable

    unless($hash_comp{"mature_struct"}){
	return 0;
    }

    unless($hash_comp{"star_struct"}){
	return 0;
    }

    unless($hash_comp{"loop_struct"}){
   	return 0;
    }

    return 1;
}



sub no_bifurcations_precursor{

    #tests whether there are bifurcations in the hairpin

    #assembles the potential precursor sequence and structure from the expected Dicer products
    #this is the expected biological precursor, in contrast with 'pri_seq' that includes
    #some genomic flanks on both sides

    my $pre_struct;
    my $pre_seq;
    if($hash_comp{"mature_arm"} eq "first"){
	$pre_struct.=$hash_comp{"mature_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"star_struct"};
	$pre_seq.=$hash_comp{"mature_seq"}.$hash_comp{"loop_seq"}.$hash_comp{"star_seq"};
    }else{
	$pre_struct.=$hash_comp{"star_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"mature_struct"};
	$pre_seq.=$hash_comp{"star_seq"}.$hash_comp{"loop_seq"}.$hash_comp{"mature_seq"};
    }

    #read into hash
    $hash_comp{"pre_struct"}=$pre_struct;
    $hash_comp{"pre_seq"}=$pre_seq;

    #simple pattern matching checks for bifurcations
    unless($pre_struct=~/^((\.|\()+..(\.|\))+)$/){
	return 0;
    }

    return 1;
}


sub bp_duplex{

    #fraction of nts of the mature sequence that are base paired in the duplex

    my $lng=$hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1;

    my $duplex_bps=0;
    my $mature_struct=$hash_comp{"mature_struct"};

    #simple pattern matching
    while($mature_struct=~/(\(|\))/g){
	$duplex_bps++;
    }
    return ($duplex_bps/$lng);
}

sub diff_lng{

    #find difference between mature and star lengths

    my $mature_lng=length $hash_comp{"mature_struct"};
    my $star_lng=length $hash_comp{"star_struct"};
    my $diff_lng=$mature_lng-$star_lng;
    return $diff_lng;
}



######################################################

sub fill_structure{

    #reads the dot bracket structure into the 'bp' hash where each key and value are basepaired

    my $struct=$hash_struct{$db_old};
    my $lng=length $struct;

    #local stack for keeping track of basepairings
    my @bps;

    for(my $pos=1;$pos<=$lng;$pos++){
	my $struct_pos=excise_struct($struct,$pos,$pos,"+");

	if($struct_pos eq "("){
	    push(@bps,$pos);
	}

	if($struct_pos eq ")"){
	    my $pos_prev=pop(@bps);
	    $hash_bp{$pos_prev}=$pos;
	    $hash_bp{$pos}=$pos_prev;
	}
    }
    return;
}




sub fill_pri{

    #fills basic specifics on the precursor into the 'comp' hash
    
    my $seq=$hash_seq{$db_old};
    my $struct=$hash_struct{$db_old};
    my $mfe=$hash_mfe{$db_old};
    my $length=length $seq;
    
    $hash_comp{"pri_id"}=$db_old;
    $hash_comp{"pri_seq"}=$seq;
    $hash_comp{"pri_struct"}=$struct;
    $hash_comp{"pri_mfe"}=$mfe;
    $hash_comp{"pri_beg"}=1;
    $hash_comp{"pri_end"}=$length;
    
    return;
}



sub fill_mature{

    #fills specifics on the mature sequence into the 'comp' hash

    my $mature_query=find_mature_query();
    my($mature_beg,$mature_end)=find_positions_query($mature_query);
    my $mature_strand=find_strand_query($mature_query);
    my $mature_seq=excise_seq($hash_comp{"pri_seq"},$mature_beg,$mature_end,$mature_strand);
    my $mature_struct=excise_struct($hash_comp{"pri_struct"},$mature_beg,$mature_end,$mature_strand);
    my $mature_arm=arm_mature($mature_beg,$mature_end,$mature_strand);

    $hash_comp{"mature_query"}=$mature_query;
    $hash_comp{"mature_beg"}=$mature_beg;
    $hash_comp{"mature_end"}=$mature_end;
    $hash_comp{"mature_strand"}=$mature_strand;
    $hash_comp{"mature_struct"}=$mature_struct;
    $hash_comp{"mature_seq"}=$mature_seq;
    $hash_comp{"mature_arm"}=$mature_arm;

    return;
}




sub fill_star{

    #fills specifics on the expected star strand into 'comp' hash ('component' hash)
    
    #if the mature sequence is not plausible, don't look for the star arm
    my $mature_arm=$hash_comp{"mature_arm"};
    unless($mature_arm){$hash_comp{"star_arm"}=0; return;}
 
    #if the star sequence is not plausible, don't fill into the hash
    my($star_beg,$star_end)=find_star();
    my $star_arm=arm_star($star_beg,$star_end);
    unless($star_arm){return;}

    #excise expected star sequence and structure
    my $star_seq=excise_seq($hash_comp{"pri_seq"},$star_beg,$star_end,"+");
    my $star_struct=excise_seq($hash_comp{"pri_struct"},$star_beg,$star_end,"+");

    #fill into hash
    $hash_comp{"star_beg"}=$star_beg;
    $hash_comp{"star_end"}=$star_end;
    $hash_comp{"star_seq"}=$star_seq;
    $hash_comp{"star_struct"}=$star_struct;
    $hash_comp{"star_arm"}=$star_arm;

    return;
}




sub fill_loop{

    #fills specifics on the loop sequence into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the loop
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $loop_beg;
    my $loop_end;

    #defining the begin and end positions of the loop from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' of the loop ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$loop_beg=$hash_comp{"mature_end"}+1;
    }else{
	$loop_end=$hash_comp{"mature_beg"}-1;
    }
    
    if($hash_comp{"star_arm"} eq "first"){
	$loop_beg=$hash_comp{"star_end"}+1;
    }else{
	$loop_end=$hash_comp{"star_beg"}-1;
    }

    #unless the positions are plausible, do not fill into hash
    unless(test_loop($loop_beg,$loop_end)){return;}

    my $loop_seq=excise_seq($hash_comp{"pri_seq"},$loop_beg,$loop_end,"+");
    my $loop_struct=excise_struct($hash_comp{"pri_struct"},$loop_beg,$loop_end,"+");

    $hash_comp{"loop_beg"}=$loop_beg;
    $hash_comp{"loop_end"}=$loop_end;
    $hash_comp{"loop_seq"}=$loop_seq;
    $hash_comp{"loop_struct"}=$loop_struct;

    return;
}




###############################3



sub find_mature_query{

    #finds the query with the highest frequency of reads and returns it
    #is used to determine the positions of the potential mature sequence

    my @queries=sort {$hash_query{$b}{"freq"} <=> $hash_query{$a}{"freq"}} keys %hash_query;
    my $mature_query=$queries[0];
    return $mature_query;
}




sub find_positions_query{

    #subroutine to find the begin and end positions for a given query

    my $query=shift;
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
    return ($beg,$end);
}




sub find_strand_query{

    #subroutine to find the strand for a given query

    my $query=shift;
    my $strand=$hash_query{$query}{"strand"};
    return $strand;
}



sub arm_mature{
 
    #tests whether the mature sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end,$strand)=@_;
 
    #mature and star sequences should alway be on plus strand
    if($strand eq "-"){return 0;}

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,$strand);
    if(defined($struct) and $struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif(defined($struct) and $struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}


sub arm_star{

    #tests whether the star sequence is in the 5' ('first') or 3' ('second') arm of the potential precursor

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    #no overlap between the mature and the star sequence
    if($hash_comp{"mature_arm"} eq "first"){
	($hash_comp{"mature_end"}<$beg) or return 0;
    }elsif($hash_comp{"mature_arm"} eq "second"){
	($end<$hash_comp{"mature_beg"}) or return 0;
    }

    #there should be no bifurcations and minimum one base pairing
    my $struct=excise_seq($hash_comp{"pri_struct"},$beg,$end,"+");
    if($struct=~/^(\(|\.)+$/ and $struct=~/\(/){
	return "first";
    }elsif($struct=~/^(\)|\.)+$/ and $struct=~/\)/){
	return "second";
    }
    return 0;
}

sub find_star{

    #uses the 'bp' hash to find the expected star begin and end positions from the mature positions

    #the -2 is for the overhang
    my $mature_beg=$hash_comp{"mature_beg"};
    my $mature_end=$hash_comp{"mature_end"}-2;
    my $mature_lng=$mature_end-$mature_beg+1;

    #in some cases, the last nucleotide of the mature sequence does not form a base pair,
    #and therefore does not basepair with the first nucleotide of the star sequence.
    #In this case, the algorithm searches for the last nucleotide of the mature sequence
    #to form a base pair. The offset is the number of nucleotides searched through.
    my $offset_star_beg=0;
    my $offset_beg=0;

    #the offset should not be longer than the length of the mature sequence, then it
    #means that the mature sequence does not form any base pairs
    while(not $offset_star_beg and $offset_beg<$mature_lng){
	if($hash_bp{$mature_end-$offset_beg}){
	    $offset_star_beg=$hash_bp{$mature_end-$offset_beg};
	}else{
	    $offset_beg++;
	}
    }
    #when defining the beginning of the star sequence, compensate for the offset
    my $star_beg=$offset_star_beg-$offset_beg;

    #same as above
    my $offset_star_end=0;
    my $offset_end=0;
    while(not $offset_star_end and $offset_end<$mature_lng){
	if($hash_bp{$mature_beg+$offset_end}){
	    $offset_star_end=$hash_bp{$mature_beg+$offset_end};
	}else{
	    $offset_end++;
	}
    }
    #the +2 is for the overhang
    my $star_end=$offset_star_end+$offset_end+2;

    return($star_beg,$star_end);
}



sub test_loop{

    #tests the loop positions

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}





###############################################################

  

sub reset_variables{

    #resets the hashes for the next potential precursor

    %hash_query=();
    %hash_comp=();
    %hash_bp=();

    return;
}




sub parse_file_struct{
 
    #parses the output from RNAfold and reads it into hashes

    my($file) = @_;
    my($id,$desc,$seq,$struct,$mfe) = ();

    open (FILE_STRUCT, "<$file") or die "can not open $file\n";
    while (<FILE_STRUCT>)
    {
        chomp;
        if (/^>(\S+)\s*(.*)/)
	{
	    $id          = $1;
	    $desc        = $2;
	    $seq         = "";
	    $struct      = "";
	    $mfe         = "";
	    while (<FILE_STRUCT>){
                chomp;
                if (/^>(\S+)\s*(.*)/){
		    $hash_desc{$id}   = $desc;
		    $hash_seq{$id}    = $seq;
		    $hash_struct{$id} = $struct;
		    $hash_mfe{$id}    = $mfe;

		    $id          = $1;
		    $desc        = $2;
		    $seq         = "";
		    $struct      = "";
		    $mfe         = "";

		    next;
                }
		if(/^\w/){
		    tr/uU/tT/;
		    $seq .= $_;
		}if(/((\.|\(|\))+)/){
		    $struct .=$1;
		}
		if(/\((\s*-\d+\.\d+)\)/){
		    $mfe = $1;
		}
	    
	    }
        }
    }

    $hash_desc{$id}        = $desc;
    $hash_seq{$id}         = $seq;
    $hash_struct{$id}      = $struct;
    $hash_mfe{$id}         = $mfe;

    close FILE_STRUCT;
    return;
}


sub revcom{
    
    #reverse complement

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}

sub rev{

    #reverses the order of nucleotides in a sequence

    my($sequence)=@_;

    my $rev=reverse $sequence;   

    return $rev;
}

sub com{

    #the complementary of a sequence

    my($sequence)=@_;

    $sequence=~tr/acgtuACGTU/TGCAATGCAA/;   
 
    return $sequence;
}
  
sub find_freq{

    #finds the frequency of a given read query from its id.

    my($query)=@_;

    if($query=~/x(\d+)/){
	my $freq=$1;
	return $freq;
    }else{
	return 0;
    }
}


sub excise_seq{

    #excise sub sequence from the potential precursor

    my($seq,$beg,$end,$strand)=@_;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $db_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($seq)){return 0;}
 
    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_seq=substr($seq,$beg-1,$end-$beg+1);

    #if on the minus strand, the reverse complement should be returned
    if($strand eq "-"){
	$sub_seq=revcom($sub_seq);
    }

    return $sub_seq;

}



sub excise_struct{

    #excise sub structure

    my($struct,$beg,$end,$strand)=@_;
    my $lng=length $struct;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $db_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($struct)){return 0;}

    #the blast parsed format is 1-indexed, substr is 0-indexed
    my $sub_struct=substr($struct,$beg-1,$end-$beg+1);
 
    return $sub_struct;
}



