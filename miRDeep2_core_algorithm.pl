#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;



################################# MIRDEEP #################################################

################################## USAGE ##################################################


my $usage=
"$0 file_signature file_structure

This is the core algorithm of miRDeep. It takes as input a file in blastparsed format with
information on the positions of reads aligned to potential precursor sequences (signature).
It also takes as input an RNAfold output file, giving information on the sequence, structure
and mimimum free energy of the potential precursor sequences.

Extra arguments can be given. -s specifies a fastafile containing the known mature miRNA
sequences that should be considered for conservation purposes. -t prints out the potential
precursor sequences that do _not_ exceed the cut-off (default prints out the sequences that
exceeds the cut-off). -u gives limited output, that is only the ids of the potential precursors
that exceed the cut-off. -v varies the cut-off. -x is a sensitive option for Sanger sequences
obtained through conventional cloning. -z consider the number of base pairings in the lower
stems (this option is not well tested).

-h print this usage
-s fasta file with known miRNAs
-t print filtered
-u limited output (only ids)
-v cut-off (default 1)
-x sensitive option for Sanger sequences
-y file with randfold p-values
-z consider Drosha processing
";





############################################################################################

################################### INPUT ##################################################


#signature file in arf format
my $file_arf=shift or die $usage;

#structure file outputted from RNAfold
my $file_struct=shift or die $usage;

#options
my %options=();
getopts("hs:tuv:xy:zl:",\%options);






#############################################################################################

############################# GLOBAL VARIABLES ##############################################


#parameters
my $seed_lng=7;
my $mature_lng_max=25;

my $score_star=3.9;
my $score_star_not=-1.3;
my $score_seed=3;
my $score_seed_not=-0.6;
my $score_randfold=1.6;
my $score_randfold_not=-2.2;
my @scores_stem=(-3.1,-2.3,-2.2,-1.6,-1.5,0.1,0.6,0.8,0.9,0.9,0);
my $score_min=1;
if($options{v}){$score_min=$options{v};}
if($options{x}){$score_min=-5;}

my $e=2.718281828;

#hashes
my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;
my %hash_seeds;
my %hash_mirs;
my %hash_randfold;
my %hash_query;
my %hash_comp;
my %hash_bp;
my %hash_stars;


#other variables
my $db_old;
my $lines;
my $out_of_bound;




##############################################################################################

################################  MAIN  ######################################################


#print help if that option is used
if($options{h}){die $usage;}

#parse structure file outputted from RNAfold
parse_file_struct($file_struct);

#if conservation is scored, the fasta file of known miRNA sequences is parsed
if($options{s}){create_hash_seeds($options{s})};

#if randfold p-value is scored, a randfold output file of the precursors is parsed
if($options{y}){parse_file_randfold($options{y})};

#parse signature file in arf format and resolve each potential precursor
parse_file_arf($file_arf);

exit;




##############################################################################################

############################## SUBROUTINES ###################################################



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

	    #save the signature information
	    $lines.=$line;

	    $db_old=$db;
	}
    }
    resolve_potential_precursor();
}

sub resolve_potential_precursor{

#    dissects the potential precursor in parts by filling hashes, and tests if it passes the
#    initial filter and the scoring filter

#    binary variable whether the potential precursor is still viable
    my $ret=1;

    fill_structure();

    fill_pri();

    fill_mature();

    fill_star();

    fill_loop();

    fill_lower_flanks();

#    do_test_assemble();

#    this is the actual classification
    unless(pass_filtering_initial() and pass_threshold_score()){$ret=0;}

    print_results($ret);

    reset_variables();

    return;

}



sub print_results{

    my $ret=shift;

#    print out if the precursor is accepted and accepted precursors should be printed out
#    or if the potential precursor is discarded and discarded potential precursors should
#    be printed out

    if((not $options{t} and $ret) or ($options{t} and not $ret)){
	#full output
	unless($options{u}){
	    print_hash_comp();
	    #print $lines."\n\n";
	    return;
	}
	#limited output (only ids)
	my $id=$hash_comp{"pri_id"};
	print "$id\n";
    }
}







sub pass_threshold_score{

#    this is the scoring

    #minimum free energy of the potential precursor
    my $score_mfe=score_mfe($hash_comp{"pri_mfe"});

    #count of reads that map in accordance with Dicer processing
    my $score_freq=score_freq($hash_comp{"freq_total"});

    #basic score
    my $score=$score_mfe+$score_freq;

    #scoring of conserved seed/seed (optional)
    if($options{s}){

	#if the seed is conserved
	if(test_seed_conservation()){

	    #seed from position 2-8
	    my $seed=substr($hash_comp{"mature_seq"},1,$seed_lng);

	    #resolve DNA/RNA ambiguities
	    $seed=~tr/[T]/[U]/;

	    #print score contribution
	    $hash_comp{"score_seed"}=$score_seed;

	    #add to score
	    $score+=$score_seed;

	#if the seed is not conserved
	}else{
	    #print (negative) score contribution
	    $hash_comp{"score_seed"}=$score_seed_not;

	    #add (negative) score contribution
	    $score+=$score_seed_not;
	}
    }

    #if the majority of potential star reads fall as expected from Dicer processing
    if($hash_comp{"star_read"}){
	$hash_comp{"score_star"}=$score_star;
	$score+=$score_star;
    }else{
	$hash_comp{"score_star"}=$score_star_not;
	$score+=$score_star_not;
    }

    #score lower stems for potential for Drosha recognition (highly optional)
    if($options{z}){
	my $stem_bp=$hash_comp{"stem_bp"};
	my $score_stem=$scores_stem[$stem_bp];
	$score+=$score_stem;
	$hash_comp{"score_stem"}=$score_stem;
    }


    #score for randfold (optional)
    if($options{y}){

	#randfold p-value not defined for this potential precursor
	if(not(defined $hash_randfold{$hash_comp{"pri_id"}})){$hash_comp{"score_randfold"}="0";}

	#randfold value<0.05
	elsif(test_randfold()){$score+=$score_randfold;$hash_comp{"score_randfold"}=$score_randfold;}

	#randfold value>0.05
	else{$score+=$score_randfold_not;$hash_comp{"score_randfold"}=$score_randfold_not;}
    }

    #round off values to one decimal
    my $round_mfe=round($score_mfe*10)/10;
    my $round_freq=round($score_freq*10)/10;
    my $round=round($score*10)/10;

    #print scores
    $hash_comp{"score_mfe"}=$round_mfe;
    $hash_comp{"score_freq"}=$round_freq;
    $hash_comp{"score"}=$round;

    #return 1 if the potential precursor is accepted, return 0 if discarded
    unless($score>=$score_min){return 0;}
    return 1;
}

sub test_randfold{

    my $pri_id=$hash_comp{"pri_id"};

    my $p_value=$hash_randfold{$pri_id};

    if($p_value<=0.05){return 1;}

    return 0;
}


sub print_file{

    #print string to file

    my($file,$string)=@_;

    open(FILE, ">$file");
    print FILE "$string";
    close FILE;
}


sub test_seed_conservation{

    #test if seed is identical to seed from known miRNA, return 1 or 0

    my $seed=substr($hash_comp{"mature_seq"},1,$seed_lng);
    $seed=~tr/[T]/[U]/;
    if($hash_seeds{$seed}){

	my $seed_family=$hash_mirs{$seed};
	if($seed_family=~/^(\S+)/){
	    #example of a miRNA with the same seed
	    $hash_comp{"seed_family"}=$1;
	}

	return 1;
    }

    return 0;
}



sub pass_filtering_initial{

    #test if the structure forms a plausible hairpin
    unless(pass_filtering_structure()){$hash_comp{"problem_structure"}="The candidate has been discarded because the structure is inconsistent with Dicer processing:\n"; return 0;}

    #test if >90% of reads map to the hairpin in consistence with Dicer processing
    unless(pass_filtering_signature()){$hash_comp{"problem_signature"}="The candidate has been discarded because the signature is inconsistent with Dicer processing:\n";return 0;}

    return 1;

}


sub pass_filtering_signature{

    #if the putative mature sequence is longer than the designated number of nts, discard
    my $mature_lng=$hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1;
    if($mature_lng>$mature_lng_max){$hash_comp{"problem_mature_lng"}="the candidate mature seq is > $mature_lng_max nts\n"; return 0;}

    #number of reads that map in consistence with Dicer processing
    my $consistent=0;

    #number of reads that map inconsistent with Dicer processing
    my $inconsistent=0;

#   number of potential star reads map in good consistence with Drosha/Dicer processing
#   (3' overhangs relative to mature product)
    my $star_perfect=0;

#   number of potential star reads that do not map in good consistence with 3' overhang
    my $star_fuzzy=0;

    my $freq_mature=0;
    my $freq_loop=0;
    my $freq_star=0;


    #sort queries (deep sequences) by their position on the hairpin
    my @queries=sort {$hash_query{$a}{"db_beg"} <=> $hash_query{$b}{"db_beg"}} keys %hash_query;

    foreach my $query(@queries){

	#number of reads that the deep sequence represents
	unless(defined($hash_query{$query}{"freq"})){next;}
	my $query_freq=$hash_query{$query}{"freq"};

	#test which Dicer product (if any) the deep sequence corresponds to
	my $product=test_query($query);

	#and add the appropriate read counts
	if($product eq "mature"){$freq_mature+=$query_freq;}
	if($product eq "loop"){$freq_loop+=$query_freq;}
	if($product eq "star"){$freq_star+=$query_freq;}

	#if the deep sequence corresponds to a Dicer product, add to the 'consistent' variable
	if($product){$consistent+=$query_freq;}

	#if the deep sequence do not correspond to a Dicer product, add to the 'inconsistent' variable
	else{$inconsistent+=$query_freq;}

	#test a potential star sequence has good 3' overhang
	if($product eq "star"){
	    if(test_star($query)){$star_perfect+=$query_freq;}
	    else{$star_fuzzy+=$query_freq;}
	}
    }

#   if the majority of potential star sequences map in good accordance with 3' overhang
#    score for the presence of star evidence
    if($star_perfect>$star_fuzzy){$hash_comp{"star_read"}=1;}
    #find the most frequent star positions (this is opposed to star_expm which are the
    #positions that would be expected from the positions of the mature sequence and
    #the model of Dicer processing
    if(0<$star_perfect+$star_fuzzy){observed_star()};

    #number of reads that correspond to mature, loop, star and the sum of these
    $hash_comp{"freq_mature"}=$freq_mature;
    $hash_comp{"freq_loop"}=$freq_loop;
    $hash_comp{"freq_star"}=$freq_star;
    $hash_comp{"freq_total"}=$consistent;

    unless($consistent+$inconsistent>0){$hash_comp{"problem_read_cnt"}="no reads map to the candidate\n"; return 0;}

    #unless >90% of the reads map in consistence with Dicer processing, the hairpin is discarded
    my $inconsistent_fraction=$inconsistent/($inconsistent+$consistent);
    unless($inconsistent_fraction<=0.1){$hash_comp{"problem_Dicer_signature"}="more than 10% of the reads map inconsistent with Dicer processing\n"; return 0;}

    #the hairpin is retained
    return 1;
}

sub test_star{

    #test if a deep sequence maps in good consistence with 3' overhang

    my $query=shift;

    #5' begin and 3' end positions
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
    my $freq=$hash_query{$query}{"freq"};

    $hash_stars{$beg}{$end}+=$freq;

    #the difference between observed and expected begin positions must be 0 or 1
    my $offset=$beg-$hash_comp{"star_beg_exp"};
    if($offset==0 or $offset==1 or $offset==-1){return 1;}

    return 0;
}


sub observed_star{

    my $beg_max;
    my $end_max;
    my $freq_max=0;

    my @begs=sort {$a<=>$b} keys %hash_stars;
    foreach my $beg(@begs){

	my @ends=sort {$b<=>$a} keys %{$hash_stars{$beg}};
	foreach my $end(@ends){

	    my $freq=$hash_stars{$beg}{$end};
	    if($freq_max<$freq){

		$beg_max=$beg;
		$end_max=$end;
		$freq_max=$freq;
	    }
	}
    }

    $hash_comp{"star_beg_obs"}=$beg_max;
    $hash_comp{"star_end_obs"}=$end_max;
    $hash_comp{"star_freq_consensus"}=$freq_max;

    return;
}





sub test_query{

    #test if deep sequence maps in consistence with Dicer processing

    my $query=shift;

    #begin, end, strand and read count
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
    my $strand=$hash_query{$query}{"strand"};
    my $freq=$hash_query{$query}{"freq"};

    #should not be on the minus strand (although this has in fact anecdotally been observed for known miRNAs)
    if($strand eq '-'){return 0;}

    #the deep sequence is allowed to stretch 2 nt beyond the expected 5' end
    my $fuzz_beg=2;
    #the deep sequence is allowed to stretch 5 nt beyond the expected 3' end
    my $fuzz_end=5;

    #if in accordance with Dicer processing, return the type of Dicer product
    if(contained($beg,$end,$hash_comp{"mature_beg"}-$fuzz_beg,$hash_comp{"mature_end"}+$fuzz_end)){return "mature";}
    if(contained($beg,$end,$hash_comp{"star_beg_exp"}-$fuzz_beg,$hash_comp{"star_end_exp"}+$fuzz_end)){return "star";}
    if(contained($beg,$end,$hash_comp{"loop_beg"}-$fuzz_beg,$hash_comp{"loop_end"}+$fuzz_end)){return "loop";}

    #if not in accordance, return 0
    return 0;
}


sub pass_filtering_structure{

    #The potential precursor must form a hairpin with miRNA precursor-like characteristics

    #return value
    my $ret=1;

    #potential mature, star, loop and lower flank parts must be identifiable
    unless(test_components()){return 0;}

    #no bifurcations
    unless(no_bifurcations_precursor()){$ret=0;}

    #minimum 60% base pairings in duplex
    unless(bp_duplex()>=0.6){$ret=0; $hash_comp{"problem_duplex_bp"}="less than 60% of the nts in the candidate mature seq are base paired in the duplex\n";}

    #not more than 6 nt difference between mature and star length
    unless(-6<diff_lng() and diff_lng()<6){$ret=0; $hash_comp{"problem_duplex_lng"}="the difference between the candidate mature and star sequence is > 6 nts\n";}

    return $ret;
}






sub test_components{

    #tests whether potential mature, star, loop and lower flank parts are identifiable

    unless($hash_comp{"mature_struct"}){
	$hash_comp{"problem_struct_mature"}="no mature sequence could be identified\n";
	return 0;
    }

    unless($hash_comp{"star_struct"}){
	$hash_comp{"problem_struct_star"}="the candidate mature sequence does not form part of a likely miRNA duplex\n";
	return 0;
    }

    unless($hash_comp{"loop_struct"}){
	$hash_comp{"problem_struct_loop"}="no loop sequence could be identified\n";
   	return 0;
    }

    unless($hash_comp{"flank_first_struct"}){
	$hash_comp{"problem_struct_us_flanks"}="no upstream flanking sequence could be identified\n";
   	return 0;
    }

    unless($hash_comp{"flank_second_struct"}){
	$hash_comp{"problem_struct_ds_flanks"}="no downstream flanking sequence could be identified\n";
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
	$hash_comp{"problem_struct_bifurcation"}="there are bifurcations in the candidate precursor structure\n";
	return 0;
    }

    return 1;
}

sub bp_precursor{

    #total number of bps in the precursor

    my $pre_struct=$hash_comp{"pre_struct"};

    #simple pattern matching
    my $pre_bps=0;
    while($pre_struct=~/\(/g){
	$pre_bps++;
    }
    return $pre_bps;
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



sub do_test_assemble{

#    not currently used, tests if the 'pri_struct' as assembled from the parts (Dicer products, lower flanks)
#    is identical to 'pri_struct' before disassembly into parts

    my $assemble_struct;

    if($hash_comp{"flank_first_struct"} and $hash_comp{"mature_struct"} and $hash_comp{"loop_struct"} and $hash_comp{"star_struct"} and $hash_comp{"flank_second_struct"}){
	if($hash_comp{"mature_arm"} eq "first"){
	    $assemble_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"star_struct"}.$hash_comp{"flank_second_struct"};
	}else{
	    $assemble_struct.=$hash_comp{"flank_first_struct"}.$hash_comp{"star_struct"}.$hash_comp{"loop_struct"}.$hash_comp{"mature_struct"}.$hash_comp{"flank_second_struct"};
	}
	unless($assemble_struct eq $hash_comp{"pri_struct"}){
	    $hash_comp{"test_assemble"}=$assemble_struct;
	    print_hash_comp();
	}
    }
    return;
 }



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



sub fill_star{

    #fills specifics on the expected star strand into 'comp' hash ('component' hash)

    #if the mature sequence is not plausible, don't look for the star arm
    my $mature_arm=$hash_comp{"mature_arm"};
    unless($mature_arm){$hash_comp{"star_arm"}=0; return;}

    #if the star sequence is not plausible, don't fill into the hash
    my($star_beg_exp,$star_end_exp)=find_star();
    my $star_arm=arm_star($star_beg_exp,$star_end_exp);
    unless($star_arm){return;}

    #excise expected star sequence and structure
    my $star_seq=excise_seq($hash_comp{"pri_seq"},$star_beg_exp,$star_end_exp,"+");
    my $star_struct=excise_seq($hash_comp{"pri_struct"},$star_beg_exp,$star_end_exp,"+");

    #fill into hash
    $hash_comp{"star_beg_exp"}=$star_beg_exp;
    $hash_comp{"star_end_exp"}=$star_end_exp;
    $hash_comp{"star_seq"}=$star_seq;
    $hash_comp{"star_struct"}=$star_struct;
    $hash_comp{"star_arm"}=$star_arm;

    return;
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
    my $star_beg_exp=$offset_star_beg-$offset_beg;

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
    my $star_end_exp=$offset_star_end+$offset_end+2;

    return($star_beg_exp,$star_end_exp);
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
	$loop_beg=$hash_comp{"star_end_exp"}+1;
    }else{
	$loop_end=$hash_comp{"star_beg_exp"}-1;
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


sub fill_lower_flanks{

    #fills specifics on the lower flanks and unpaired strands into the 'comp' hash

    #unless both mature and star sequences are plausible, do not look for the flanks
    unless($hash_comp{"mature_arm"} and $hash_comp{"star_arm"}){return;}

    my $flank_first_end;
    my $flank_second_beg;

    #defining the begin and end positions of the flanks from the mature and star positions
    #excision depends on whether the mature or star sequence is 5' in the potenitial precursor ('first')
    if($hash_comp{"mature_arm"} eq "first"){
	$flank_first_end=$hash_comp{"mature_beg"}-1;
    }else{
	$flank_second_beg=$hash_comp{"mature_end"}+1;
    }

    if($hash_comp{"star_arm"} eq "first"){
	$flank_first_end=$hash_comp{"star_beg_exp"}-1;
    }else{
	$flank_second_beg=$hash_comp{"star_end_exp"}+1;
    }

    #unless the positions are plausible, do not fill into hash
    unless(test_flanks($flank_first_end,$flank_second_beg)){return;}

    $hash_comp{"flank_first_end"}=$flank_first_end;
    $hash_comp{"flank_second_beg"}=$flank_second_beg;
    $hash_comp{"flank_first_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    $hash_comp{"flank_second_seq"}=excise_seq($hash_comp{"pri_seq"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");
    $hash_comp{"flank_first_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"pri_beg"},$hash_comp{"flank_first_end"},"+");
    $hash_comp{"flank_second_struct"}=excise_struct($hash_comp{"pri_struct"},$hash_comp{"flank_second_beg"},$hash_comp{"pri_end"},"+");

    if($options{z}){
	fill_stems_drosha();
    }

    return;
}


sub fill_stems_drosha{

    #scores the number of base pairings formed by the first ten nt of the lower stems
    #in general, the more stems, the higher the score contribution
    #warning: this options has not been thoroughly tested

    my $flank_first_struct=$hash_comp{"flank_first_struct"};
    my $flank_second_struct=$hash_comp{"flank_second_struct"};

    my $stem_first=substr($flank_first_struct,-10);
    my $stem_second=substr($flank_second_struct,0,10);

    my $stem_bp_first=0;
    my $stem_bp_second=0;

    #find base pairings by simple pattern matching
    while($stem_first=~/\(/g){
	$stem_bp_first++;
    }

    while($stem_second=~/\)/g){
	$stem_bp_second++;
    }

    my $stem_bp=min2($stem_bp_first,$stem_bp_second);

    $hash_comp{"stem_first"}=$stem_first;
    $hash_comp{"stem_second"}=$stem_second;
    $hash_comp{"stem_bp_first"}=$stem_bp_first;
    $hash_comp{"stem_bp_second"}=$stem_bp_second;
    $hash_comp{"stem_bp"}=$stem_bp;

    return;
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


sub test_loop{

    #tests the loop positions

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}


sub test_flanks{

    #tests the positions of the lower flanks

    my ($beg,$end)=@_;

    #unless the begin and end positions are plausible, test negative
    unless($beg>0 and $beg<=$hash_comp{"pri_end"} and $end>0 and $end<=$hash_comp{"pri_end"} and $beg<=$end){return 0;}

    return 1;
}


sub comp{

    #subroutine to retrive from the 'comp' hash

    my $type=shift;
    my $component=$hash_comp{$type};
    return $component;
}


sub find_strand_query{

    #subroutine to find the strand for a given query

    my $query=shift;
    my $strand=$hash_query{$query}{"strand"};
    return $strand;
}


sub find_positions_query{

    #subroutine to find the begin and end positions for a given query

    my $query=shift;
    my $beg=$hash_query{$query}{"db_beg"};
    my $end=$hash_query{$query}{"db_end"};
    return ($beg,$end);
}



sub find_mature_query{

    #finds the query with the highest frequency of reads and returns it
    #is used to determine the positions of the potential mature sequence

    my @queries=sort {$hash_query{$b}{"freq"} <=> $hash_query{$a}{"freq"}} keys %hash_query;
    my $mature_query=$queries[0];
    return $mature_query;
}




sub reset_variables{

    #resets the hashes for the next potential precursor

    %hash_query=();
    %hash_comp=();
    %hash_bp=();
    %hash_stars=();

    $lines=();

    return;
}



sub excise_seq{

    #excise sub sequence from the potential precursor

    my($seq,$beg,$end,$strand)=@_;

    #begin can be equal to end if only one nucleotide is excised
    unless($beg<=$end){print STDERR "begin can not be smaller than end for $db_old\n";exit;}

    #rarely, permuted combinations of signature and structure cause out of bound excision errors.
    #this happens once appr. every two thousand combinations
    unless($beg<=length($seq)){$out_of_bound++;return 0;}

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


sub create_hash_seeds{

    #parses a fasta file with sequences of known miRNAs considered for conservation purposes
    #reads the seeds into a hash

    my ($file) = @_;
    my ($id, $desc, $sequence, $seed) = ();

    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $sequence = "";
	    $seed  = "";
	    while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){
		    $seed                = substr($sequence,1,$seed_lng);
		    $seed                =~ tr/[T]/[U]/;
		    $hash_mirs{$seed}   .="$id\t";
		    $hash_seeds{$seed} += 1;

		    $id               = $1;
		    $desc             = $2;
		    $sequence         = "";
		    $seed          = "";
		    next;
                }
		$sequence .= $_;
            }
        }
    }
    $seed                = substr($sequence,1,$seed_lng);
    $seed                =~ tr/[T]/[U]/;
    $hash_mirs{$seed}   .="$id\t";
    $hash_seeds{$seed} += 1;
    close FASTA;
}


sub parse_file_randfold{

    my $file=shift;

    open (FILE, "<$file") or die "can not open $file\n";
    while (<FILE>){

	if(/^(\S+)\s+(\S+)\s+(\S+)/){
	    my $id=$1;
	    my $randfold=$3;
	    $hash_randfold{$id}=$randfold;
	}
    }
    close FILE;
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



sub find_freq{

    #finds the frequency of a given read query from its id.

    my($query)=@_;

    if($query=~/_x(\d+)/){
	my $freq=$1;
	return $freq;
    }else{
	return 0;
    }
}

sub print_hash_comp{
    my @reads;
    my $lr;
    my $bef;

    my $after;
    my $col1_width = 40;
	$col1_width=$options{'l'} if($options{'l'});
	$col1_width+=5;
    my $rseq;

    my @sread;
    $hash_comp{"pri_seq"} =~ tr/Tt/Uu/;
    my @pseq = split(//,lc $hash_comp{"pri_seq"});

    my $spacer = " " x ($col1_width - length("pri_seq"));
    my $shift;
    my %reads_hash;

    print ">$hash_comp{\"pri_id\"}\n";
    my $spacer2 = 14;
    my $spacer2s=0;

    if(not $hash_comp{"problem_structure"} and not $hash_comp{"problem_signature"}){
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score"}));
        print "score total\t\t$spacer2s$hash_comp{\"score\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_star"}));
        print "score for star read(s)\t$spacer2s$hash_comp{\"score_star\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_freq"}));
        print "score for read counts\t$spacer2s$hash_comp{\"score_freq\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_mfe"}));
        print "score for mfe\t\t$spacer2s$hash_comp{\"score_mfe\"}\n";

        if($options{'y'}){
            $spacer2s = " " x ($spacer2 - length($hash_comp{"score_randfold"}));
            print "score for randfold\t$spacer2s$hash_comp{\"score_randfold\"}\n";
        }
	if($options{'s'}){
        $spacer2s = " " x ($spacer2 - length($hash_comp{"score_seed"}));
        print "score for cons. seed\t$spacer2s$hash_comp{\"score_seed\"}\n";
    }
	if($options{'s'} and defined $hash_comp{"seed_family"}){
        $spacer2s = " " x ($spacer2 - length($hash_comp{"seed_family"}));
        print "miRNA with same seed\t$spacer2s$hash_comp{\"seed_family\"}\n";

    }
        $spacer2s = " " x ($spacer2 - length($hash_comp{"freq_total"}));
        print "total read count\t$spacer2s$hash_comp{\"freq_total\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"freq_mature"}));
        print "mature read count\t$spacer2s$hash_comp{\"freq_mature\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"freq_loop"}));
        print "loop read count\t\t$spacer2s$hash_comp{\"freq_loop\"}\n";
        $spacer2s = " " x ($spacer2 - length($hash_comp{"freq_star"}));
        print "star read count\t\t$spacer2s$hash_comp{\"freq_star\"}\n";


    }else{
	#structure problems
	if(defined $hash_comp{"problem_structure"}){ print $hash_comp{"problem_structure"};}
	if(defined $hash_comp{"problem_mature_lng"}){ print $hash_comp{"problem_mature_lng"};}
	if(defined $hash_comp{"problem_duplex_bp"}){ print $hash_comp{"problem_duplex_bp"};}
	if(defined $hash_comp{"problem_duplex_lng"}){ print $hash_comp{"problem_duplex_lng"};}
	if(defined $hash_comp{"problem_struct_mature"}){ print $hash_comp{"problem_struct_mature"};}
	if(defined $hash_comp{"problem_struct_star"}){ print $hash_comp{"problem_struct_star"};}
	if(defined $hash_comp{"problem_struct_loop"}){ print $hash_comp{"problem_struct_loop"};}
	if(defined $hash_comp{"problem_struct_us_flanks"}){ print $hash_comp{"problem_struct_us_flanks"};}
	if(defined $hash_comp{"problem_struct_ds_flanks"}){ print $hash_comp{"problem_struct_ds_flanks"};}
	if(defined $hash_comp{"problem_struct_bifurcation"}){ print $hash_comp{"problem_struct_bifurcation"};}

	#signature problems
	if(defined $hash_comp{"problem_signature"}){ print $hash_comp{"problem_signature"};}
	if(defined $hash_comp{"problem_read_cnt"}){ print $hash_comp{"problem_read_cnt"};}
	if(defined $hash_comp{"problem_Dicer_signature"}){ print $hash_comp{"problem_Dicer_signature"};}
    }


    #### printing alignment;
    print "exp";

    print " " x ($col1_width-3);

    if($hash_comp{"problem_structure"}){
	print "n" x length($hash_comp{"pri_seq"});
    }else{
	print "f" x $hash_comp{"flank_first_end"} ;
	if($hash_comp{"mature_beg"}  < $hash_comp{"loop_beg"}){
	    $bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
	    print $bef;
	    $bef = "l" x ($hash_comp{"loop_end"}-$hash_comp{"loop_beg"}+1);
	    print $bef;
	    $bef = "S" x ($hash_comp{"star_end_exp"}-$hash_comp{"star_beg_exp"}+1);
	    print $bef;
	}else{
	    $bef = "S" x ($hash_comp{"star_end_exp"}-$hash_comp{"star_beg_exp"}+1);
	    print $bef;
	    $bef = "l" x ($hash_comp{"loop_end"}-$hash_comp{"loop_beg"}+1);
	    print $bef;

	    $bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
	    print $bef;
	}
	$bef = "f" x ($hash_comp{"pri_end"}-$hash_comp{"flank_second_beg"}+1);
	print $bef;
    }

    if(defined $hash_comp{"star_beg_obs"}){

	print "\nobs";

	print " " x ($col1_width-3);

	if($hash_comp{"problem_structure"}){
	    print "n" x length($hash_comp{"pri_seq"});
	}else{
	    if($hash_comp{"mature_beg"}  < $hash_comp{"loop_beg"}){

		$bef="f" x ($hash_comp{"mature_beg"}-1);
		print $bef;
		$bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
		print $bef;
		$bef = "l" x ($hash_comp{"star_beg_obs"}-$hash_comp{"mature_end"}-1);
		print $bef;
		$bef = "S" x ($hash_comp{"star_end_obs"}-$hash_comp{"star_beg_obs"}+1);
		print $bef;
		$bef="f" x ($hash_comp{"pri_end"}-$hash_comp{"star_end_obs"});
		print $bef;

	    }else{

		$bef="f" x ($hash_comp{"star_beg_obs"}-1);
		print $bef;
		$bef = "S" x ($hash_comp{"star_end_obs"}-$hash_comp{"star_beg_obs"}+1);
		print $bef;
		$bef = "l" x ($hash_comp{"mature_beg"}-$hash_comp{"star_end_obs"}-1);
		print $bef;
		$bef = "M" x ($hash_comp{"mature_end"}-$hash_comp{"mature_beg"}+1);
		print $bef;
		$bef="f" x ($hash_comp{"pri_end"}-$hash_comp{"mature_end"});
		print $bef;
	    }
	}
    }

    print "\npri_seq$spacer",lc $hash_comp{"pri_seq"},"\n";
    $spacer = " " x ($col1_width - length("pri_struct"));
    print "pri_struct$spacer$hash_comp{\"pri_struct\"}\t#MM\n";
    @reads = split(/\n/,$lines);
    $lr = scalar @reads;
    my $rrc = 0;
    foreach(@reads){
        if(/^(\S+)\s+\d+\s+\d+\s+\d+\s+(\S+)\s+\S+\s+\d+\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\d+).+$/){
            $rrc++;
            $reads_hash{$rrc}{"id"}=$1;
            $reads_hash{$rrc}{"seq"}=$2;
            $reads_hash{$rrc}{"beg"}=$3;
            $reads_hash{$rrc}{"end"}=$4;
            $reads_hash{$rrc}{"mm"}=$5;
        }
    }

    ## sorted keys by begin postion
    my @skeys = sort { $reads_hash{$a}{"beg"} <=> $reads_hash{$b}{"beg"} } keys %reads_hash;
    my @elist; # final sorted array

    my $first = $reads_hash{$skeys[0]}{"beg"};  ## all keys that have same begin position should match this value
    my %rorder;                                 ## temporary hash to store all keys with same begin position

    for(my $j = 0; $j < scalar @skeys; $j++){
        if($reads_hash{$skeys[$j]}{"beg"} eq $first){
            $rorder{$j} = $reads_hash{$skeys[$j]}{"end"};  ## insert key and end position to hash
        }else{                                             ## if new begin position
            $first = $reads_hash{$skeys[$j]}{"beg"};
            for(sort {$rorder{$a} <=> $rorder{$b}} keys %rorder){ ## sort hash keys by end position
                push(@elist,$skeys[$_]);                          ## attend keys to elist
            }
            for(keys %rorder){delete $rorder{$_};}                ## delete hash
            $rorder{$j} = $reads_hash{$skeys[$j]}{"end"};
        }
    }

    for(sort {$rorder{$a} <=> $rorder{$b}} keys %rorder){
        push(@elist,$skeys[$_]);
    }

    foreach(@elist){                                                       ## output elist.
        $rseq  = lc $reads_hash{$_}{'seq'};
        $rseq =~ tr/t/u/;
        $bef="." x ($reads_hash{$_}{'beg'}-1);
        $after = "." x ($hash_comp{'pri_end'} - $reads_hash{$_}{"end"});
        $spacer = " " x ($col1_width - length($reads_hash{$_}{'id'}));
        @sread = split(//,$rseq);


        my $bshift = 0;
        $rseq = "";
        for(my $i=0; $i < scalar @sread; $i++){
            if(not $pseq[$i+$reads_hash{$_}{'beg'}-1]){                        ### read is longer than sequence
            }else{
                if($pseq[$i+$reads_hash{$_}{'beg'}-1] eq $sread[$i]){
                    $rseq .=  lc $sread[$i];
                }else{

                        $sread[$i] = uc $sread[$i];
                        $rseq .= uc $sread[$i];
                }

            }
        }
        print "$reads_hash{$_}{'id'}$spacer$bef$rseq$after\t$reads_hash{$_}{'mm'}\n";
    }
    print "\n\n\n";
}




sub print_hash_bp{

    #prints the 'bp' hash

    my @keys=sort {$a<=>$b} keys %hash_bp;
    foreach my $key(@keys){
	my $value=$hash_bp{$key};
	print "$key\t$value\n";
    }
    print "\n";
}



sub contained{

    #Is the stretch defined by the first positions contained in the stretch defined by the second?

    my($beg1,$end1,$beg2,$end2)=@_;

    testbeginend($beg1,$end1,$beg2,$end2);

    if($beg2<=$beg1 and $end1<=$end2){
	return 1;
    }else{
	return 0;
    }
}


sub testbeginend{

    #Are the beginposition numerically smaller than the endposition for each pair?

    my($begin1,$end1,$begin2,$end2)=@_;

    unless($begin1<=$end1 and $begin2<=$end2){
	print STDERR "beg can not be larger than end for $db_old\n";
	exit;
    }
}


sub rev_pos{

#   This subroutine reverses the begin and end positions

    my($beg,$end,$lng)=@_;

    my $new_end=$lng-$beg+1;
    my $new_beg=$lng-$end+1;

    return($new_beg,$new_end);
}

sub round {

    #rounds to nearest integer

    my($number) = shift;
    return int($number + .5);

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

sub revcom{

    #reverse complement

    my($sequence)=@_;

    my $revcom=rev(com($sequence));

    return $revcom;
}


sub max2 {

    #max of two numbers

    my($a, $b) = @_;
    return ($a>$b ? $a : $b);
}

sub min2  {

    #min of two numbers

    my($a, $b) = @_;
    return ($a<$b ? $a : $b);
}



sub score_freq{

#   scores the count of reads that map to the potential precursor
#   Assumes geometric distribution as described in methods section of manuscript

    my $freq=shift;

    #parameters of known precursors and background hairpins
    my $parameter_test=0.999;
    my $parameter_control=0.6;

    #log_odds calculated directly to avoid underflow
    my $intercept=log((1-$parameter_test)/(1-$parameter_control));
    my $slope=log($parameter_test/$parameter_control);
    my $log_odds=$slope*$freq+$intercept;

    #if no strong evidence for 3' overhangs, limit the score contribution to 0
    unless($options{x} or $hash_comp{"star_read"}){$log_odds=min2($log_odds,0);}

    return $log_odds;
}



sub score_mfe{

#   scores the minimum free energy in kCal/mol of the potential precursor
#   Assumes Gumbel distribution as described in methods section of manuscript

    my $mfe=shift;

    #numerical value, minimum 1
    my $mfe_adj=max2(1,-$mfe);

    #parameters of known precursors and background hairpins, scale and location
    my $prob_test=prob_gumbel_discretized($mfe_adj,5.5,32);
    my $prob_background=prob_gumbel_discretized($mfe_adj,4.8,23);

    my $odds=$prob_test/$prob_background;
    my $log_odds=log($odds);

    return $log_odds;
}



sub prob_gumbel_discretized{

#   discretized Gumbel distribution, probabilities within windows of 1 kCal/mol
#   uses the subroutine that calculates the cdf to find the probabilities

    my ($var,$scale,$location)=@_;

    my $bound_lower=$var-0.5;
    my $bound_upper=$var+0.5;

    my $cdf_lower=cdf_gumbel($bound_lower,$scale,$location);
    my $cdf_upper=cdf_gumbel($bound_upper,$scale,$location);

    my $prob=$cdf_upper-$cdf_lower;

    return $prob;
}


sub cdf_gumbel{

#   calculates the cumulative distribution function of the Gumbel distribution

    my ($var,$scale,$location)=@_;

    my $cdf=$e**(-($e**(-($var-$location)/$scale)));

    return $cdf;
}

