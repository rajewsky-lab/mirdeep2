#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Std;
use POSIX qw(ceil floor);



my $usage =
"$0 output_mirdeep

Options:

-a file  file outputted by controls
-b file  mature miRNA fasta reference file for the species
-c file  signature file
-d int   read stack height necessary for triggering excision
";

my $file_out=shift or die $usage;

my $rounds_controls=100;


my %options=();
getopts("a:b:c:d:",\%options);

my %hash_con;

my %hash_ref;

my %hash_sig;

my %hash_out;


#MAIN

if(($options{b} and not $options{c}) or (not $options{b} and $options{c})){print STDERR "options -b and -c must be used in conjunction\n"; exit;}

if($options{b} and $options{c}){

    parse_file_ref(\$options{b},\%hash_ref);

    parse_file_sig(\$options{c},\%hash_sig);

}

if($options{a}){parse_file_mrd_permuted(\$options{a},\%hash_con)}

parse_file_mrd(\$file_out,\%hash_out);

print_header();

survey();



#iterate over all scores from -10 to 10
sub survey{
    
    for(my $score=10; $score>=-10; $score--){
        
        print "$score";
	
	if($options{b}){
            
	    survey_hairpins($score);
	    
            survey_known($score);
	}
	
		if($options{a}){survey_signal_to_noise($score)};
	
		if($options{d}){
            
			my $read_stack_min=$options{d};
			
			print "\t$read_stack_min";
		}
		
        print "\n";
    }
    return;
}


#summary statistics for mature miRNAs
sub survey_known{

    my $score=shift;

    #mature miRNAs for the species in reference (miRBase) file

    my $matures_cnt=scalar (keys %hash_ref);

    print "\t$matures_cnt";

    #matures present in data

    my $matures_present=scalar (keys %hash_sig);

    print "\t$matures_present";
    
    #matures recovered

    my $matures_recov_cnt=0;
    
    my @matures_present=keys %hash_sig;

    foreach my $mature_present(@matures_present){
	
	my $score_mature_present=$hash_sig{$mature_present};
	
	if($score<=$score_mature_present){$matures_recov_cnt++;}
    }
    
    print "\t$matures_recov_cnt";
    
    #matures recovered in percent

    my $percent=100*round_decimals(div($matures_recov_cnt,$matures_present),2);	    

    print " ($percent%)";
    
    return;
}


sub survey_signal_to_noise{

    my $score=shift;
    
    #total hairpins reported
    my $hairpins_total=hairpins_cnt("total",$score);
    
#    print "\t$hairpins_total";
  
    my($hairpins_total_fp_mean,$hairpins_total_fp_sd,$estimated_total_true_mean,$estimated_total_true_sd,$percent_total_mean,$percent_total_sd)=mean_sd("total",$score,$hairpins_total);
    
    my $hairpins_total_fp_mean_round=round_decimals($hairpins_total_fp_mean,0);
    
    my $hairpins_total_fp_sd_round=round_decimals($hairpins_total_fp_sd,0);
    
#        print "\t$hairpins_total_fp_mean_round +/- $hairpins_total_fp_sd_round";
    
    my $signal_to_noise_total;
    if($hairpins_total_fp_mean eq 0){
	$signal_to_noise_total = 0;
    }else{
		$signal_to_noise_total=round_decimals(div($hairpins_total,$hairpins_total_fp_mean),1);
    }
    
    print "\t$signal_to_noise_total";
    
    
    return;
}






#summary statistics for miRNA hairpin precursors
sub survey_hairpins{

    my $score=shift;

    #partitioning of hairpins into known and novel
         
    my $hairpins_known=hairpins_cnt("known",$score);
    
#    print "\t$hairpins_known";
    
    my $hairpins_novel=hairpins_cnt("novel",$score);
    
    print "\t$hairpins_novel";
    
    #estimation of false positives for the set of novel hairpins
    if($options{a}){
		
		my ($hairpins_fp_mean,$hairpins_fp_sd,$estimated_true_mean,$estimated_true_sd,$percent_mean,$percent_sd)=(0,0,0,0,0,0);
		if($hairpins_novel){ ## check if novel hairpins detected 
			($hairpins_fp_mean,$hairpins_fp_sd,$estimated_true_mean,$estimated_true_sd,$percent_mean,$percent_sd)=mean_sd("novel",$score,$hairpins_novel);
		}
		my $hairpins_fp_mean_round=round_decimals($hairpins_fp_mean,0);
		
		my $hairpins_fp_sd_round=round_decimals($hairpins_fp_sd,0);
	
		print "\t$hairpins_fp_mean_round +/- $hairpins_fp_sd_round";
		
		my $estimated_true_mean_round=round_decimals($estimated_true_mean,0);
		
		print "\t$estimated_true_mean_round";
		
		my $estimated_true_sd_round=round_decimals($estimated_true_sd,0);
		
		print " +/- $estimated_true_sd_round";
		
		my $percent_mean_round=round_decimals($percent_mean,0);
		
		print " ($percent_mean_round";
		
		my $percent_sd_round=round_decimals($percent_sd,0);
		
		print " +/- $percent_sd_round%)";
		
    }
    return;
}


#report number of hairpins of a given type ("total","known" or "novel") for a given score cut-off
sub hairpins_cnt{

    my($type,$score)=@_;

    my $hairpins_cnt=0;

    my @scores_cnt=sort keys %{$hash_out{$type}};
    
    foreach my $score_cnt(@scores_cnt){
	
	if($score<=$score_cnt){
	    
	    my $count=$hash_out{$type}{$score_cnt};

	    $hairpins_cnt+=$count;
	}
    }
    return $hairpins_cnt;
}


#calculate all summary statistics (estimated number of false and true positives, standard deviations, percentages of total hairpins that are true)
#these statistics are calculated from the permuted controls 
sub mean_sd{

    my($type,$score,$norm)=@_;

    my @counts;

    my @counts_above_norm;

    my @percentages;

    my @permutations=sort keys %{$hash_con{$type}};

    #iterate over runs of permuted controls. Each run is separated by 'permutation int' in the output_permuted.mrd file.
    foreach my $permutation(@permutations){

		#number of hairpins that exceed the score cut-off for this run of permuted controls
		my $cnt_total=0;
		
		my @scores_cnt=sort keys %{$hash_con{$type}{$permutation}};
		
		foreach my $score_cnt(@scores_cnt){

			if($score<=$score_cnt){
				
				my $count=$hash_con{$type}{$permutation}{$score_cnt};
				
				$cnt_total+=$count;
				
			}
		}
		
		push(@counts,$cnt_total);
	}
	## error
    my $fp_mean=mean(@counts);
	


    my $fp_sd=sd(@counts);

    #iterate over runs of permuted controls again, using the @counts array just generated
    foreach my $count(@counts){
	
		#estimated number of true positives t for this run
		#t = number of hairpins reported - estimated number of false positives for this run
		my $cnt_total_above_norm=max2(0,$norm-$count);
		
		my $percentage=max2(0,100*(div($cnt_total_above_norm,$norm)));
		
		push(@counts_above_norm,$cnt_total_above_norm);
		
		push(@percentages,$percentage);
    }
    
    #calculate mean, sd, and percentages
	## error
    my $true_mean=mean(@counts_above_norm);

    my $true_sd=sd(@counts_above_norm);

    my $percentage_mean=mean(@percentages);

    my $percentage_sd=sd(@percentages);

    return ($fp_mean,$fp_sd,$true_mean,$true_sd,$percentage_mean,$percentage_sd);
}



sub print_header{

    print "miRDeep2 score";
    
    if($options{b}){

	print "\tnovel miRNAs reported by miRDeep2";

	if($options{a}){
	    
	    print "\tnovel miRNAs, estimated false positives";

	    print "\tnovel miRNAs, estimated true positives";

	}

	print "\tknown miRNAs in species";

	print "\tknown miRNAs in data";

	print "\tknown miRNAs detected by miRDeep2";
    }

    if($options{a}){print "\testimated signal-to-noise";}

    if($options{d}){print "\texcision gearing";}

    print "\n";

    return;
}


#parse the output.mrd file. Every time a '>' is encountered, all variables are resolved via a subroutine 
sub parse_file_mrd{

    my($file,$hash)=@_;

    my($score,@refs)=();

    open (FILE, "<$$file") or die "can not open $$file\n";
    
    while (my $line=<FILE>){

        

	if($line=~/^score\s+total\s+(\S+)/){   
	    
	    $score=$1;

	#the hairpin has reference miRBase miRNA mapping to it. This is used to partition known and novel hairpins.    
	}elsif($line=~/^(\S+)/ and defined $hash_sig{$1}){
	    
	    push(@refs,$1);
	    
	}elsif($line=~/^>/ and defined $score){    
	    
	    resolve_entry_file_mrd(\$score,\@refs,\$hash);
	}
    }
    resolve_entry_file_mrd(\$score,\@refs,\$hash);

    return;
}


#parse the output_permuted.mrd file. Every time a '>' is encountered, all variables are resolved via a subroutine
sub parse_file_mrd_permuted{

    my($file,$hash)=@_;

    my($score,$permutation,@refs)=();

    open (FILE, "<$$file") or die "can not open $$file\n";
    
    while (my $line=<FILE>){
		
		if($line=~/^permutation\s+(\d+)/){
			
			$permutation=$1;
			
		}elsif($line=~/^score\s+total\s+(\S+)/){   
			
			$score=$1;
			
			#the hairpin has reference miRBase miRNA mapping to it. This is used to partition known and novel hairpins.
		}elsif($line=~/^(\S+)/ and defined $hash_sig{$1}){
			
			push(@refs,$1);
			
		}elsif($line=~/^>/ and defined $score){    
			
			resolve_entry_file_mrd_permuted(\$score,\$permutation,\@refs,\$hash);
		}
    }
    resolve_entry_file_mrd_permuted(\$score,\$permutation,\@refs,\$hash);

    return;
}


#read all variables for one hairpin into hash
sub resolve_entry_file_mrd_permuted{

    my($score,$permutation,$refs,$hash)=@_;
 
    unless(defined $$permutation){
		print STDERR "The $options{a} file is not properly formatted.\nMaybe it does not contain the lines with \"permutation int\"?\n"; 
		exit;
	}
   
    #all scores with the same integer part are pooled (e.g. score 1.2 and 1.0).
    my $floor=floor($$score);
    
    $$$hash{"total"}{$$permutation}{$floor}++;
    
    #if the hairpin has one or more reference miRNAs mapping, it is treated as known hairpin for purposes of the controls
    if(@$refs){
		$$$hash{"known"}{$$permutation}{$floor}++;

    }else{
		
		$$$hash{"novel"}{$$permutation}{$floor}++;
    }
    
    ($$score,@$refs)=();
    
    return;
}


#read all variables for one hairpin into hash
sub resolve_entry_file_mrd{

    my($score,$refs,$hash)=@_;

    #all scores with the same integer part are pooled (e.g. score 1.2 and 1.0).
    my $floor=floor($$score);
    
    $$$hash{"total"}{$floor}++;
    
    #if the hairpin has one or more reference miRNAs mapping, it is a known hairpin
    if(@$refs){
	
	$$$hash{"known"}{$floor}++;

	foreach my $ref(@$refs){

	    #if the score for this hairpin is higher than the score for previous hairpins
	    #with the same miRBase mature mapping, use this score instead
	    if($hash_sig{$ref}<$floor){$hash_sig{$ref}=$floor;}
	}	

    }else{
	
	$$$hash{"novel"}{$floor}++;
    }
    
    ($$score,@$refs)=();
    
    return;
}



sub parse_file_ref{
    
    my ($file,$hash) = @_;
    my ($id, $desc, $seq) = ();

    open (FASTA, "<$$file") or die "can not open $$file\n";
    while (<FASTA>)
    {

        chomp;
        if (/^>(\S+)(.*)/)
	{
	    $id       = $1;
	    $desc     = $2;
	    $seq      = "";
	    while (<FASTA>){
                    
                chomp;
                if (/^>(\S+)(.*)/){
		    $$hash{$id} = $seq;
		    $id         = $1;
		    $desc       = $2;
		    $seq        = "";
		    next;
                }
		$seq  .= $_;
            }
        }
    }
    $$hash{$id} = $seq;

    close FASTA;
}

#initialize hash with the best scores for all mature miRBase miRNAs mapping to hairpins
sub parse_file_sig{

    my($file_sig,$hash_sig)=@_;
    
    open (FILE, "<$$file_sig") or die "can not open $$file_sig\n";

    while (my $line=<FILE>){

	if($line=~/^(\S+)/ and defined $hash_ref{$1}){

	    #dummy value
	    $hash_sig{$1}=-50;
	}
    }
    return;
}


#for testing purposes
sub print_hash_sig{

    my @mirs=sort keys %hash_sig;

    foreach my $mir(@mirs){
	
	my $score=$hash_sig{$mir};

	print "$mir\t$score\n";
    }
    return;
}


#calculate mean of an array
sub mean{
	
    my @array=@_;

    my $sum = 0;

    ($sum+=$_) for @array;

    my $array=scalar @array;

    return ($sum/$array);
}


#calculate standard deviation of an array
sub sd{

    my @array=@_;

    my $mean=mean(@array);

    my $sum_squares=0;

    foreach my $var(@array){

	my $dif=$mean-$var;

	my $square=$dif**2;

	$sum_squares+=$square;
    }

    my $freedom=(scalar @array)-1;
	my $sd=0;
	if($freedom != 0){$sd=sqrt ($sum_squares/$freedom)};

    return $sd;
}
	


sub round_decimals{

    my($number,$decimals)=@_;

    return round(10**$decimals*$number)/10**$decimals;
}


sub round {
    
    my($number) = shift;
    return int($number + .5);
    
}



sub div{

    my($num,$den)=@_;

    if($den==0){return "0";}

    else{return $num/$den;}
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
