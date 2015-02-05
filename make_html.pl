#!/usr/bin/perl


use strict;
use Getopt::Std;
use Cwd;


use PDF::API2;              ## needed to create a PDF file
use Math::Trig;             ## needed for sinus function
use List::Util qw(min max); ## needed to provide min and max function for lists  
use File::Basename;


Usage() if(not $ARGV[0] or $ARGV[0] =~ /-*-he*l*p*/);

my $n_count = 0;
my $k_count = 0;
my $e_count = 0;
my $sig = 0;

my $infile;         ## miRDeep output file
my $pdfs = 0;       ## force pdf creation
my $threshold= 0;   ## hairpins have to have score above threshold in order to be reported
my $csv = 0;

my $xposshift = 40;
my $lstruct_multi;
my $sb_obs;

## read in available organisms at the end of the script

my %organisms;
while(<DATA>){
    chomp;
    my $tmp;
    $tmp=$_;
    $_ =~ s/\s+//g;
    $organisms{$_}=$tmp;
}
my $known;

## read in miRDeep output file
my $id;
my %hash;   ## takes up all entries from the miRDeep2 module


my $created = 1;

my $in;
my $counter=0;

my %struct; ## counts how often a nt is covered reads

my $i;

my $offset = 0;

my $me=0;     ## mature end coordinate
my @desc;

my $lflank1; ## length(string of left flank)
my $fl1;    ## 
my $lflank2; ## length string of right flank

my $fl2b=0;   ## right flank begin  
my $lloop;   ## string of loop
my $lb=0;     ## starting position of loop

my $lstar;   ## string of star sequence
my $sb=0;     ## starting 
my $lmature; ## string of mature 

my $mb=0;     ## mature begin
my $struct; ## structure string
my $pri_seq;## pri-cursor sequence
my $lenstr=0; 

my $pdf;    ## pdf descriptor
my $page;   ## page descriptor
my $gfx;    ## graphic variable

my $trb;    ## fontvariable
my $text;
my $text2;


my $aligned;                          ## reads reading 
my %hash2;                            ## begin of read in precursor

my %hash2c;                           ## number of reads per read sequence 
my %hash2key;
my %hash2mm;                          ## number of mismatches
my %hash2order;                       ## output order saved
my %hash2seq;
my %hash2sample;

my %order;                            ## stores begin coordinates of fl1,m,l,s,fl2 sequences
my $multiplier = 3.6;#4.825;               ## minimal distance between two letters


## calculate predefined pdf loci for alignment characters
my %position_hash;
$counter = 0;
for(my $i=0;$i < 200; $i++){
    $position_hash{$counter} = $xposshift+$multiplier+20+$i*$multiplier;
    $counter++;
}


my $yorig = 500; ## 500
my $downy = 50;

my $dline;                            ## line graphic handler

my $first=1;
my $lastx;
my $lasty;

my $final;                            ## final output string of a read
my @pseq;                             ## precursor sequence  
my @rseq;                             ## read sequence

my $totalreads = 0;
                     
my %assign_str;                       ## color assigned to position where letter is drawn
my %assign_str_exp;

my $bpo1=-10;                             ## left nt pos in first bp 
my $bpo2=-10;                             ## right nt pos in first bp 
my $bpo1r=-10;                            ## left nt pos in second bp 
my $bpo2r=-10;                            ## right nt pos in second bp 


my $ffe=0;                              ## first flank end position
my $ff2b=0;                             ## second flank begin position

my @sorted;                           ## array that stores sorted order of fl1,m,l,s,fl2 
my $y=$yorig;                         ## y coordinate


my ($minx,$miny,$maxx,$maxy);        ## min and max x,y coordinates of rna sequence
my @rna;                 ## rna sequence
my @rna_d;
my %xc;                  ## holds x cooridnate of each nt
my %yc;                  ## holds y coordinate of each nt
my $sid="";        


my %pres_coords;


## pdf histogram colors
my $col_star_exp = 'lightskyblue';
my $col_star_obs = 'darkviolet';
my $col_mature = 'red';
my $col_loop = 'orange';


## options
my %options=();
getopts("ug:v:f:ck:os:t:er:q:dx:zy:ab:p:V:E",\%options);

if($options{'u'}){
    print "\n\nAvailable species organisms were:\n\n";
   for(keys %organisms){
        print "$_\n";
    }
    print "\n\n\n";
    exit;
}



my %pres_coords;
if($options{'p'}){
	get_precursor_pos();
}


my $time = $options{'y'} or die "no timestamp given with parameter y\n";

if($options{'x'} and not $options{'q'}){
    die "\nError:\n\toption -x can only be used together with option -q\n\n";
}

## obtain current working directory
my $cwd = cwd;

## order output by sample (give -o option) or just by beginning position (no -o option)


## organism parameter
my $org=$organisms{$options{'t'}};


## some quantifier variables
my $mirbase = 0;
my %mature2hairpin;
my %hairpin2mature;  ## some hairpins have more than 1 mature assigned, circumvent this problem
my %hash_q; ## takes up all entries from the quantifier module
 
my $blast="http://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY=";
my $blast_query = "&db=nucleotide&QUERY_FROM=&QUERY_TO=&QUERYFILE=&GENETIC_CODE=1&SUBJECTS=&stype=nucleotide&SUBJECTS_FROM=&SUBJECTS_TO=&SUBJECTFILE=&DBTYPE=gc&DATABASE=nr&EQ_MENU=&NUM_ORG=1&EQ_TEXT=&BLAST_PROGRAMS=blastn&PHI_PATTERN=&MAX_NUM_SEQ=100&SHORT_QUERY_ADJUST=on&EXPECT=10&WORD_SIZE=7&MATRIX_NAME=PAM30&MATCH_SCORES=2,-3&GAPCOSTS=5+2&COMPOSITION_BASED_STATISTICS=0&FILTER=L&REPEATS=repeat_9606&FILTER=m&TEMPLATE_LENGTH=0&TEMPLATE_TYPE=0&PSSM=&I_THRESH=&SHOW_OVERVIEW=true&SHOW_LINKOUT=true&GET_SEQUENCE=auauauaauauauauauauuauaa&FORMAT_OBJECT=Alignment&FORMAT_TYPE=HTML&ALIGNMENT_VIEW=Pairwise&MASK_CHAR=2&MASK_COLOR=1&DESCRIPTIONS=100&ALIGNMENTS=100&NEW_VIEW=true&OLD_BLAST=false&NCBI_GI=false&SHOW_CDS_FEATURE=false&NUM_OVERVIEW=100&FORMAT_EQ_TEXT=&FORMAT_ORGANISM=&EXPECT_LOW=&EXPECT_HIGH=&QUERY_INDEX=&CLIENT=web&SERVICE=plain&CMD=request&PAGE=Nucleotides&PROGRAM=blastn&MEGABLAST=&RUN_PSIBLAST=&TWO_HITS=&DEFAULT_PROG=megaBlast&WWW_BLAST_TYPE=&DB_ABBR=&SAVED_PSSM=&SELECTED_PROG_TYPE=blastn&SAVED_SEARCH=true&BLAST_SPEC=&QUERY_BELIEVE_DEFLINE=&DB_DIR_PREFIX=&USER_DATABASE=&USER_WORD_SIZE=&USER_MATCH_SCORES=&USER_FORMAT_DEFAULTS=&NO_COMMON=&NUM_DIFFS=2&NUM_OPTS_DIFFS=1&UNIQ_DEFAULTS_NAME=A_SearchDefaults_1Mn7ZD_2Sq4_1Z58HQ5Jb_23tpbD_167y9p&PAGE_TYPE=BlastSearch&USER_DEFAULT_PROG_TYPE=blastn&USER_DEFAULT_MATCH_SCORES=3.";

## get mature positions if options a is set
my %mature_pos_hash;
if($options{'a'} and -f "mirdeep_runs/run_$time/tmp/signature.arf"){
    get_mature_pos();
}

my %confident =();
if($options{b}){
    open IN,"<$options{'b'}" or die "file with confidential ids not found\n";
    while(<IN>){
        next if(/\#/);
        chomp;
        my $r = $_;
        $r =~ s/\|/_/g;
        $confident{$r} = 1;
    }
    close IN;
}


## this option is only set if the quantifier module is used alone!!!
if($options{'z'}){
    PrintKnownnotfound();
    CloseHTML();
    if(not $options{'d'}){
        $mirbase = 1;
        CreateStructurePDF(%hash_q);
    }

    exit;
}



## scan program parameters
if($options{'f'}){
	$infile = $options{'f'};                 
	
}else{
    print STDERR "Error: no miRDeep2 output file specified by -f option\n\n";
	
    Usage();
}


if(not $options{'t'}){
    print STDERR "\nNo organism specified by switch -t [organism] so UCSC browser BLAT will not be used.

To get a list of possible species names type make_html.pl -u 

Specify organism Reads are coming from by  -t [organism] to get links in html file
\n";
}


if($options{'v'}){
	$threshold = $options{'v'};
}



## read in known miRNA identifier
my %knownones;
my %starknownones;
my %seen;
if($options{'k'}){
    open KN,"<$options{'k'}" or die "Error: $options{'k'} file not found\n";
    my $seq1;
    while(<KN>){
        if(/\>(\S+)/){
            $knownones{$1};
            $seq1= lc <KN>;
			chomp $seq1;
            $seq1=~ tr /t/u/;
            $knownones{$1}=$seq1;
        }
    }
    close KN;
}
 
## create csv file    
if($options{'c'}){
	$csv = 1;
    open CSV,">$cwd/result_${time}.csv" or die "Error: cannot create csv file\n";
}


my $reads=0;
my @line;

## open output files
open HTML,">$cwd/result_${time}.html" or die "Error: cannot create html file\n";

## Create HTML file
CreateHTML();

## print paramters and files used to html file
ReadInParameters();

my $spacer;   ## length of longest entry
my $spaces;   ## string of spaces to fill up spacer


## Specificity value hash for each score cutoff value
my %SP_values;


## read in survey file
if($options{'s'}){
    print "making survey now\n";
    $spacer = length("known hairpins (estimated false positives):       ");
    $spaces = 0;
    
    ## is for extended output
    if($options{'e'}){

        my $p1 ='<th><a href="#" class="tooltip">';
        my $p11 ='<th><a href="#" class="tooltip2">';
        #my $p2='<span class="tooltip"><span class="top"></span><span class="middle">';
        #my $q ='</span><span class="bottom"></span></span></a></th>';

        my $p2='<span>';
        my $q ='</span></a></th>';

        ## hash for mouse overs
        my %ha;
		$ha{1}{1} = 'miRDeep2 score';
$ha{1}{2} = 'for details on how the log-odds score is calculated, see Friedlander et al., Nature Biotechnology, 2008.';
$ha{2}{1} = 'predicted by miRDeep2';
$ha{2}{2} = 'novel miRNA hairpins are here defined by not having any of the reference mature miRNAs mapping perfectly (full length, no mismatches). The numbers show how many novel miRNA hairpins have a score equal to or exceeding the cut-off.';
$ha{3}{1} = 'estimated false positives';
$ha{3}{2} = 'number of false positive miRNA hairpins predicted at this cut-off, as estimated by the miRDeep2 controls (see Friedlander et al., Nature Biotechnology, 2008). Mean and standard deviation is estimated from 100 rounds of permuted controls.';
$ha{4}{1} = 'estimated true positives';
$ha{4}{2} = 'the number of true positive miRNA hairpins is estimated as t = total novel miRNAs - false positive novel miRNAs. The percentage of the predicted novel miRNAs that is estimated to be true positives is calculated as p = t / total novel miRNAs. The number of false positives is estimated from 100 rounds of permuted controls. In each of the 100 rounds, t and p are calculated, generating mean and standard deviation of t and p. The variable p can be used as an estimation of miRDeep2 positive predictive value at the score cut-off. ';
$ha{5}{1} = 'in species';
$ha{5}{2} = 'number of reference mature miRNAs for that species given as input to miRDeep2.';
$ha{6}{1} = 'in data';
$ha{6}{2} = 'number of reference mature miRNAs for that species that map perfectly (full length, no mismatches) to one or more of precursor candidates that have been excised from the genome by miRDeep2.';
$ha{7}{1} = 'detected by miRDeep2';
$ha{7}{2} = 'number of reference mature miRNAs for that species that map perfectly (full length, no mismatches) to one or more of predicted miRNA hairpins that have a score equal to or exceeding the cut-off. The percentage of reference mature miRNAs in data that is detected by miRDeep2 is calculated as s = reference mature miRNAs detected / reference mature miRNAs in data. s can be used as an estimation of miRDeep2 sensitivity at the score cut-off.';
$ha{8}{1} = 'estimated signal-to-noise';
$ha{8}{2} = 'for the given score cut-off, the signal-to-noise ratio is estimated as r = total miRNA hairpins reported / mean estimated false positive miRNA hairpins over 100 rounds of permuted controls.';
$ha{9}{1} = 'excision gearing';
$ha{9}{2} = 'this is the minimum read stack height required for excising a potential miRNA precursor from the genome in this analysis.';


## hash finished

        

        open SURVEY,"<$options{'s'}" or die "Error: $options{'s'} file not found $!\n";
        my @tmp = split(/\t/,<SURVEY>);
        my $reduced = 0;
        if(scalar @tmp < 4){
            $reduced =1;
        }
        close SURVEY;

        open SURVEY,"<$options{'s'}" or die "Error: $options{'s'} file not found $!\n";
        print HTML "
<font face=\"Times New Roman\" size=\"4\">
<b>Survey of miRDeep2 performance for score cut-offs -10 to 10</b>\n
<font face=\"Times New Roman\" size=\"3\">
<table border=\"1\" >\n";
        my $i = 0;
        my %sh;


        if(not $reduced){
            
            print HTML "<tr><th></th><TH colspan=\"3\">novel miRNAs</th><th colspan=\"3\">known miRBase miRNAs</th><th></th><th></th></tr>\n";


            for(sort {$a <=> $b} keys %ha){
                if($_ eq 1){
                    print HTML "<th><a href=\"http://www.nature.com/nbt/journal/v26/n4/abs/nbt1394.html\" target=\"blank\" class=\"tooltip\">$ha{$_}{1}$p2$ha{$_}{2}$q</a>\n\n";

                
                }
                elsif($_ ne 9){
                    print HTML "$p1$ha{$_}{1}$p2$ha{$_}{2}$q\n\n";
                }else{
                    print HTML "$p11$ha{$_}{1}$p2$ha{$_}{2}$q\n\n";
                }
            }
        }else{
            print HTML "<th><a href=\"http://www.nature.com/nbt/journal/v26/n4/abs/nbt1394.html\" target=\"blank\" class=\"tooltip\">$ha{1}{1}$p2$ha{1}{2}$q</a>\n\n";
            print HTML "$p1$ha{8}{1}$p2$ha{8}{2}$q\n\n";
            print HTML "$p1$ha{9}{1}$p2$ha{9}{2}$q\n\n";

        }




		
        while(<SURVEY>){
			my @line = split(/\t/);
            if($csv and not $i){
                print CSV "$_";
            }
			
			$i++;
            if($i>1){
                if($csv){
                    print CSV "$_";
                }
                $_ =~ s/\t/<td>/g;
				if($line[0] < 0){
					$sh{$i} = "<tr bgcolor=silver><td>$_</td></tr>\n";
				}else{
					$sh{$i} = "<tr><td>$_</td></tr>\n";
				}
            }
        }
        if($csv){
            print CSV "\n\n\n";
        }
        
        ##print in reverse order
        for(sort {$a <=> $b } keys %sh){
            $sh{$_} =~ s/\+\/-/\&\#177/g;
            print HTML $sh{$_};
        }
		
        print HTML "</table><br><br>\n";
        close SURVEY; 
    }
    
    
    open SURVEY,"<$options{'s'}" or die "Error: $options{'s'} file not found $!\n";
    my @header=split(/\t/,<SURVEY>);
    my @survey;
    my @std_survey;
	
    ## also report survey for score cutoff
    if($options{'g'}){
        print HTML "
<font face=\"Times New Roman\" size=\"3\">
<table border=\"1\" >\n
    <th>Survey of run</th>\n";
    }
	
    while(<SURVEY>){
        @survey = split(/\t/);
        my $lcol;  
        if($survey[3] =~ /\((.+)\)/){
            $lcol = $1;
            $lcol =~ s/\+\/-/\&\#177/g;
        }
        
#         ## get some numbers of miRBase miRNAs
#         $aknown = $survey[4];
#         my @trx = split(/\s+/,$survey[4]) ;
#         if($survey[3] - $trx[0] < $anotscored){
#             $anotscored = $survey[3] - $trx[0];
#             $aknown = $survey[3];
#             $anovel = $survey[9];    
#         }
        

       
        #if($lcol>= 75 and $threshold eq "na" ){
        #    $threshold = $survey[0];
        #}
		
        if($survey[0] eq 0){
            @std_survey = @survey;
        }

        $SP_values{$survey[0]} = $lcol;


        ## if line in survey file is reached that hits the threshold then print out stuff
        if($survey[0] eq $threshold){
            for(my $i = 0; $i < scalar @survey; $i++){

                    $spaces = $spacer - (length($header[$i])+1);
                    if($i eq 0){
						#print CSV "Survey of run\t\n$header[$i] used\t$survey[$i]\n" if($csv);
                        if($options{'g'}){
                            print HTML "<tr><td>$header[$i] used</td><td>$survey[$i]</td></tr>\n"; 
                        }
                    }else{
						#print CSV "$header[$i]\t$survey[$i]\n" if($csv);
                        if($options{'g'}){
                            print HTML "<tr><td>$header[$i]</td><td>$survey[$i]</td></tr>\n"; 
                        }
                    }

            }

            my $expected = $SP_values{$survey[0]};
            my $novel;
            if($survey[6] =~ /(\d+)\(\d+/){
                $novel = int($1+0.5);
            }
			#print CSV "novel hairpins, expected true\t$novel ($expected%)\n" if($csv);
            if($options{'g'}){
                print HTML "<tr><td>novel hairpins, expected true</td><td>$novel ($expected%)</td></tr>\n"; 
            }
        }
    }
    close SURVEY;

    if($threshold eq "na"){
        $threshold = $std_survey[0];
        for(my $i = 0; $i < scalar @std_survey; $i++){
            #if($sh{$i}){
                $spaces = $spacer - (length($header[$i])+1);
                if($options{'g'}){
                    print HTML "<tr><td>$header[$i]</td><td>$std_survey[$i]</td></tr>\n"; 
                }
            #}
        }
    }
    if($options{'g'}){
        print HTML "</table></font>\n";
    }
}else{
    die "Error: no survey file specified with option -s. Please also make sure that your survey file is not empty\n\n";
}



#########################################################
##
##
## here the output.mrd file given by switch -f is parsed
##
##
#########################################################

open IN,"<$infile" or die "Error: cannot open $infile\n";
my $oid;
my $star_exp_hit = 0;
my %star_exp_hit_pos;

while(<IN>){
	if(/^\>\s*(\S+)/){
		$id = $1;
        $oid = $1;
        $id =~ tr/|/_/;
        $hash{$id}{"oid"}=$oid;
		$hash{$id}{"id"}=$id;
        $counter =0;
	}
	elsif(/^score_mfe\s+(\S+)/){
		$hash{$id}{"score_mfe"}=$1;
	}
	elsif(/^score for star read.+\s+(\S+)/){
		$hash{$id}{"score_star"}=$1;
	}
	elsif(/^score for read counts\s+(\S+)/){
		$hash{$id}{"score_read"}=$1;
	}
	elsif(/^score for mfe\s+(\S+)/){
		$hash{$id}{"score_mfe"}=$1;
	}
	elsif(/^score total\s+(\S+)/){
		$hash{$id}{"score"}=$1;
	}
	elsif(/^total read count\s+(\S+)/){
		$hash{$id}{"freq_total"}=$1;
	}
	elsif(/^mature read count\s+(\S+)/){
		$hash{$id}{"freq_mature"}=$1;
	}
	elsif(/^loop read count\s+(\S+)/){
		$hash{$id}{"freq_loop"}=$1;
	}
	elsif(/^star read count\s+(\S+)/){
		$hash{$id}{"freq_star"}=$1;
	}
	elsif(/^pri_seq\s+(\S+)/){
		$hash{$id}{"pri_seq"}=$1;
        my @d;
        if($hash{$id}{'obs'}){
            @d = split(//,$hash{$id}{'obs'});
        }else{
            @d = split(//,$hash{$id}{'exp'});
        }
        my @s = split(//,$hash{$id}{"pri_seq"});
        my $mseq="";
        my $sseq="";
        my $lseq='';

        for(my $i=0; $i < length($1); $i++){
            if($d[$i] ne "f"){ $hash{$id}{"ucsc_seq"} .= $s[$i];}
            if($d[$i] eq "M"){
                $mseq.= $s[$i];
            }elsif($d[$i] eq "S"){
                $sseq.= $s[$i];
            }elsif($d[$i] eq "l"){
                $lseq.= $s[$i];
            }else{}
        }
        my $sseq_obs="";
        if($hash{$id}{'obs'}){
            @d = split(//,$hash{$id}{'obs'});
            for(my $i=0; $i < length($1); $i++){
                $sseq_obs.= $s[$i] if($d[$i] eq "S");
            }
        }


        $hash{$id}{'mat_seq'} = $mseq; 
        $hash{$id}{'loop_seq'} = $lseq;
        $hash{$id}{"star_seq"} = $sseq;
        $hash{$id}{"star_seq_obs"} = $sseq_obs;
	}
    elsif(/^score for randfold\s+(\S+)/){
        $hash{$id}{"rand"}=$1;
        if($1>0){
            $hash{$id}{"randfold"}="yes";
        }else{
            $hash{$id}{"randfold"}="no";
        }
	}
    elsif(/^score for cons.\s+seed\s+(\S+)/){
        $hash{$id}{"score_cons"}=$1;
        if($1>0){
            $hash{$id}{"cons_seed"}="yes";
        }else{
            $hash{$id}{"cons_seed"}="";
        }
	}
    elsif(/^miRNA with same seed\s+(\S+)/){
        $hash{$id}{"cons_seed"}=$1;
    }

	elsif(/^exp\s+(\S+)/){
		$hash{$id}{'exp'}=$1;
	}
	elsif(/^obs\s+(\S+)/){
		$hash{$id}{'obs'}=$1;
	}
	elsif(/^pri_struct\s+(\S+)/){
		$hash{$id}{"pri_struct"}=$1;
		$reads=1;
        $counter = 0;

		next;
	}
	elsif(/^(\S\S\S)(\S+)\s+(\S+)\s+(\S+)$/ and $reads){  ## do not read in known miRNAs
	    $counter++;
	    $hash{$id}{"reads"}{$1}{$counter}{"rid"} = "$1$2";
	    $hash{$id}{"reads"}{$1}{$counter}{"seq"} = $3;
	    $hash{$id}{"reads"}{$1}{$counter}{"mm"}  = $4;
	    
	    if($knownones{"$1$2"}){
			$seen{"$1$2"} = 1;
	    }


	    ## check here if expected sequence is hit partially by at least one read ## maybe full coverage should be check => later
	    if($star_exp_hit eq 0 and $hash{$id}{"reads"}{$1}{$counter}{"seq"} =~ /$hash{$id}{"star_seq"}/i){
			$star_exp_hit = 1;
	    }

        
        #open MMU,">>recovered_miRNA";
        if($knownones{"$1$2"} and $reads){
        #    print MMU "$1\n";
            if(not $hash{$id}{"known"}){
				#die "$id $1$2";
                $hash{$id}{"known"} = "$1$2";
		
		## now check if mature of miRDeep is same as mirbase or not
		if(index($knownones{"$1$2"},substr($hash{$id}{'mat_seq'},4,12))>=0){
		    $hash{$id}{'mirbase'} = 'TRUE';
			$hash{$id}{'mirbase-out'} = "$1$2";
			## if star sequence mirbase mature then output star
		}elsif(index($knownones{"$1$2"},substr($hash{$id}{'star_seq'},4,12))>=0){
		    $hash{$id}{'mirbase'} = 'STAR';
			$hash{$id}{'mirbase-out'} = "$1$2";

			## if mature matches mirbase star then output star sequence, else output NA and precursor would be a novel prediction
			# if(index($starknownones{"$1$2"},substr($hash{$id}{'mat_seq'},4,12))>=0){
			# 	$hash{$id}{'mirbase-out'} = "$1$2";
			# }else{
			# 	$hash{$id}{'mirbase-out'} = "NA";
			# }
		}else{## this case should not happen
		    $hash{$id}{'mirbase'} = 'NA';
			$hash{$id}{'mirbase-out'} = "NA";
		}

            }
            next;
        }
        #close MMU;
        
        if(/^\s+$/ and $reads){
            $reads = 0;
            next;
        }
    }else{}
}
close IN;
print STDERR "parsing input file finished\n";

##########################################################
#
#
# Check precursors for Rfam hits
#
#
##########################################################

if($options{'r'}){
    print STDERR "checking Rfam for hits to precursors\n";
    check_Rfam(%hash);
}else{
    print STDERR "no Rfam file specified\tskipping\n";
}


my $percentage=0;
my @klist = sort {$b <=> $a} keys %SP_values;
my $klen = scalar @klist;

my $blat;


## count number of novel miRNAs
my $novelc = 0;
for(keys %hash){
    next if($hash{$_}{"known"});
    $novelc++;
}


##################################################################
#
#
# Print to HTML file the novel miRNAs
#
#
##################################################################

if(not $novelc){
    print HTML " <br><h2>no novel miRNAs detected</h2><br>";
}else{
    PrintHtmlTableHeader("novel miRNAs predicted by miRDeep2");

    if($csv){
        PrintHtmlTableHeader("novel miRNAs predicted by miRDeep2",1);
    }

    for $id (sort { $hash{$b}{"score"} <=> $hash{$a}{"score"} } keys %hash){
        
        next if($hash{$id}{"known"});
        if($options{'b'}){
            next if($confident{$id} != 1);
        }

        if($hash{$id}{"score"} >= $threshold){
            

            $hash{$id}{"pdf"} = 1;

            if($hash{$id}{"score"} >= $klist[0]){
                $percentage = $SP_values{$klist[0]};
            }else{
                for(my $i = 0; $i < $klen; $i++){
                    if($hash{$id}{"score"} >= $klist[$i]){
                        $percentage = $SP_values{$klist[$i]};
                        last;
                    }
                }
            }

            
            ## print to CSV file
            if($csv){
                my $s_star=$hash{$id}{'star_seq'};
                $s_star=$hash{$id}{'star_seq_obs'} if($hash{$id}{'star_seq_obs'});
                my $p =$percentage;
                $p =~ s/&\#177/+\/-/;
                
                print CSV "$hash{$id}{'oid'}\t$hash{$id}{'score'}\t$p\t";

                if($hash{$id}{'rfam'}){print CSV "$hash{$id}{'rfam'}\t";}else{print CSV "-\t";}
                
                print CSV "$hash{$id}{'freq_total'}\t$hash{$id}{'freq_mature'}\t$hash{$id}{'freq_loop'}\t$hash{$id}{'freq_star'}\t";
                if($hash{$id}{'randfold'}){print CSV "$hash{$id}{'randfold'}\t-\t";} else{ print CSV "-\t-\t";}
                if($hash{$id}{'cons_seed'}){print CSV "$hash{$id}{'cons_seed'}\t-\t-\t";} else{ print CSV "-\t-\t-\t";}
                print CSV "$hash{$id}{'mat_seq'}\t$s_star\t$hash{$id}{'ucsc_seq'}";

				# if($options{'p'}){
				# 	my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
				# 	if($pres_coords{$id}{'strand'} eq '-'){
				# 		print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}";
				# 	}else{
				# 		print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset,"..",$pres_coords{$id}{'s'}+$offset+length($hash{$id}{'ucsc_seq'}),":$pres_coords{$id}{'strand'}";
				# 	}



				# 	print CSV "\t$pres_coords{$id}{'chr'}:$pres_coords{$id}{'s'}..$pres_coords{$id}{'e'}:$pres_coords{$id}{'strand'}";
				# }

				
				if($options{'p'}){
					my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
					if($pres_coords{$id}{'strand'} eq '-'){
						print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}";
					}else{
						print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset-1,"..",$pres_coords{$id}{'s'}+$offset+length($hash{$id}{'ucsc_seq'})-1,":$pres_coords{$id}{'strand'}";
					}
				}

				print CSV "\n";
            }
            
            if($org eq ""){
                $blat = '<td></td>';
            }else{
                $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash{$id}{'ucsc_seq'}\" target=\"_blank\">blat</a></td>";
            }
            my $s_star=$hash{$id}{'star_seq'};
            $s_star=$hash{$id}{'star_seq_obs'} if($hash{$id}{'star_seq_obs'});
            

            $known="<td nowrap=\"nowrap\"></td>";
            
            #if($hash{$id}{"score"} > 10 ){ $percentage = ">$percentage";}
            
            $n_count++;

            my $science= getEvalue($hash{$id}{"score"},1);
	#		print STDERR $percentage,"\n";
			if($percentage =~ /^(\d+)\s+(\S+)\s+(\d+)\%$/){
		#		print STDERR "$1\t$2\t$3\n";
				my $one = $1/100;
				my $mid=$2;
				my $two = $3/100;
#				die "$3\t$two\n";
				if(length($one) < 4){
					$one .= "0" x (4-length($one));
				}
				if($one =~ /0000/){$one = '0.00';}
				
				if(length($two) < 4){
					$two .= "0" x (4-length($two));
				}
				if($two =~ /0000/){$two = '0.00';}
				
				$percentage = "$one $mid $two";
			}
			#die $percentage,"\n";
        

			#print STDERR "==$percentage\n";

            print HTML <<EOF;
            <tr><td><a href="pdfs_$time/$id.pdf">$hash{$id}{'oid'}</a></td>
                <td nowrap="nowrap">$science</td>
                <td>$percentage</td>
                <td>$hash{$id}{'rfam'}</td>
                <td>$hash{$id}{"freq_total"}</td>
                <td>$hash{$id}{"freq_mature"}</td>
                <td>$hash{$id}{"freq_loop"}</td>
                <td>$hash{$id}{"freq_star"}</td>
                <td>$hash{$id}{"randfold"}</td>
                $known
                <td nowrap="nowrap">$hash{$id}{"cons_seed"}</td>
                $blat
                <td><a href=${blast}$hash{$id}{"ucsc_seq"}&JOB_TITLE=$hash{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
                <td>$hash{$id}{"mat_seq"}</td>
                <td>$s_star</td>
                <td>$hash{$id}{"ucsc_seq"}</td>                
EOF


if($options{'p'}){
	my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
	if($pres_coords{$id}{'strand'} eq '-'){
		print HTML "<td>$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}</td>\n";
	}else{
		print HTML "<td>$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset-1,"..",$pres_coords{$id}{'s'}+$offset+length($hash{$id}{'ucsc_seq'})-1,":$pres_coords{$id}{'strand'}</td>\n";
	}
}

# if($options{'a'}){
#     print HTML "<td>$mature_pos_hash{$id}{'s'}-$mature_pos_hash{$id}{'e'}:$mature_pos_hash{$id}{'strand'}</td>\n";
# }
print HTML  "</tr>";

#&Lucky=\"I'm feeling lucky\"
#            <td>$hash{$id}{"pri_seq"}</td>
                
            }
    }

    ## now print the unconfident ones
    if($options{'b'}){
        print HTML  "<tr><td>--</td></tr>";
    
    
		for $id (sort { $hash{$b}{"score"} <=> $hash{$a}{"score"} } keys %hash){
			
			next if($hash{$id}{"known"});
			if($options{'b'}){
				next if($confident{$id} == 1);
			}
			
			if($hash{$id}{"score"} >= $threshold){
				
				
				$hash{$id}{"pdf"} = 1;
				
				if($hash{$id}{"score"} >= $klist[0]){
					$percentage = $SP_values{$klist[0]};
				}else{
					for(my $i = 0; $i < $klen; $i++){
						if($hash{$id}{"score"} >= $klist[$i]){
							$percentage = $SP_values{$klist[$i]};
							last;
						}
					}
				}
				
				
            ## print to CSV file
				if($csv){
					my $s_star=$hash{$id}{'star_seq'};
					$s_star=$hash{$id}{'star_seq_obs'} if($hash{$id}{'star_seq_obs'});
					my $p =$percentage;
					$p =~ s/&\#177/+\/-/;
					
					print CSV "$hash{$id}{'oid'}\t$hash{$id}{'score'}\t$p\t";
					
					if($hash{$id}{'rfam'}){print CSV "$hash{$id}{'rfam'}\t";}else{print CSV "-\t";}
					
					print CSV "$hash{$id}{'freq_total'}\t$hash{$id}{'freq_mature'}\t$hash{$id}{'freq_loop'}\t$hash{$id}{'freq_star'}\t";
					if($hash{$id}{'randfold'}){print CSV "$hash{$id}{'randfold'}\t-\t";} else{ print CSV "-\t-\t";}
					if($hash{$id}{'cons_seed'}){print CSV "$hash{$id}{'cons_seed'}\t-\t-\t";} else{ print CSV "-\t-\t-\t";}
					print CSV "$hash{$id}{'mat_seq'}\t$s_star\t$hash{$id}{'ucsc_seq'}";

					# if($options{'p'}){
	# 					my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
	# 					if($pres_coords{$id}{'strand'} eq '-'){
	# 						print CSV "<td>$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}</td>\n";
	# 					}else{
	# 						print CSV "<td>$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset,"..",$pres_coords{$id}{'s'}+$offset+length($hash{$id}{'ucsc_seq'}),":$pres_coords{$id}{'strand'}</td>\n";
	# }

	# 				}else{
	# 					print CSV "\tna";
	# 				}


					if($options{'p'}){
						my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
						if($pres_coords{$id}{'strand'} eq '-'){
							print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}";
						}else{
							print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset,"..",$pres_coords{$id}{'s'}+$offset+length($hash{$id}{'ucsc_seq'}),":$pres_coords{$id}{'strand'}";
						}
					}




					print CSV "\n";
					
				}
				
				if($org eq ""){
					$blat = '<td></td>';
				}else{
					$blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash{$id}{'ucsc_seq'}\" target=\"_blank\">blat</a></td>";
				}
				my $s_star=$hash{$id}{'star_seq'};
				$s_star=$hash{$id}{'star_seq_obs'} if($hash{$id}{'star_seq_obs'});
				
				
				$known="<td nowrap=\"nowrap\"></td>";
				
				#if($hash{$id}{"score"} > 10 ){ $percentage = ">$percentage";}
				
				$n_count++;
				
				my $science= getEvalue($hash{$id}{"score"},1);
			
				#print STDERR "==$percentage\n";
				
				if($percentage =~ /^(\d+)\s+(\S+)\s+(\d+)\%$/){
					my $one = $1/100;
					my $mid=$2;
					my $two = $3/100;
					if(length($one) < 4){
					$one .= "0" x (4-length($one));
				}
					if($one =~ /0000/){$one = '0.00';}
					
					if(length($two) < 4){
						$two .= "0" x (4-length($two));
					}
					if($two =~ /0000/){$two = '0.00';}
					
					$percentage = "$one $mid $two";
				}
				
				print HTML <<EOF;
				<tr><td><a href="pdfs_$time/$id.pdf">$hash{$id}{'oid'}</a></td>
                <td nowrap="nowrap">$science</td>
                <td>$percentage</td>
                <td>$hash{$id}{'rfam'}</td>
                <td>$hash{$id}{"freq_total"}</td>
                <td>$hash{$id}{"freq_mature"}</td>
                <td>$hash{$id}{"freq_loop"}</td>
                <td>$hash{$id}{"freq_star"}</td>
                <td>$hash{$id}{"randfold"}</td>
                $known
                <td nowrap="nowrap">$hash{$id}{"cons_seed"}</td>
                $blat
                <td><a href=${blast}$hash{$id}{"ucsc_seq"}&JOB_TITLE=$hash{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
                <td>$hash{$id}{"mat_seq"}</td>
                <td>$s_star</td>
                <td>$hash{$id}{"ucsc_seq"}</td>                
EOF
				
if($options{'p'}){
    print HTML "<td>$pres_coords{$id}{'chr'}:$pres_coords{$id}{'s'}..$pres_coords{$id}{'e'}:$pres_coords{$id}{'strand'}</td>\n";
}else{
	print HTML "<td>na</td>\n";
}


print HTML  "</tr>";

            }
    }
    } ## again end 2

    print HTML "</table></font>";
}


##################################################################
#
#
# Print to HTML file the known miRNAs
#
#
##################################################################

if($options{'k'}){
    PrintHtmlTableHeader("mature miRBase miRNAs detected by miRDeep2");


    if($csv){
        print CSV "\n\n\n";
        PrintHtmlTableHeader("mature miRBase miRNAs detected by miRDeep2",1);
    }
    
    
    for $id (sort { $hash{$b}{"score"} <=> $hash{$a}{"score"} } keys %hash){
        next if(not $hash{$id}{"known"});

		$hash{$id}{"pdf"} = 1;

		if($hash{$id}{"score"} >= $klist[0]){
			$percentage = $SP_values{$klist[0]};
		}else{
            for(my $i = 0; $i < $klen; $i++){
                if($hash{$id}{"score"} >= $klist[$i]){
                    $percentage = $SP_values{$klist[$i]};
                    last;
                }
            }
		}
		
        
		
		if($hash{$id}{"known"}){
			$known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$hash{$id}{\"known\"}\" target=\"_blank\"> $hash{$id}{\"known\"}</a></td>";
		}else{
			$known="<td nowrap=\"nowrap\"></td>";
		}

	 	if($csv){
            my $s_star=$hash{$id}{'star_seq'};
            $s_star=$hash{$id}{'star_seq_obs'} if($hash{$id}{'star_seq_obs'});
            my $p =$percentage;
            $p =~ s/&\#177/+\/-/;
            
            if(not $novelc){
                print CSV "$hash{$id}{'oid'}\t$hash{$id}{'score'}\t-\t";
            }else{
                print CSV "$hash{$id}{'oid'}\t$hash{$id}{'score'}\t$p\t";
            }
            
            if($hash{$id}{'rfam'}){print CSV "$hash{$id}{'rfam'}\t";}else{print CSV "-\t";}
            
            print CSV "$hash{$id}{'freq_total'}\t$hash{$id}{'freq_mature'}\t$hash{$id}{'freq_loop'}\t$hash{$id}{'freq_star'}\t";
            if($hash{$id}{'randfold'}){print CSV $hash{$id}{'randfold'},"\t";} else{ print CSV "-\t";}
            if($hash{$id}{'known'}){print CSV "$hash{$id}{'known'}\t";} else{ print CSV "-\t";}
            if($hash{$id}{'cons_seed'}){print CSV "$hash{$id}{'cons_seed'}\t-\t-\t";} else{ print CSV "-\t-\t-\t";}
            print CSV "$hash{$id}{'mat_seq'}\t$s_star\t$hash{$id}{'ucsc_seq'}";
			
			# if($options{'p'}){
			# 	print CSV "\t$pres_coords{$id}{'chr'}:$pres_coords{$id}{'s'}..$pres_coords{$id}{'e'}:$pres_coords{$id}{'strand'}";
			# }else{
			# 	print CSV "\tna";
			# }

			
if($options{'p'}){
	my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
	if($pres_coords{$id}{'strand'} eq '-'){
		print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}";
	}else{
		print CSV "\t$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset-1,"..",$pres_coords{$id}{'s'}+$offset-1+length($hash{$id}{'ucsc_seq'}),":$pres_coords{$id}{'strand'}";
	}
}




			print CSV "\n";
		}
		
        if($org eq ""){
            $blat = '<td></td>';
        }else{
            $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash{$id}{'ucsc_seq'}\"  target=\"_blank\">blat</a></td>";
        }
        my $s_star=$hash{$id}{'star_seq'};
        $s_star=$hash{$id}{'star_seq_obs'} if($hash{$id}{'star_seq_obs'});
		
        $k_count++;
        
        my $science= getEvalue($hash{$id}{"score"},1);
        
        if(not $novelc){
            $percentage = "-";
        }
        
        
        if($percentage =~ /^(\d+)\s+(.+)\s+(\d+)\%$/){
            my $one = $1/100;
            my $mid=$2;
            my $two = $3/100;
             if(length($one) < 4){
                 $one .= "0" x (4-length($one));
             }
             if($one =~ /0000/){$one = '0.00';}
            
             if(length($two) < 4){
                 $two .= "0" x (4-length($two));
             }
             if($two =~ /0000/){$two = '0.00';}

            $percentage = "$one $mid $two";
        }
        
        
        print HTML <<EOF;
		<tr><td><a href="pdfs_$time/$id.pdf">$hash{$id}{'oid'}</a></td>
			<td nowrap="nowrap">$science</td>
            <td>$percentage</td>
            <td>$hash{$id}{'rfam'}</td>
	    <td>$hash{$id}{'mirbase'}</td>
			<td>$hash{$id}{"freq_total"}</td>
			<td>$hash{$id}{"freq_mature"}</td>
            <td>$hash{$id}{"freq_loop"}</td>
            <td>$hash{$id}{"freq_star"}</td>
            <td>$hash{$id}{"randfold"}</td>
			$known
			<td nowrap="nowrap">$hash{$id}{"cons_seed"}</td>
            $blat
            <td><a href=${blast}$hash{$id}{"ucsc_seq"}&JOB_TITLE=$hash{$id}{'oid'}${blast_query}  target="_blank">blast</a></td>
            <td>$hash{$id}{"mat_seq"}</td>
            <td>$s_star</td>
            <td>$hash{$id}{"ucsc_seq"}</td>
EOF
## <td nowrap="nowrap">$hash{$id}{"cons_seed"}</td>

if($options{'p'}){
	my $offset=index($hash{$id}{'pri_seq'},$hash{$id}{'ucsc_seq'});
	if($pres_coords{$id}{'strand'} eq '-'){
		print HTML "<td>$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'e'}-$offset-length($hash{$id}{'ucsc_seq'}),"..",$pres_coords{$id}{'e'}-$offset,":$pres_coords{$id}{'strand'}</td>\n";
	}else{
		print HTML "<td>$pres_coords{$id}{'chr'}:",$pres_coords{$id}{'s'}+$offset-1,"..",$pres_coords{$id}{'s'}-1+$offset+length($hash{$id}{'ucsc_seq'}),":$pres_coords{$id}{'strand'}</td>\n";
	}
}else{
	print HTML "<td>na</td>";
}

	print HTML  "</tr>";

    }
    close CSV if($csv);
    close IN;
    print HTML "</table></font>";
}

##printing to HTML file finished for novel and known miRNAs


##################################################################
#
#
# Print to HTML file the miRNAs in data but not detected by miRDeep2
#
#
##################################################################

if($options{'q'}){
#    PrintHtmlTableHeader("mature miRBase miRNAs not detected by miRDeep2");
    PrintKnownnotfound();

}else{
    close CSV if($csv);
    close IN;
}
CloseHTML();




##################################################################
#
#
# Create PDF files for all entrys in the HTML file
#
#
##################################################################

if(not $options{'d'}){
    CreateStructurePDF(%hash);
    if($options{'q'}){
        $mirbase = 1;
        CreateStructurePDF(%hash_q);
    }
}

# print "


# novel ones printed:\t $n_count

# known ones printed:\t $k_count

# not scored printed:\t $sig

# star_c > mature_c: \t $e_count


# ";

exit;



###########################################################################
####
#### Subroutines
####
###########################################################################


sub CreateStructurePDF{
    my %hash = @_;
    my $filename;
    print STDERR "creating PDF files\n";
    for(sort { $hash{$b}{"score"} <=> $hash{$a}{"score"} } keys %hash){

        next if(not $hash{$_}{'pdf'}); 
		next if(not $hash{$_}{'freq_total'});
        $sid = $_;
        $sid =~ tr/\|/_/;
        %star_exp_hit_pos =();
        $filename = $sid;

        if($mirbase){
            $filename = $sid; #$hairpin2mature{$sid};
        }

        next if ($seen{$filename});
		next if(-f "$cwd/pdfs_$time/$filename.pdf");
#		next if($hash{$sid}{"score"} < $threshold); ## skip if threshold is not reached not used anymore
        ## reinit variables;
        $i=0;



        $offset = 0;

        $me=0;     ## mature end coordinate
        @desc;
        
        $lflank1 = 0; ## length(string of left flank)
        $fl1 = 0;    ## 
        $lflank2 = 0; ## length string of right flank
        $fl2b=-1;   ## right flank begin  
        $lloop = 0;   ## string of loop
        $lb=-1;     ## starting position of loop
        $lstar = 0;   ## string of star sequence
        $sb=-1;     ## starting 
        $lmature = 0; ## string of mature 
        $mb=-1;     ## mature begin
        $struct = 0; ## structure string
        $pri_seq="";## pri-cursor sequence
        $lenstr=0; 

        $pdf;    ## pdf descriptor
        $page;   ## page descriptor
        $gfx;    ## graphic variable
        $trb;    ## fontvariable

        %hash2 = ();
        %hash2c = ();
        %hash2mm = ();
        %hash2order = ();
        %order = ();

        $yorig = 500;
        $downy = 50;

        $dline;                            ## line graphic handler

        $first=1;
        $lastx=0;
        $lasty=0;

        $final="";                         ## final output string of a read
        @pseq;                             ## precursor sequence  
        @rseq;                             ## read sequence

        $totalreads = 0;

        %assign_str = ();
        %assign_str_exp = ();

        %struct = ();
        

        $bpo1=-10;                             ## left nt pos in first bp 
        $bpo2=-10;                             ## right nt pos in first bp 
        $bpo1r=-10;                            ## left nt pos in second bp 
        $bpo2r=-10;                            ## right nt pos in second bp 


        $ffe=0;                                ## first flank end position
        $ff2b=0;                               ## second flank begin position
        
        @sorted;                               ## array that stores sorted order of fl1,m,l,s,fl2 
        $y=$yorig;                             ## y coordinate


        ($minx,$miny,$maxx,$maxy);             ## min and max x,y coordinates of rna sequence
        @rna;                                  ## rna sequence

        %xc = ();
        %yc = ();


#### main program;       
        $pri_seq =  $hash{$sid}{"pri_seq"};
        @rna = split(//,$pri_seq);
        
        chomp $pri_seq;

        if($hash{$sid}{'obs'}){
            my @desc2 =split(//,$hash{$sid}{'exp'});
            for (my $i=0;$i < scalar @desc2; $i++){
                if ($desc2[$i] eq "f"){   ## assign_str now starts at 0 not at one
                    $assign_str_exp{$i} = "black";
                }elsif ($desc2[$i] eq "M"){
                    $assign_str_exp{$i} = $col_mature;
                } elsif ($desc2[$i] eq "l"){
                    $assign_str_exp{$i} = $col_loop;
                } elsif ($desc2[$i] eq "S"){
                    $assign_str_exp{$i} = $col_star_exp;
                } else {
                    print STDERR "something went wrong while parsing alignment in output.mrd file $!\n";
                }   
            }
        }

        ## if observed star strand then use this one otherwise use the expected one
        my @desc_tmp;
        if($hash{$sid}{'obs'}){
            @desc = split(//,$hash{$sid}{'obs'});
            @desc_tmp = split(//,$hash{$sid}{'exp'});
        }else{
            @desc = split(//,$hash{$sid}{'exp'});
        }
 
        my $run_var_loop =0;
        for (my $i=0;$i < scalar @desc; $i++){
            #print STDERR "-> $i\t$desc_tmp[$i]\t$desc[$i]\t$star_exp_hit\n";

            if ($desc[$i] eq "f"){   ## assign_str now starts at 0 not at one
                    $assign_str{$i} = "black";

                if ($mb eq -1){           ## if in first flank
                    $ffe = $i if($desc[$i+1] ne "f");
                    $lflank1++;

                }
                if ($mb ne -1){           
                    

                    ## second flank   ## dependent on exp/obs both or not
                    if($star_exp_hit and $desc_tmp[$i] eq "f"){
                        $fl2b = $i if($fl2b eq -1);
                        $order{"f"} = $fl2b if(not $order{"f"});
                        $lflank2++;
                    }elsif($star_exp_hit and $desc_tmp[$i] eq "S"){
                        $star_exp_hit_pos{$i}=$i;                            ##star hit expression
                        next;
                    }else{
                        $fl2b = $i if($fl2b eq -1);
                        $order{"f"} = $fl2b if(not $order{"f"});
                        $lflank2++;
                    }
                }
            }elsif ($desc[$i] eq "M"){
                $mb = $i if($mb eq -1);
                $me = $i if ($desc[$i+1] ne "M");
                $order{"m"} = $i if(not $order{"m"});
               
                $assign_str{$i} = $col_mature;
                
            
                $lmature++;
                
            }elsif ($desc[$i] eq "l"){
                $lb = $i if($lb eq -1);
                $order{"l"} = $i if(not $order{"l"});
                $assign_str{$i} = $col_loop;
                $lloop++;
                if($star_exp_hit and $desc_tmp[$i] eq "S"){
                   # print STDERR "########## $i\t$desc_tmp[$i]\n";
                    $star_exp_hit_pos{$i}=$i;                  ##star hit expression
                }

            } elsif ($desc[$i] eq "S"){


                if($sb eq -1){
                    $sb = $i;
                    $sb_obs = $i;
                    my $run  = 1;
                    while(1){
                   #  for(my $run = 1;1;$run++){
                         if($desc_tmp[$i-$run] eq "S"){
                             $star_exp_hit_pos{$i-$run}=$i-$run;              ##star hit expression
                             $sb = $i-$run;
                             $run++;
                         
                         }else{
                             last;
                         }
                     }
                }
 
                
                $order{"s"} = $i if(not $order{"s"});
                if($hash{$sid}{'obs'}){
                    $assign_str{$i} = $col_star_obs;
                }else{
                    $assign_str{$i} = $col_star_exp;
                }
                $lstar++;
            } else {
                print STDERR "something went wrong while parsing alignment in output.mrd file $!\n";
            }   
        }
        

        if(not -d "$cwd/pdfs_$time"){
            mkdir "$cwd/pdfs_$time";
        }
        
        
        open FOLD,">$cwd/pdfs_$time/$filename.tmp" or die "Error: cannot create tmp file$!\n";
        
        $struct = $hash{$sid}{"pri_struct"};
        $lstruct_multi = ((length($struct)+2)*$multiplier);

        ## only for quantitation module
        if($mirbase){
            print FOLD "5${pri_seq}3\n";
        }else{ ## this is for the miRDeep module, where only the assumed precursor without flanks is folded

            if ($mb < $lb ){

                my $err = substr($pri_seq,$mb,$fl2b-$mb);
                print FOLD "5${err}3\n";
                $offset = $mb;
            } else {

                my $err = substr($pri_seq,$sb,$fl2b-$sb);

                print FOLD "5${err}3\n";

                $offset = $sb;

            }
        }
        close FOLD;
        
        for ($i=0; $i < length($struct);$i++)
        {
            $struct{$i} = 0;
        }


        ## new approach
        for my $tag(keys %{$hash{$sid}{"reads"}}){ ## for each tag
            for my $read(sort keys %{$hash{$sid}{"reads"}{$tag}}){
#                print "$tag x \t $read r \t $sid s \n";
                if ($hash{$sid}{"reads"}{$tag}{$read}{"seq"} =~ /^(\.*)(\w+)\.*$/){
                    my $v1=$1;
                    my $v2=$2;
                    
                    $hash2{$read}=length($v1); ## begin of read in precursor
                    if($hash{$sid}{"reads"}{$tag}{$read}{"rid"} =~ /_x(\d+)/){
                        my $dc = $1;
                        $totalreads+= $dc;
                    
                        $hash2c{$tag}{$v2}+=$1;     ## number of reads with same sequence
                        $hash2key{$read}=$v2;
                        $hash2order{$read} = $read;
                        $hash2mm{$read}=$hash{$sid}{"reads"}{$tag}{$read}{"mm"}; ## number of mismatches with precursor sequence
                        $hash2seq{$read} = $hash{$sid}{"reads"}{$tag}{$read}{"seq"};
                        
                        
                        $hash2sample{$read} = $tag;

                        for ($i=length($v1); $i < (length($v1)+length($v2)); $i++)
                        {            
                            #$struct{$i+1}+= $dc; ## saves how often a nt in precursor is covered by a read
                            $struct{$i}+= $dc; ## saves how often a nt in precursor is covered by a read
                        }
                    }
                    
                }
            }
        }                           ## end of reading aligned sequences
        $y = $yorig-$downy;
        chdir "./pdfs_$time";

        CreatePDF(%hash);
##############################################################################################
##############################################################################################
## insert secondary structure now
        DrawStructure($filename);
        if($totalreads ne '0'){
            CreateHistogram(%hash);
        }
        

        CreateAlignment(%hash);
        $y -=20;
        
### here the Frequency histogram is drawn
        
        
        
        ClosePDF($filename);
        unlink("$cwd/pdfs_$time/${filename}_ss.ps");
        unlink("$cwd/pdfs_$time/$filename.tmp");
        unlink("$cwd/pdfs_$time/rna.ps");
        chdir "..";
        print STDERR "creating pdf for $filename finished\n";
    }
    unlink("$cwd/pdfs_$time/tmp");
}

## create a PDF
sub CreatePDF{
    my %hash = @_;
    $pdf=PDF::API2->new; 
    #$spacer = length($sid);
    $spacer = 10;
	$pdf->mediabox('A4');
    $page=$pdf->page;
    $gfx=$page->gfx;
    $text=$gfx;
    $trb=$pdf->corefont('Times-Roman', -encode=>'latin1');
 
    
    ## move everything except the structure downwards if $mirbase is set
    my $madd = 0;
    if($mirbase){
        $madd  = 60;
    }
    $gfx->textlabel($xposshift+20,$y+300+$downy,$trb,8,"Provisional ID",-color=>'black');
    $gfx->textlabel($xposshift+100,$y+300+$downy,$trb,8,": $sid",-color=>'black');

    ##only print for discovered miRNAs
    if(not $mirbase){
    $spaces = " " x ($spacer - length($hash{$sid}{score}));
    $gfx->textlabel($xposshift+20,$y+290+$downy,$trb,8,"Score total",-color=>'black');            
    $gfx->textlabel($xposshift+100,$y+290+$downy,$trb,8,": $spaces$hash{$sid}{\"score\"}",-color=>'black');

    $spaces = " " x ($spacer - length($hash{$sid}{score_star}));
    $gfx->textlabel($xposshift+20,$y+280+$downy,$trb,8,"Score for star read(s)",-color=>'black'); 
    $gfx->textlabel($xposshift+100,$y+280+$downy,$trb,8,": $spaces$hash{$sid}{\"score_star\"}",-color=>'black');

    $spaces = " " x ($spacer - length($hash{$sid}{score_read}));
    $gfx->textlabel($xposshift+20,$y+270+$downy,$trb,8,"Score for read counts",-color=>'black');  
    $gfx->textlabel($xposshift+100,$y+270+$downy,$trb,8,": $spaces$hash{$sid}{\"score_read\"}",-color=>'black');

    $spaces = " " x ($spacer - length($hash{$sid}{score_mfe}));
    $gfx->textlabel($xposshift+20,$y+260+$downy,$trb,8,"Score for mfe",-color=>'black');          
    $gfx->textlabel($xposshift+100,$y+260+$downy,$trb,8,": $spaces$hash{$sid}{\"score_mfe\"}",-color=>'black');

    $spaces = " " x ($spacer - length($hash{$sid}{rand}));
    $gfx->textlabel($xposshift+20,$y+250+$downy,$trb,8,"Score for randfold",-color=>'black');     
    $gfx->textlabel($xposshift+100,$y+250+$downy,$trb,8,": $spaces$hash{$sid}{\"rand\"}",-color=>'black');

    $spaces = " " x ($spacer - length($hash{$sid}{score_cons}));
    $gfx->textlabel($xposshift+20,$y+240+$downy,$trb,8,"Score for cons. seed",-color=>'black');   
    $gfx->textlabel($xposshift+100,$y+240+$downy,$trb,8,": $spaces$hash{$sid}{\"score_cons\"}",-color=>'black');
}
    

    $spaces = " " x ($spacer - length($hash{$sid}{"freq_total"}));
    $gfx->textlabel($xposshift+20,$y+230+$madd+$downy,$trb,8,"Total read count",-color=>'black');       
    $gfx->textlabel($xposshift+100,$y+230+$madd+$downy,$trb,8,": $spaces$hash{$sid}{\"freq_total\"}",-color=>'black');

    $spaces = " " x ($spacer - length($hash{$sid}{"freq_mature"}));
    $gfx->textlabel($xposshift+20,$y+220+$madd+$downy,$trb,8,"Mature read count",-color=>'black');      
    $gfx->textlabel($xposshift+100,$y+220+$madd+$downy,$trb,8,": $spaces$hash{$sid}{\"freq_mature\"}",-color=>'black');

    if(not $mirbase){
        $spaces = " " x ($spacer - length($hash{$sid}{"freq_loop"}));
        $gfx->textlabel($xposshift+20,$y+210+$madd+$downy,$trb,8,"Loop read count",-color=>'black');        
        $gfx->textlabel($xposshift+100,$y+210+$madd+$downy,$trb,8,": $spaces$hash{$sid}{\"freq_loop\"}",-color=>'black');
    
    $spaces = " " x ($spacer - length($hash{$sid}{"freq_star"}));
    $gfx->textlabel($xposshift+20,$y+200+$madd+$downy,$trb,8,"Star read count",-color=>'black');        
    $gfx->textlabel($xposshift+100,$y+200+$madd+$downy,$trb,8,": $spaces$hash{$sid}{\"freq_star\"}",-color=>'black');
    }else{
        $spaces = " " x ($spacer - length($hash{$sid}{"freq_star"}));
        $gfx->textlabel($xposshift+20,$y+210+$madd+$downy,$trb,8,"Star read count",-color=>'black');        
        $gfx->textlabel($xposshift+100,$y+210+$madd+$downy,$trb,8,": $spaces$hash{$sid}{\"freq_star\"}",-color=>'black');
    }
    $trb=$pdf->corefont('Courier', -encode=>'latin1');
}



sub ClosePDF{
	my $file = shift;
	$file = "output" if($file eq"");
	$pdf->saveas("$cwd/pdfs_$time/$file.pdf");
}



## draw a line in PDF 
sub Line{
	my ($x1,$y1,$x2,$y2,$col,$width) = @_;
	$dline->linewidth($width);
	$dline->strokecolor($col);
	$dline->move($x1,$y1);
	$dline->line($x2,$y2);
	$dline->stroke;
}

## draw a letter in PDF
sub Base{
	my ($x1,$y1,$base,$col,$size) = @_;
	$trb=$pdf->corefont('Courierbold', -encode=>'latin1');
	$gfx->textlabel($x1,$y1,$trb,$size,$base,-color=>$col);
}	


## draw alignment between precursor sequence and read sequences
sub CreateAlignment{
    my %hash = @_;

	## draw left flank of sequence 
    $gfx->textlabel($xposshift+(18+(($mb+1)*$multiplier)),$y+40,$trb,6,$mb-$lflank1+1,-color=>'black');

    $gfx->textlabel($xposshift+(20+(($mb+1)*$multiplier)),$y+10,$trb,12,'Mature',-color=>$col_mature);
    
		
    $gfx->textlabel($xposshift+(18+(($lb+1)*$multiplier)),$y+40,$trb,6,$lb-$lflank1+1,-color=>'black');


#	$gfx->textlabel((20+(($lb)*$multiplier)),$y+10,$trb,15,'L',-color=>$col_loop);
    
    $gfx->textlabel($xposshift+(18+(($sb+1)*$multiplier)),$y+40,$trb,6,$sb-$lflank1+1,-color=>'black');


    #if(not $mirbase){
	if($hash{$sid}{'obs'} =~ /S/){
		$gfx->textlabel($xposshift+(20+(($sb_obs+1)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_obs);
	}elsif($hash{$sid}{'exp'} =~ /S/){
		$gfx->textlabel($xposshift+(20+(($sb+1)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_exp);
	}else{}


    #}
    $gfx->textlabel($xposshift+(18+(($fl2b+1)*$multiplier)),$y+40,$trb,6,$fl2b-$lflank1+1,-color=>'black');




	
    for(my $i=0; $i < scalar @rna; $i++){
	    $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});   
	}

    $gfx->textlabel($xposshift+25+ $lstruct_multi,$y,$trb,6,'-3\'' ,-color=>'black');
    $gfx->textlabel($xposshift+10,$y,$trb,6,'5\'-' ,-color=>'black');

    if($hash{$sid}{'obs'}){
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'obs' ,-color=>'black');
    }else{
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'exp' ,-color=>'black');
    }
    
    if($hash{$sid}{'obs'}){
        $y -= 10;
        for(my $i=0; $i < scalar @rna; $i++){
            $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str_exp{$i}); 
        }
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'exp' ,-color=>'black');
        
    }

	if($hash{$sid}{'known'}){
        $y -= 10;
		my $mid = $hash{$sid}{"known"};
		my $mse=  $knownones{$mid};
		

		## now check if mature of miRDeep is same as mirbase or not
		my $mbaseb = index($hash{$sid}{'pri_seq'},$mse);
		if($mbaseb>=0){
			##draw flank first
			$gfx->textlabel($xposshift+20+$multiplier,$y,$trb,6,substr($hash{$sid}{'pri_seq'},0,$mbaseb) ,-color=>'black');
			## draw mature part
			$gfx->textlabel($xposshift+20+(($mbaseb+1)*$multiplier),$y,$trb,6,$mse ,-color=>'green');
			##draw right flank last
			$gfx->textlabel($xposshift+20+((length($mse)+$mbaseb+1)*$multiplier),$y,$trb,6,substr($hash{$sid}{'pri_seq'},$mbaseb+length($mse)) ,-color=>'black');
		}else{
			print STDERR "mature sequence not found \n";
		}
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'known' ,-color=>'green');
    }    
    $y -= 10;
	my @structx = split(//,$struct);
	my $sadd = 0;
    
    $gfx->textlabel($position_hash{0},$y,$trb,6,$struct,-color=>'black');
    
    $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,'reads' ,-color=>'black');
    $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,'mm' ,-color=>'black');
    $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,'sample' ,-color=>'black');
    $y -= 10;    
    if($options{'o'}){
        for my $tag(keys %{$hash{$sid}{"reads"}}){
            for my $k(sort { $hash2order{$a} <=> $hash2order{$b} } keys %hash2order){
                next if($hash2sample{$k} ne $tag);
                
                $gfx->textlabel($position_hash{0},$y,$trb,6,$hash2seq{$k},-color=>'black');
                
                ## matches and read numbers
                $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,$hash2c{$tag}{$hash2key{$k}} ,-color=>'black');
                $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,$hash2mm{$k} ,-color=>'black');
                $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,$hash2sample{$k} ,-color=>'black');
                
                $y -= 10;
                if($y < 100){
                    #last;#***
					$page=$pdf->page();
					$pdf->mediabox('A4');
					$text2=$page->text();
					$gfx = $text2;
					$y=800;
					for(my $i=0; $i < scalar @rna; $i++){
						$gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});
					}
					$gfx->textlabel($xposshift+(20+(($mb+1)*$multiplier)),$y+10,$trb,12,'Mature',-color=>$col_mature);
					#$gfx->textlabel((20+(($lb)*$multiplier)),$y+10,$trb,12,'L',-color=>$col_loop);
					
					if($hash{$sid}{'obs'} =~ /S/){
						$gfx->textlabel($xposshift+(20+(($sb_obs+1)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_obs);
					}elsif($hash{$sid}{'exp'} =~ /S/){
						$gfx->textlabel($xposshift+(20+(($sb+1)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_exp);
					}else{}
					#  if(not $mirbase){
#                 if($hash{$sid}{'obs'}){
#                     $gfx->textlabel($xposshift+(20+(($sb)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_obs);
#                 }else{
#                     $gfx->textlabel($xposshift+(20+(($sb)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_exp);
#                 }
					$y-=10;
					#}
				}
                
            }
            $y-=10;
        }
    }else{
        for my $k(sort { $hash2order{$a} <=> $hash2order{$b} } keys %hash2order){
            my $tag = $hash2sample{$k};
            $gfx->textlabel($position_hash{0},$y,$trb,6,$hash2seq{$k},-color=>'black');
            
            ## matches and read numbers
            $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,$hash2c{$tag}{$hash2key{$k}} ,-color=>'black');
            $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,$hash2mm{$k} ,-color=>'black');
            $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,$hash2sample{$k} ,-color=>'black');
            
            $y -= 10;
            if($y < 100){
                #last;#***
                $page=$pdf->page();
                $pdf->mediabox('A4');
                $text2=$page->text();
                $gfx = $text2;
                $y=800;
                for(my $i=0; $i < scalar @rna; $i++){
                    $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});
				}
                $gfx->textlabel($xposshift+(20+(($mb+1)*$multiplier)),$y+10,$trb,12,'Mature',-color=>$col_mature);
                #$gfx->textlabel((20+(($lb)*$multiplier)),$y+10,$trb,12,'L',-color=>$col_loop);
                
                if($hash{$sid}{'obs'} =~ /S/){
                    $gfx->textlabel($xposshift+(20+(($sb+1)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_obs);
                }elsif($hash{$sid}{'exp'} =~ /S/){
                    $gfx->textlabel($xposshift+(20+(($sb+1)*$multiplier)),$y+10,$trb,12,'Star',-color=>$col_star_exp);
                }else{}
                $y-=10;
            }
        }
    }
}

## write description of graph
sub CreateDescription{
	my $trbb=$pdf->corefont('Courier-Bold', -encode=>'latin1');

	$gfx->textlabel(20,$y-10,$trbb,14,'FIGURE',-color=>'black');
	$gfx->textlabel(20,$y-20,$trb,8,"The plot shows the relative frequency of nucleotides of the precursor sequence in mapped read",-color=>'black');
	$gfx->textlabel(20,$y-30,$trb,8,"sequences. Below is shown the precursor sequence with mapped reads. Mismatches between",-color=>'black');
	
	$gfx->textlabel(20,$y-40,$trb,8,"read and precursor sequence are indicated by lowercase letter in the read sequence while matches",-color=>'black');
	$gfx->textlabel(20,$y-50,$trb,8,"are indicated by uppercase letters. The big letters M, L and S (red, blue, green) indicate",-color=>'black');
	
	$gfx->textlabel(20,$y-60,$trb,8,"the starting positions of the mature miRNA sequence, the loop sequence and the star sequence.",-color=>'black');
	$gfx->textlabel(20,$y-70,$trb,8,"reads is the number of reads that have the same sequence on the left. mm is the number of",-color=>'black');
        $gfx->textlabel(20,$y-80,$trb,8,"mismatches between the alignment of the read sequence and the precursor sequence.",-color=>'black');
}


sub CreateHistogram{
    my %hash = @_;
    $dline->linewidth(2);
    $y = $yorig-$downy;
	
    $dline = $page->gfx;
    
    ##draw axes
    $dline->strokecolor('black');
    
	
    $dline->move($xposshift+20,$y+160);
    $dline->line($xposshift+20,$y+50-1);
    $dline->stroke;
    $dline->move($xposshift+20,$y+50);
    $dline->strokecolor('grey');
    $dline->line($xposshift+20+$lstruct_multi,$y+50);
    $dline->stroke;
    $dline->strokecolor('black');
    $dline->move($xposshift+17+$lstruct_multi,$y+53);
    $dline->line($xposshift+20+$lstruct_multi,$y+50);
    $dline->line($xposshift+17+$lstruct_multi,$y+47);

    $dline->move($xposshift+17,$y+157);
    $dline->line($xposshift+20,$y+160);
    $dline->line($xposshift+23,$y+157);

    $dline->move($xposshift+17,$y+150); 
    $dline->line($xposshift+23,$y+150);

    $gfx->textlabel($xposshift+12,$y+165,$trb,6,"freq." ,-color=>'black'); 
    $gfx->textlabel($xposshift+$lstruct_multi,$y+40,$trb,8,"length" ,-color=>'black'); 

    $gfx->textlabel($xposshift+10,$y+148,$trb,6,"1" ,-color=>'black'); 


    $dline->move($xposshift+17,$y+125); ##.75
    $dline->line($xposshift+23,$y+125);
    $gfx->textlabel($xposshift+2,$y+122,$trb,6,"0.75" ,-color=>'black'); 

    $dline->move($xposshift+17,$y+100); ## .5
    $dline->line($xposshift+23,$y+100);
    $gfx->textlabel($xposshift+6,$y+98,$trb,6,"0.5" ,-color=>'black'); 

    $dline->move($xposshift+17,$y+75); ## .25
    $dline->line($xposshift+23,$y+75); ## .25
    $gfx->textlabel($xposshift+2,$y+73,$trb,6,"0.25",-color=>'black'); 
    
    $gfx->textlabel($xposshift+12,$y+48,$trb,6,"0" ,-color=>'black'); 
    $dline->stroke;
	
    ## draw flank1
    $dline->strokecolor('black');
    $dline->move($xposshift+20,$y+50);
	
## example case for mmu-mir-497
    
    for($i = 0; $i <= $lflank1; $i++){ #0..12 means to 13th char in string ## print one char further to have the transition
        $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
#        print "$i\t$position_hash{$i}\t",(($struct{$i}/$totalreads)*100)+$y+50,"\n";
    }
    
	$i-=1; ## $i is 13 now where mature starts
 

    $lastx = $position_hash{$i+1};
    $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
	
	
    $dline->stroke;
    $dline->strokecolor('black');
    $dline->linewidth(2);
    ## 
    @sorted = sort { $order{$a} <=> $order{$b} } keys %order; 
    #print @sorted;
    
	
    foreach my $k(@sorted){
        if($k eq "m"){
            $dline->strokecolor($col_mature);
            $dline->move($lastx,$lasty);
 
            ## graphical output corrections
            if($mb > $sb){
                for($i = ($mb+1); $i < $mb+$lmature ; $i++){
                
                    $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
                #print "$struct{$i}\t$totalreads\tvalue=",(($struct{$i}/$totalreads)*100),"\ty=",((($struct{$i}/$totalreads)*100)+$y+50),"\n";
                }
            }else{
                for($i = $mb; $i < $mb+$lmature ; $i++){
                    
                    $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
                    #print "$struct{$i}\t$totalreads\tvalue=",(($struct{$i}/$totalreads)*100),"\ty=",((($struct{$i}/$totalreads)*100)+$y+50),"\n";
                }
            }

            $i-=1; ## $i is 34 now

            $dline->stroke;
            $lastx = $position_hash{$i+1};
            $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
            $dline->strokecolor('black');

        }elsif($k eq "s"){
            if($hash{$sid}{'obs'}){
                $dline->strokecolor($col_star_obs);
            }else{
                $dline->strokecolor($col_star_exp);
            }

            $dline->move($lastx,$lasty);

            
            ## graphical output corrections
            if($sb > $mb){
                for($i = ($sb+1); $i < $sb+$lstar ; $i++){
                    $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
                }
            }else{
                for($i = $sb; $i < $sb+$lstar ; $i++){
                    $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
                }
            }

            $i-=1;
            $dline->stroke;
            $lastx = $position_hash{$i+1};
            $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
            $dline->strokecolor('black');
        }elsif($k eq "l"){
            $dline->strokecolor('orange');
            $dline->move($lastx,$lasty);
            for($i = $lb+1; $i <= $lb+$lloop ; $i++){
                $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
            }

            $i-=1;
            $dline->stroke;

            $lastx = $position_hash{$i+1};
            $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
            
            $dline->strokecolor('black');
        }elsif($k eq "f"){
            #$dline->linewidth(1);
            $dline->strokecolor('black');
            $dline->move($lastx,$lasty);
            for($i = $fl2b+1; $i <= $fl2b+$lflank2 ; $i++){
                $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
            }

            $i-=1;
            $dline->stroke;
            $lastx = $position_hash{$i+1};
            $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
            $dline->strokecolor('black');
        }else{
        }
    }

#    print "\n\t$lflank1\nmt\t$mb\nmature\t$lmature\n";
#    for(sort {$a <=> $b} keys %position_hash){
#        print "$_\t$struct{$_}\t$position_hash{$_}\t",(($struct{$_}/$totalreads)*100),"\n";
#    }
}



sub Shifting{
	$minx = min(values %xc);
	$maxx = max(values %xc);
	$miny = min(values %yc);
	$maxy = max(values %yc);


	## now place sec-structure in upper right corner
	my $shiftx = abs($minx)+10;
	my $shifty = abs($miny)+10;

	#shift everything to printable area;
	if($minx < 0){
		for(my $i=0; $i < scalar @rna_d-1; $i++){
			$xc{$i}+=$shiftx;
		}
	}
	if($miny < 0){
		for(my $i=0; $i < scalar @rna_d-1; $i++){
			$yc{$i}+=$shifty;
		}
	}
	
	if($maxx > 600){
		for(my $i=0; $i < scalar @rna_d-1; $i++){
			$xc{$i}-=$minx+10;
		}
	}
	if($maxy > 800){
		for(my $i=0; $i < scalar @rna_d-1; $i++){
			$yc{$i}-=$miny+10;
		}
	}	
}


sub DrawStructure{
    my $filename = shift;
    #return 
    $dline = $page->gfx;
    ## run RNAplot to create secondary structure ps file
    my $cw =cwd;
    #print STDERR "$sid\tcwd\t$cw\n";
    system("RNAfold -d 0 < $cwd/pdfs_$time/${filename}.tmp > $cwd/pdfs_$time/tmp");

    #exit;
    my $in_pos=0;
    my $in_pairs=0;
    
    my $count=0;                ## counter


    my %bp;                     ## which nt pair

    my $line;                   ## input line of file
    open PS,"<$cwd/pdfs_$time/rna.ps" or die "Error: cannot open $cwd/pdfs_$time/rna.ps\t$!\n";

    my ($minx, $miny) = 10000;
    my ($maxx,$maxy)  = 0;

    my $centering_x= 0; ## base in precursor sequences that is choosen as center point for rotation
    my $twisted=0; ## if mature sequence comes after loop twisted is 1
    my $sums=0;                 ## dif between 2bp nts

    while (<PS>)
    {
        if (/\/sequence/)       ## if sequence matched in rna.ps
        {
            $line = <PS>;       ## read in rna sequence
            chomp $line;
            #$line =~ s/U/T/g;		
            @rna_d=split(//,$line); ## read in to @rna_d
            next;
        }

        if (/\/coor\s*\[/) ## if nt coordinates section in ps comes now
        {
            $in_pos = 1;
            next;
        }
        if ($in_pos and /\[(\S+)\s+(\S+)\]/) 
        {
            $xc{$count} = $1;   ## x cooridnate
            $yc{$count} = $2;   ## y coordinate

            $count++;
            next;
        }

        if (/\/pairs/)          ## read in base pairs
        {
            $count=0;
            $in_pos=0;
            $in_pairs = 1;
            next;
        }
        
        $twisted = 1 if($mb > $lb ); ## mature begin is after loop begin
        

        if ($in_pairs and /\[(\S+)\s+(\S+)\]/)
        {
            if ($twisted)
            {
                $bpo2r=$1;
                $bpo1r=$2;
            } else {
                $bpo2r=$2;
                $bpo1r=$1;
            }

            ## determine two subsequent bases in mature having subsequent paired bases
            if ($bpo1r >= $mb-$offset and $bpo1r < $me-$offset and $centering_x==0)
            {
                if ($twisted)
                {
                    $sums = -1;
                } else {
                    $sums=1;
                }
                if (($bpo1r-$bpo1) == $sums and ($bpo2-$bpo2r ) == $sums)
                {
                    if($twisted){
                        $centering_x=$bpo1r;
                    }else{
                        $centering_x=$bpo1r-2;
                    }
                }
            }
            
            $bpo1= $bpo1r;
            $bpo2 = $bpo2r;
            

            $bp{$bpo1r-1}=$bpo2r-1; ## saving nt pairs in hash %bp
            next;
        }
        if ($in_pairs and /\]\s*def/) ## end of nt pairs in ps file 
        {
            $in_pairs = 0;
            last;
        }
    }
    close PS;

    Shifting();

    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);
    
    ##determine if mirror or not
    my $mir=0;
    
    $mir = 1 if($twisted);

    my $yshift=0;
    ########### mirror sequence so that loop is on the right hand side

    my $cshift = 3; ## determines the nt in mature sequence which should be the rotating center
    if ($mir)
    {
        for (my $i=0; $i < scalar @rna_d; $i++)
        {
            $xc{$i} = -$xc{$i};
        }
    }
    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);


    #translate back to positive area
    if ($mir)
    {
        for (my $i=0; $i < scalar @rna_d; $i++)
        {
            $xc{$i} += abs($minx)+10;
        }
    }

    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);


    my $ax = $xc{$centering_x};
    my $ay = $yc{$centering_x};


    ## point relativ to center
    my $bx;
    my $by;
    #print "twisted $twisted\n";
    if ($twisted)
    {
        $bx = $xc{$centering_x-1};
        $by = $yc{$centering_x-1};
    } else {
        $bx = $xc{$centering_x+1};
        $by = $yc{$centering_x+1};
    }


    my $gk = $by-$ay;
    my $ak = $bx-$ax;


    my $r = sqrt(($ak**2)+($gk**2));       

    my $phi = asin($gk/$r);
    
    if ($bx < $ax and $by > $ay)
    {
		$phi = 3.141593-$phi;
    }
    if ($bx <= $ax and $by <= $ay)
    {
		$phi *= (-1);
		$phi += 3.141593;
    }

    my $alpha;
    my $do_rot = 1;
    if ($do_rot)
    {
        my $last = $xc{0};
        ### now rotate every point in a designated angle of phi
        for (my $i=0; $i < scalar @rna_d -1; $i++)
        {
            next if ($i == $centering_x);

            $bx = $xc{$i};
            $by = $yc{$i};

            $gk = $by-$ay;
            $ak = $bx-$ax;             
            $r = sqrt(($ak**2)+($gk**2));

            $alpha = asin($gk/$r);

            if ($bx < $ax and $by > $ay)
            {
                $alpha = 3.141593-$alpha;
            }
            if ($bx <= $ax and $by <= $ay)
            {
                $alpha *= (-1);
                $alpha += 3.141593;
            }

            $alpha -= $phi;
            
            $xc{$i} = $ax + $r*cos($alpha);
            $yc{$i} = $ay + $r*sin($alpha);

            my $dif =  ($xc{$i} - $last);
            $last=$xc{$i};
        }
    }
    
    
    my $reduce = 0;
    my $red_dist= abs($xc{$mb+$cshift}-$xc{$mb-1+$cshift});



    Shifting();
    
 #    ## scaling now structure

#     my ($scx,$scy) = (250,250); ## scaling center coordinates on page
#     my ($tx,$ty);               ## translated x and y coordinates
#     my $scfactor = 0.6;         ## scaling factor


#     for (my $i=0; $i < scalar @rna_d -1; $i++)
#     {
#         $tx = $xc{$i}-$scx;
#         $ty = $yc{$i}-$scy;
#         $xc{$i} = $tx*$scfactor + $scx;
#         $yc{$i} = $ty*$scfactor + $scx;
#     }


    ## check if to mirror horizontally again because RNAfold does not take care of sequence input direction
    if (not $twisted)
    {
        my @bpkeys = sort keys %bp;
        $maxy = max(values %yc);
        if ($yc{$bpkeys[0]} < $yc{$bp{$bpkeys[0]}})
        {
            for (my $i=0; $i < scalar @rna_d-1; $i++)
            {
                $yc{$i} = -$yc{$i}+$maxy+10;
            }
        }
    }



    ## draw structure
    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);

    $y = $yorig+300;
    my $x = 550;

    my $rx = 300;

    if ($maxx < $x)
    {
        $x = $x-$maxx; 
    } else {
        $x = -abs($x-$maxx); 
    }
    
    if ($maxy < $y)
    {
        $y = $y-$maxy; 
    } else {
        $y = -abs($y-$maxy); 
    }

 #   print "$minx+$x\t$maxx+$x\n";




 ## scaling now structure again

    my ($scx,$scy) = (250,250); ## scaling center coordinates on page
    my ($tx,$ty);               ## translated x and y coordinates
    my $scfactor = ($rx/($maxx-$minx));         ## scaling factor
 #   print $scfactor,"\n";
    
    if($scfactor < 1){

        for (my $i=0; $i < scalar @rna_d -1; $i++)
        {
            $tx = $xc{$i}-$scx;
            $ty = $yc{$i}-$scy;
            $xc{$i} = $tx*$scfactor + $scx;
            $yc{$i} = $ty*$scfactor + $scx;
        }
    }
    
    $minx = min(values %xc);
    $maxx = max(values %xc);
    $miny = min(values %yc);
    $maxy = max(values %yc);

    $y = $yorig+280;#300;
    my $x = 550;


    if ($maxx < $x)
    {
        $x = $x-$maxx; 
    } else {
        $x = -abs($x-$maxx); 
    }
    
    if ($maxy < $y)
    {
        $y = $y-$maxy; 
    } else {
        $y = -abs($y-$maxy); 
    }



    if(($minx+$x) < ($x-$rx)){ ## left most x-coord for a base 
        for(keys %xc){
            $xc{$_} += ($x-$rx)-($minx+$x);
        }
    }




    for (my $i=0; $i < scalar @rna_d-1; $i++)
    {
        if ($i > 0)
        {
            Line($xc{$i-1}+$x+2.5,$yc{$i-1}+2+$y,$xc{$i}+2.5+$x,$yc{$i}+2+$y,"grey",0.5);
        }
        
    }

    Base($xc{0}+$x,$yc{0}+$y,"$rna_d[0]'",'black',8);
    for (my $i=1; $i < (scalar @rna_d)-2; $i++)
    {
        
        if($star_exp_hit_pos{$i-1+$offset}){
            Base($xc{$i}+$x,$yc{$i}+$y,$rna_d[$i],$col_star_exp,8);
        }else{
            Base($xc{$i}+$x,$yc{$i}+$y,$rna_d[$i],$assign_str{$i-1+$offset},8);
        }
    }
    Base($xc{(scalar @rna_d)-2}+$x,$yc{(scalar @rna_d)-2}+$y,"$rna_d[(scalar @rna_d)-2]'",'black',8);

    ## drawing bp lines
    my $scfactorl= 0.4;         #scaling factor

    my $fx;                     ## from x coordinate
    my $tox;                    ## from y coordinate

    my $fy;                     ## from y cooridinate
    my $toy;                    ## to y coordinate
    
    my $dx;                     ## xlength
    my $dy;                     ## y length

    my $dx1;   ## difference between orig x length and scaled x length
    my $dy1;   ## difference between orig y length and scaled y length


    for (keys %bp)
    {
        $dx = abs($xc{$_} - $xc{$bp{$_}});
        $dy = abs($yc{$_} - $yc{$bp{$_}});

        $dx1 = ($dx-$scfactorl*$dx)/2;
        $dy1 = ($dy-$scfactorl*$dy)/2;


        if ($xc{$_} > $xc{$bp{$_}})
        {
            $fx = $xc{$_} - $dx1;
            $tox = $xc{$bp{$_}} + $dx1;
        } else {
            $fx = $xc{$_} + $dx1;
            $tox = $xc{$bp{$_}} - $dx1;

        }

        if ($yc{$_} > $yc{$bp{$_}})
        {
            $fy = $yc{$_} - $dy1;
            $toy = $yc{$bp{$_}} + $dy1;
        } else {
            $fy = $yc{$_} + $dy1;
            $toy = $yc{$bp{$_}} - $dy1;
        }

        Line($fx+2.5+$x,$fy+2+$y,$tox+2.5+$x,$toy+2+$y,"black",0.5);	
    }
}

sub CreateHTML{
## print html
print HTML <<EOF;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN"
    "http://www.w3.org/TR/html4/strict.dtd">
    <html>
    <head>
    <title>miRDeep2</title>
    <!-- CSS code -->
    <style type="text/css">
    body{
      font: .75em Arial,sans-serif; 
      background: #FFFFFF; 
      color: #333333
    }
div#
container
{
  width: 1000px; 
  margin:0 auto
  }
h1{ 
  color: #F60; 
  margin: 1em 0 0; 
    letter-spacing: -2px; 
}
p{
  margin: 0 0 1.7em; 
}

a.tooltip{
	text-decoration:none;
	color:black;
        font-weight:bold;
}

a.tooltip:hover
{    position: relative;               
     background: transparent;          
	
}  

a.tooltip span  
{    position: absolute;               
     visibility: hidden;               

     width: 20em;                      
     top: 1.5em;             
     background: #ffffdd;                    
}

a.tooltip:hover span  
{    visibility: visible;   }


a.tooltip2{
	text-decoration:none;
	color:black;
        font-weight:bold;
}

a.tooltip2:hover
{    position: relative;               
     background: transparent;          
	
}  

a.tooltip2 span  
{    position: absolute;               
     visibility: hidden;               

     width: 20em;                      
     top: 1.5em;left:-10em;             
     background: #ffffdd;                    
}

a.tooltip2:hover span  
{    visibility: visible;   }


a.tooltip3{
	text-decoration:none;
	color:black;
        font-weight:bold;
}

a.tooltip3:hover
{    position: relative;               
     background: transparent;          
	
}  

a.tooltip3 span  
{    position: absolute;               
     visibility: hidden;               

     width: 20em;                      
     top: 3em;             
     background: #ffffdd;                    
}

a.tooltip3:hover span  
{    visibility: visible;   }


a.tooltip4{
	text-decoration:none;
	color:black;
        font-weight:bold;
}

a.tooltip4:hover
{    position: relative;               
     background: transparent;          
	
}  

a.tooltip4 span  
{    position: absolute;               
     visibility: hidden;               

     width: 20em;                      
     top:3em;left:-10em;             
     background: #ffffdd;                    
}

a.tooltip4:hover span  
{    visibility: visible;   }



</style>
        
    </head>
    <body>
    <table border="0" width="100%">
    <colgroup>
    <col width="5*">
    <col width="5*">
    </colgroup>
    <tr height="200" valign="top">
    <td><font face="Times New Roman" size="8">
    <b><a href="http://www.mdc-berlin.de/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/index.html" target="_blank" style="color:Black;text-decoration:none" title="miRDeep2 homepage at MDC-Berlin" >miRDeep2</a></b>
    
    
    <br>
    <br>
    </font></td>
    <td> <a href="http://www.mdc-berlin.de/en/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/index.html" target="_blank" >
    <img src="http://www.mdc-berlin.de/de/research/research_teams/systems_biology_of_gene_regulatory_elements/projects/miRDeep/mirdeep.jpg" style="border-style: none" name="precursor miRNA" title="miRDeep2 homepage at MDC-Berlin" align=right alt="miRDeep home"/></a></td>
    </tr>
    </table>
EOF
}

sub CloseHTML{
print HTML <<EOF;
</table>
</body>
</html>
EOF
}

sub ReadInParameters(){
	if(-f "mirdeep_runs/run_${time}/run_${time}_parameters"){
		print HTML "<h2>Parameters used</h2>\n";
		print HTML " <table border=\"0\">\n";
		print HTML "<tr><td>miRDeep2 version</td><td>$options{'V'}</td></tr><br>\n";
		open PAR,"<mirdeep_runs/run_${time}/run_${time}_parameters" or print STDERR "File with parameters could not be openend\n";
		while(<PAR>){
			if(/args\s+(\S.+)/){
				print HTML "<tr><td>Program call</td><td>$1</td><tr>\n";
				next;
			}elsif(/^dir with tmp files\s+(\S+)/){
				print HTML "<tr><td>Working directory</td><td>$1<td></tr><br>\n";
				next;
			}elsif(/^file_reads\s+(\S+)/){
				print HTML "<tr><td>Reads</td><td>$1</td></tr>\n";
				next;
			}elsif(/^file_genome\s+(\S+)/){
				print HTML "<tr><td>Genome</td><td>$1</td></tr>\n";
				next;
			}elsif(/^file_reads_vs_genome\s+(\S+)/){
				print HTML "<tr><td>Mappings</td><td>$1</td></tr>\n";
				next;
			}elsif(/^file_mature_ref_this_species\s+(\S+)/){
				print HTML "<tr><td>Reference mature miRNAs</td><td>$1</td></tr>\n";
				next;
			}elsif(/^file_mature_ref_other_species\s+(\S+)/){
				print HTML "<tr><td>Other mature miRNAs</td><td>$1</td></tr>\n";
				next;
			}elsif(/option{t}\s*=\s*(\S+)/){
				print HTML "<tr><td>Species</td><td>$1</td></tr>\n";
				next;
			}elsif(/option{q}\s*=\s*(\S+)/){
				print HTML "<tr><td>Expression analysis raw file</td><td>$1</td></tr>\n";
				next;
			}elsif(/option{a}\s*=\s*(\S+)/){
				print HTML "<tr><td>minimum read stack height</td><td>$1</td></tr>\n";
				next;
			}elsif(/option{b}\s*=\s*(\S+)/){
				print HTML "<tr><td>minimum score for novel precursors shown in table</td><td>$1</td></tr>\n";
				next;
			}elsif(/option{c}\s*=\s*(\S+)/){
				print HTML "<tr><td>randfold analysis disabled</td><td>yes</td></tr>\n";
			}elsif(/option{v}\s*=\s*(\S+)/){
				print HTML "<tr><td>remove temporary files</td><td>yes</td></tr>\n";
				next;				
			}elsif(/started:\s*(\S+)/){
				print HTML "<tr><td>Start</td><td>$1</td></tr>\n";
				my $tmp=<PAR>;
				if($tmp =~ /ended:\s*(\S+)/){
					print HTML "<tr><td>End</td><td>$1</td></tr>\n";
				}
				$tmp=<PAR>;
				if($tmp =~ /total:\s*(\S+)/){
					print HTML "<tr><td>Total</td><td>$1</td></tr>\n";
				}
				next;			   
			}else{

			}
		}
		close PAR;
		print HTML "</table><br><br>";
	}else{
		print STDERR "File mirdeep_runs/run_${time}/run_${time}_parameters not found\n";
	}	
}



sub PrintHtmlTableHeader{
    my ($hl,$csv) = @_;

    my %h;


    ## divide string by linebreaks every x characters
    my $p1 ='<th><a href="" class="tooltip3">';
    my $p11='<th><a href="" class="tooltip4">';

    my $p2='<span>';
    my $q ='</span></a></th>';

    if($hl =~ /novel/i){
$h{1}{1} = 'provisional id';
$h{1}{2} = 'this is a provisional miRNA name assigned by miRDeep2. The first part of the id designates the chromosome or genome contig on which the miRNA gene is located. The second part is a running number that is added to avoid identical ids. The running number is incremented by one for each potential miRNA precursor that is excised from the genome. Clicking this field will display a pdf of the structure, read signature and score breakdown of the reported miRNA.';
$h{2}{1} = 'miRDeep2 score';
$h{2}{2} = 'the log-odds score assigned to the hairpin by miRDeep2';
$h{3}{1} = 'estimated probability that the miRNA candidate is a true positive';
$h{3}{2} = 'the estimated probability that a predicted novel miRNA with a score of this or higher is a true positive. To see exactly how this probability is estimated, mouse over the \'novel miRNAs, true positives\' in the table at the top of the webpage.';
$h{4}{1} = 'rfam alert';
$h{4}{2} = 'this field indicates if the predicted miRNA hairpin has sequence similarity to reference rRNAs or tRNAs. Warnings in this field should overrule the estimated probability that a reported miRNA is a true positive (previous field).';
$h{5}{1} = 'total read count';
$h{5}{2} = 'this is the sum of read counts for the predicted mature, loop and star miRNAs.';
$h{6}{1} = 'mature read count';
$h{6}{2} = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted mature miRNA, including 2 nts upstream and 5 nts downstream.';
$h{7}{1} = 'loop read count';
$h{7}{2} = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted miRNA loop, including 2 nts upstream and 5 nts downstream.';
$h{8}{1} = 'star read count';
$h{8}{2} = 'this is the number of reads that map to the predicted miRNA hairpin and are contained in the sequence covered by the predicted star miRNA, including 2 nts upstream and 5 nts downstream.';
$h{9}{1} = 'significant randfold p-value';
$h{9}{2} = 'this field indicates if the estimated randfold p-value of the excised potential miRNA hairpin is equal to or lower than 0.05 (see Bonnet et al., Bioinformatics, 2004).';
$h{10}{1} = 'miRBase miRNA';
$h{10}{2} = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA maps to the miRNA hairpin, then only the id of the reference miRBase miRNA that matches the predicted mature sequence is output.';
$h{11}{1} = 'example miRBase miRNA with the same seed';
$h{11}{2} = 'this field displays the ids of any reference mature miRNAs from related species that have a seed sequence identical to that of the reported mature miRNA. The seed is here defined as nucleotides 2-8 from the 5\' end of the mature miRNA. If more than one reference mature miRNA have identical seed, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs from related species is displayed.';
$h{12}{1} = 'UCSC browser';
$h{12}{2} = 'if a species name was input to miRDeep2, then clicking this field will initiate a UCSC blat search of the consensus precursor sequence against the reference genome.';
$h{13}{1} = 'NCBI blastn';
$h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the consensus precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
$h{14}{1} = 'consensus mature sequence';
$h{14}{2} = 'this is the consensus mature miRNA sequence as inferred from the deep sequencing reads.';
$h{15}{1} = 'consensus star sequence';
$h{15}{2} = 'this is the consensus star miRNA sequence as inferred from the deep sequencing reads.';
$h{16}{1} = 'consensus precursor sequence';
$h{16}{2} = 'this is the consensus precursor miRNA sequence as inferred from the deep sequencing reads. Note that this is the inferred Drosha hairpin product, and therefore does not include substantial flanking genomic sequence as does most miRBase precursors.';
$h{17}{1} = 'precursor coordinate';
$h{17}{2} = 'The given precursor coordinates refer do absolute position in the mapped reference sequence';
    

}elsif($hl =~ /miRBase miRNAs detected by miRDeep2/i){
$h{1}{1} = 'tag id';
$h{1}{2} = 'this is a tag id assigned by miRDeep2. The first part of the id designates the chromosome or genome contig on which the miRNA gene is located. The second part is a running number that is added to avoid identical ids. The running number is incremented by one for each potential miRNA precursor that is excised from the genome. Clicking this field will display a pdf of the structure, read signature and score breakdown of the miRNA.';
$h{2}{1} = 'miRDeep2 score';
$h{2}{2} = 'the log-odds score assigned to the hairpin by miRDeep2';
$h{3}{1} = 'estimated probability that the miRNA is a true positive';
$h{3}{2} = 'the estimated probability that a predicted miRNA with a score of this or higher is a true positive. To see exactly how this probability is estimated, mouse over the \'novel miRNAs, true positives\' in the table at the top of the webpage. For miRBase miRNAs, this reflects the support that the data at hand lends to the miRNA.';
$h{4}{1} = 'rfam alert';
$h{4}{2} = 'this field indicates if the miRNA hairpin has sequence similarity to reference rRNAs or tRNAs. Warnings in this field should overrule the estimated probability that a reported miRNA is a true positive (previous field).';
$h{5}{1} = 'total read count';
$h{5}{2} = 'this is the sum of read counts for the mature, loop and star miRNAs.';
$h{6}{1} = 'mature read count';
$h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus mature miRNA, including 2 nts upstream and 5 nts downstream.';
$h{7}{1} = 'loop read count';
$h{7}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus miRNA loop, including 2 nts upstream and 5 nts downstream.';
$h{8}{1} = 'star read count';
$h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the consensus star miRNA, including 2 nts upstream and 5 nts downstream.';
$h{9}{1} = 'significant randfold p-value';
$h{9}{2} = 'this field indicates if the estimated randfold p-value of the miRNA hairpin is equal to or lower than 0.05 (see Bonnet et al., Bioinformatics, 2004).';
$h{10}{1} = 'mature miRBase miRNA';
$h{10}{2} = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA maps to the miRNA hairpin, then only the id of the reference miRBase miRNA that matches the predicted mature sequence is output.';
$h{11}{1} = 'example miRBase miRNA with the same seed';
$h{11}{2} = 'this field displays the ids of any reference mature miRNAs from related species that have a seed sequence identical to that of the reported mature miRNA. The seed is here defined as nucleotides 2-8 from the 5\' end of the mature miRNA. If more than one reference mature miRNA have identical seed, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs from related species is displayed.';
$h{12}{1} = 'UCSC browser';
$h{12}{2} = 'if a species name was input to miRDeep2, then clicking this field will initiate a UCSC blat search of the consensus precursor sequence against the reference genome.';
$h{13}{1} = 'NCBI blastn';
$h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the consensus precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
$h{14}{1} = 'consensus mature sequence';
$h{14}{2} = 'this is the consensus mature miRNA sequence as inferred from the deep sequencing reads.';
$h{15}{1} = 'consensus star sequence';
$h{15}{2} = 'this is the consensus star miRNA sequence as inferred from the deep sequencing reads.';
$h{16}{1} = 'consensus precursor sequence';
$h{16}{2} = 'this is the consensus precursor miRNA sequence as inferred from the deep sequencing reads. Note that this is the inferred Drosha hairpin product, and therefore does not include substantial flanking genomic sequence as does most miRBase precursors.';
$h{17}{1} = 'precursor coordinate';
$h{17}{2} = 'The given precursor coordinates refer do absolute position in the mapped reference sequence';



$h{4}{3} = 'predicted mature seq. in accordance with miRBase mature seq.';
$h{4}{4} = 'If the predicted miRDeep2 sequence overlaps with the miRBase annotated mature sequence than this is indicated by \'TRUE\'. If the predicted miRDeep2 star sequence overlaps with the miRBase annotated mature sequence this is inidicated by \'STAR\'.';

if($options{'p'}){
   # $h{17}{1} = 'genomic precursor position';
}

}else{
$h{1}{1} = 'miRBase precursor id';
$h{1}{2} = 'Clicking this field will display a pdf of the structure and read signature of the miRNA.';
$h{2}{1} = '-';
$h{2}{1} .="&#160" x 5;
$h{2}{2} = '-';
$h{3}{1} = '-';
$h{3}{1} .="&#160" x 5;
$h{3}{2} = '-';
$h{4}{1} = '-';
$h{4}{1} .="&#160" x 5;
$h{4}{2} = '-';
$h{5}{1} = 'total read count';
$h{5}{2} = 'this is the sum of read counts for the mature and star miRNAs.';
$h{6}{1} = 'mature read count(s)';
$h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the mature miRNA, including 2 nts upstream and 5 nts downstream. If more than one mature sequence is given this will be a comman separated list';
$h{7}{1} = '-    ';
$h{7}{2} = '-    ';
$h{8}{1} = 'star read count';
$h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the star miRNA, including 2 nts upstream and 5 nts downstream. This field is empty unless a reference star miRNA was given as input to quantifier.pl. If more than one mature sequence is given this will be a comman separated list';
$h{9}{1} = 'remaining reads';
$h{9}{2} = 'this is the number of reads that did not map to any of the mature and star sequences';
$h{10}{1} = '-'; #'miRBase mature id';
$h{10}{2} = '-'; #'Clicking this field will link to miRBase.';
$h{11}{1} = '-';
$h{11}{2} = '-';
$h{12}{1} = 'UCSC browser';
$h{12}{2} = 'if a species name was input to miRDeep2, then clicking this field will initiate a UCSC blat search of the miRNA precursor sequence against the reference genome.';
$h{13}{1} = 'NCBI blastn';
$h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the miRNA precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
$h{14}{1} = 'miRBase mature sequence(s)';
$h{14}{2} = 'this is/are the mature miRNA sequence(s) input to quantifier.pl.';
$h{15}{1} = 'miRBase star sequence(s)';
$h{15}{2} = 'this is/are the star miRNA sequence(s) input to quantifier.pl. This field is empty unless a reference star miRNA was given as input to quantifier.pl.';
$h{16}{1} = 'miRBase precursor sequence';
$h{16}{2} = 'this is the precursor miRNA sequence input to quantifier.pl.';
$h{17}{1} = 'precursor coordinate';
$h{17}{2} = 'The given precursor coordinates refer do absolute position in the mapped reference sequence';

if($options{'p'}){
   # $h{17}{1} = 'genomic precursor position';
}

}


    if($csv){
        my $f=1;
        print CSV "$hl\n";
        for(sort {$a <=> $b} keys %h){
            if($f){
                $f =0;
                print CSV "$h{$_}{1}";
            }else{
				if($hl =~ /miRBase miRNAs detected/i and $_ == 5){
					print CSV "$h{'mirbase'}{1}";
				}
                print CSV "\t$h{$_}{1}";
            }
        }
        print CSV "\n";
        return;
    }


print HTML <<EOF;
  
<br>
    <br>
    <br><h2>$hl</h2><br>
    <font face="Times New Roman" size="2">
    <table border="1">
EOF
for(sort {$a <=> $b} keys %h){
    if($_ ne 16){
        if($_ eq 9 and $hl !~ /not/){

           print HTML "<th><a href=\"http://www.ncbi.nlm.nih.gov/entrez/utils/fref.fcgi?PrId=3051&itool=AbstractPlus-def&uid=15217813&nlmid=9808944&db=pubmed&url=http://bioinformatics.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=15217813\" target=\"blank\"class=\"tooltip3\">$h{$_}{1}$p2$h{$_}{2}$q\n\n";

       }else{ 
 	    print HTML "$p1$h{$_}{1}$p2$h{$_}{2}$q\n\n";
	   if($hl =~ /miRBase miRNAs detected/ and $_ eq 4){
	       print HTML "$p1$h{$_}{3}$p2$h{$_}{4}$q\n\n";
	   }




	  
       }
	
    }else{

        print HTML "$p11$h{$_}{1}$p2$h{$_}{2}$q\n\n";
    }
}
}


sub PrintKnownnotfound{
    my %not_seen; 
    my %signature;
    
	my %error=();

	
    my ( $mature, $path0, $extension0 ) = fileparse ( $options{'k'}, '\..*' );
        ##filename
 
    $mature .= "${extension0}_mapped.arf"; ## change .fa suffix to _mapped.arf
    
    my %exprs;
	open IN,"<expression_analyses/expression_analyses_${time}/miRNA_expressed.csv" or die "Error: File expression_analyses/expression_analyses_${time}/miRNA_expressed.csv not found\n";
	while(<IN>){
        chomp;
        next if(/precursor/);
		my @line = split(/\t/);
		
        ## here comes up the trouble when two mature map to same precursor
        $mature2hairpin{$line[0]}=$line[2]; 
		$hairpin2mature{$line[2]}=$line[0];
        $exprs{$line[0]} = $line[1];
	}
	close IN;



	## read in error messages of mirdeep2 run
	if($options{'E'}){
		my $infile=$options{'s'};
		$infile =~ s/survey.csv/error.output.mrd/;
		open IN,$infile or print STDERR "$infile not found\n";
		my $in=0;
		my $id;
		my $mess='';
		while(<IN>){
			if(/>(\S+)/){$in=1;$id=$1;$mess='';next;}
			if(/^exp/){$in=0;next;}
			if($in){
				$mess.=$_;
			}
			if(not $in){
				next if(/pri_seq/);
				next if(/pri_struct/);
				if(/^(\S+)/){

					if($knownones{$1} and not $seen{$1}){
						$error{$1}{$id}.=$mess;
					}
				}
			}
		}
		close IN;
	}


	## only performed when option z is NOT given, so when the standard make_html.pl call is done
    if(not $options{'z'}){
        ## this are all reads of known miRBase miRNAs seen in output.mrd
        
        ## put all not seen in output.mrd miRNAs in new hash %not_seen
        for(keys %knownones){ ## all miRNAs in mature_ref_this_species listed things
            if(not $seen{$_}){
				#print STDERR "$_\n"; ## precursor printed out
                $not_seen{$_} = 1;
            }
        }
        
#### read in all micrornas in data
        
        if($options{'x'}){
            open IN,"<$options{'x'}" or die "Error: Cannot open file $options{'x'}\n";

            my @mmu;
            while(<IN>){
                next if(/ID\s+read count/);
                @mmu = split(/\t/);
				
                if($knownones{$mmu[0]} or $knownones{$mature2hairpin{$mmu[0]}}){
                    $signature{$mmu[0]} = 1;
                }
            }
            close IN;
        }
    }# else{ ## close options{'z'}
    #     ## create HTML for quantifier module
    #     open HTML,">expression_${time}.html" or die "cannot create expression_${time}.html\n";
    #     CreateHTML(); ##
    #     PrintHtmlTableHeader();
    # }
    ## now read in the mature_ref_this_species mapped against precursors from the quantifier module 
	## store ids in hashes
	
   
    
    ## open miRBase.mrd from quantifier module ## precursor ids as hash
    open IN,"<$options{'q'}" or die "Error: cannot open $options{'q'}\n";
    while(<IN>){
        if(/^\>(\S+)/){
            $id = $1;
            $oid = $1;
            $id =~ tr/|/_/;
            $hash_q{$id}{"oid"}=$oid;
            $hash_q{$id}{"id"}=$id;
            $counter =0;
        }
        elsif(/^remaining read count\s*(\d*)/){
            $hash_q{$id}{"remaining_loop"}=$1;
        }

        elsif(/^total read count\s*(\d*)/){
            $hash_q{$id}{"freq_total"}=$1;
        }

        elsif(/^(\S+) read count\s*(\d*)/){ ## everything else is just read in as id with real name
            $hash_q{$id}{$1}=$2;
        }

#       elsif(/^loop read count\s*(\d*)/){
#           $hash_q{$id}{"freq_loop"}=$1;
#       }
#        elsif(/^star read count\s*(\d*)/){
#            $hash_q{$id}{"freq_star"}=$1;
#        }
        elsif(/^pri_seq\s+(\S+)/){
            $hash_q{$id}{"pri_seq"}=$1;
            my @d;
            if($hash_q{$id}{'obs'}){
                @d = split(//,$hash_q{$id}{'obs'});
            }else{
                @d = split(//,$hash_q{$id}{'exp'});
            }

            ## put precursor sequence to array 
            my @s = split(//,$hash_q{$id}{"pri_seq"});
            my $mseq="";
            my $sseq="";

            ## now set labels for mature and star seq in explanation string
            for(my $i=0; $i < length($1); $i++){
                if($d[$i] ne "f"){ $hash_q{$id}{"ucsc_seq"} .= $s[$i];}
                if($d[$i] eq "M" or $d[$i] eq "5" or $d[$i] eq "3"){     ## accoutn for everything else
                    $mseq.= $s[$i];
                }elsif($d[$i] eq "S"){
                    $sseq.= $s[$i];
                }else{}
            }
            
            ## if there is an observed star sequence use this
            my $sseq_obs="";
            if($hash_q{$id}{'obs'}){
                @d = split(//,$hash_q{$id}{'obs'});
                for(my $i=0; $i < length($1); $i++){
                    $sseq_obs.= $s[$i] if($d[$i] eq "S");
                }
            }


            $hash_q{$id}{"mat_seq"} = $mseq;
            $hash_q{$id}{"star_seq"} = $sseq;
            $hash_q{$id}{"star_seq_obs"} = $sseq_obs;
        }
        elsif(/^exp\s+(\S+)/){
            $hash_q{$id}{'exp'}=$1;
        }
        elsif(/^obs\s+(\S+)/){
            $hash_q{$id}{'obs'}=$1;
        }
        elsif(/^pri_struct\s+(\S+)/){
            $hash_q{$id}{"pri_struct"}=$1;
            $reads=1;
            $counter = 0;

            next;
        
		}elsif(/^(\S\S\S)(\S+)\s+(\S+)\s+(\S+)$/ and $reads){  ## do not read in known miRNAs
			$counter++;
			$hash_q{$id}{"reads"}{$1}{$counter}{"rid"} = "$1$2";
			$hash_q{$id}{"reads"}{$1}{$counter}{"seq"} = $3;
			$hash_q{$id}{"reads"}{$1}{$counter}{"mm"}  = $4;
#         if($id =~/-1224/){
#             print "$id\t$_\n";
#         }

#         elsif(/^(\S+)\s+(\S+)\s+(\S+)$/){  ## hash all reads
#             $counter++;
#             $hash_q{$id}{"reads"}{$counter}{"rid"} = $1;
#             $hash_q{$id}{"reads"}{$counter}{"seq"} = $2;
#             $hash_q{$id}{"reads"}{$counter}{"mm"}  = $3;
        }else{}
    }
#    print "$hash_q{'hsa-mir-548f-3'}{'freq_total'}\n";
#    print "$hash_q{'hsa-mir-548f-3'}{'freq_mature'}\n";
#    die;
    close IN;
    print STDERR "parsing miRBase.mrd file finished\n";


    ##print out not scored but in data miRNAs by miRDeep2 to html.
    my $id;
    my $sf;
    my $mf;

	## print to csv these miRNAs as well
	if($csv){
		open CSV,">>$cwd/result_${time}.csv" or die "Error: cannot create csv file\n";
		print CSV "\n#miRBase miRNAs not detected by miRDeep2\n";
	}
	
    #if make_html was called by quantifier.pl then output all precursors that are expressed
    if($options{'f'} =~ /miRBase.mrd/){
#         for my $id(sort {$hash_q{$b}{'freq_mature'} <=> $hash_q{$a}{'freq_mature'}} keys %hash_q){
#             $oid = $hairpin2mature2{$id};
# 			## set blat link
#             if($org eq ""){
#                 $blat = '<td> - </td>';
#             }else{
#                 $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash_q{$id}{'pri_seq'}\" target=\"_blank\">blat</a></td>";
#             }

#             ##here the belonging precursor is shown
#             $known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$oid</a></td>";
			
#             $sf = $hash_q{$id}{"freq_star"};
#             $mf = $hash_q{$id}{"freq_mature"};
			
#             ## print miRBase id instead of precursor
			
#             if($options{'d'}){
#                 print HTML "<tr><td nowrap=\"nowrap\">$id</td>\n";
#             }else{
# 				if($hash_q{$id}{"freq_total"} > 0){
# 					print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$id.pdf\">$id</a></td>\n";
# 				}else{
# 					print HTML "<tr><td nowrap=\"nowrap\">$id</td>\n";
# 				}
#             }

#              print HTML <<EOF;
# 			<td/></td>
# 			<td></td>
#             <td>$hash_q{$id}{'rfam'}</td>
#             <td>$hash_q{$id}{"freq_total"}</td>
            
# 			<!--			<td>$hash_q{$id}{"freq_mature"}</td> -->
# 			<td>$mf</td>
# 			<!--            <td>$exprs{$id}</td>-->
#             <td>$hash_q{$id}{"freq_loop"}</td>
# 			<!--            <td>$hash_q{$id}{"freq_star"}</td> -->
# 			<td>$sf</td>
#             <td></td>
# 			$known
# 			<td nowrap="nowrap">$hash_q{$id}{"cons_seed"}</td>
#             $blat
#             <td><a href=${blast}$hash_q{$id}{'pri_seq'}&JOB_TITLE=$hash_q{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
#             <td>$hash_q{$id}{'mat_seq'}</td>
#             <td>$hash_q{$id}{'star_seq'}</td>
#             <td>$hash_q{$id}{'pri_seq'}</td>
# 			</tr>     
# EOF
# print CSV "$id\t-\t-\t-\t$hash_q{$id}{'freq_total'}\t$hash_q{$id}{'freq_mature'}\t$hash_q{$id}{'freq_loop'}\t$hash_q{$id}{'freq_star'}\t$oid\t-\t-\t-\t$hash_q{$id}{'mat_seq'}\t$hash_q{$id}{'star_seq'}\t$hash_q{$id}{'pri_seq'}\n" if($csv);

# 		}
# 		close CSV;
    }else{ ## if processing only the output of a normal mirdeep run
		
		## testing, Parse quantifier html and print then when necessary
		my %exhtml=();
		my $start=0;
		my $mirc;
		
		my $rt=0;

		open IN,"expression_analyses/expression_analyses_$options{'y'}/expression_$options{'y'}.html" or die "No file found called expression_analyses/expression_analyses_$options{'y'}/expression_$options{'y'}.html\n";
		while(<IN>){
			if(/miRBase precursor id/){
				$start=1;
			}
			next if(!$start);
			## read in header
			if(/nowrap/){
				$start=2;
			}
			if(/^\s+<\/tr>\s*$/){
				$start=3;
			}

			if($start == 1){
				$exhtml{'header'}.=$_;
			}elsif($start == 2){
				if(/pdf/){
					if(/([a-zA-Z0-9-]+)<\/a><\/td>/){
						$mirc=$1;
						$exhtml{$mirc}.=$_;
					}
				}else{
					$exhtml{$mirc}.=$_;
				}
			}

		}
		close IN;

		print HTML "<br><br><h2>mature miRBase miRNAs not detected by miRDeep2</h2><br><font face=\"Times New Roman\" size=\"2\">\n <table border=\"1\"> $exhtml{'header'}";
		my @csvout;
		while($exhtml{'header'} =~ /tooltip\d+\">(.+)<span>/g){
			push(@csvout,$1);
		}
		print CSV join("\t",@csvout),"\n" if($csv);
		
		my %skip;
		
		for my $oid(sort {$exprs{$b} <=> $exprs{$a}} keys %exprs){
			## if miRNA is in hash not_seen and signature or quantifier module then trigger output creation
			
			next if($skip{$mature2hairpin{$oid}});
			$skip{$mature2hairpin{$oid}}=1;
#			print STDERR "$oid\t$not_seen{$mature2hairpin{$oid}}\t$signature{$oid}\n";
			if(( ($not_seen{$oid} or $not_seen{$mature2hairpin{$oid}}) and $signature{$oid}) or $options{'z'}){ 
##				die "rehe";
				$exhtml{$mature2hairpin{$oid}}=~ s/<td>\s*<\/td>/<td>-<\/td>/g;
				#die $exhtml{$mature2hairpin{$oid}};
				@csvout=split("\n",$exhtml{$mature2hairpin{$oid}});
				
				my $count;
				my $icount;
				my $sum=0;
				my $in;
				my ($notp,$line);
				
				foreach(@csvout){
					$notp=0;
					$line=$_;

					if(/>\S+<\/a><\/td>/){
						if(/pdf\s*<\/div>(\S+)<\/a><\/td>/){
							print CSV $1,"\t" if($csv);
						}
						
						if(/nobr/){
							$icount=0;
							$in=1;
							while(/<nobr>(\d+)/g){
								$notp=1; ## do not print names of matures in here
								$in++;
								next if($in % 2 ==1); 
								$sum+=$1;
								$icount++;
								last if($count== $icount);
							}
							if($line =~ /^(.+<\/a>)<\/td><td><nobr>/){
								$line = $1;
							}
						}
					}elsif($icount > 0 and /<\/table>/){
						print CSV "$sum\t" if($csv);
						print HTML "$sum</td>";
						$icount =0;
						$sum=0;
						
						

					}elsif(/norm/){
						$count=0;
						
						while(/norm/g){ $count++; }
						#$notp=1;
						$line=substr($line,0,7);
						

					}elsif(/>(\d+\.*\d*)<\/td>/){
						print CSV $1,"\t" if($csv);
					}elsif(/<td>\s*([acgtACGTuU]+)\s*<\/td>/){
						while(/<td>\s*([acgtACGTuU]+)\s*<\/td>/g){
							print CSV $1,"\t" if($csv);
						}
					}elsif(/"nowrap">(\d+)<\/td>/){
						print CSV $1,"\t" if($csv);
					}else{

					}

					if($notp){
					}else{
						if(/^\s*-\s*$/ and $csv){
							print CSV "-\t";
						}elsif(/<td>\s*-\s*<\/td>/ and $csv){
							print CSV "-\t";
						}
						print CSV "-\t" if(/blast.ncbi/ and $csv);
						print HTML $line;
					}
		   
				}
				print CSV "\n" if($csv);

				## print to HTML now
				
				# while($exhtml{$mature2hairpin{$oid}} =~ /<td>(\S+)<\/td>/g){
				# 	print "$1\n";


# 				exit;
# 				print CSV "$id\t-\t-\t-\t$hash_q{$id}{'freq_total'}\t$hash_q{$id}{'freq_mature'}\t$hash_q{$id}{'freq_loop'}\t$hash_q{$id}{'freq_star'}\t$oid\t-\t-\t-\t$hash_q{$id}{'mat_seq'}\t$hash_q{$id}{'star_seq'}\t$hash_q{$id}{'pri_seq'}\n" if($csv);





# 				next;### this was executed before

				


# 				if($oid =~ /^(\S+)-5p/){
# 					next if($exprs{"${1}-3p"} > $exprs{$oid});
# 				}elsif($oid =~ /^(\S+)-3p/){
# 					next if($exprs{"${1}-5p"} >= $exprs{$oid});
# 				}else{
# 				}
                
				
				
# 				$id = $mature2hairpin{$oid}; ## this is now the name of the precursor
# #            if($id =~ /-1224/){print "$id\n";}
# 				$hash_q{$id}{"pdf"} = 1;
# 				$sig++;
#             ## set blat link
# 				if($org eq ""){
# 					$blat = '<td> - </td>';
# 				}else{
# 					$blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash_q{$id}{'pri_seq'}\" target=\"_blank\">blat</a></td>";
# 				}
# 				my $s_star=$hash_q{$id}{'star_seq'};
# 				$s_star=$hash_q{$id}{'star_seq_obs'} if($hash_q{$id}{'star_seq_obs'});
# 				my $s_mat = $hash_q{$id}{'mat_seq'};
				
# 				##here the belonging precursor is shown
# 				$known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$oid</a></td>";
				
# 				$sf = $hash_q{$id}{"freq_star"};
# 				$mf = $hash_q{$id}{"freq_mature"};
				
# 				if($hash_q{$id}{"freq_mature"} ne $exprs{$oid} ){
# 					$mf = $exprs{$oid};
# 					if($sf){
# 						$sf = $hash_q{$id}{"freq_mature"};
# 					}
# 				}
				
# 				if($options{'z'}){
# 					my $m;
					
# 					if($oid =~ /3p/){
# 						$m = $hash_q{$id}{'mat_seq'};
# 						$s_mat = $s_star;
# 						$s_star = $m;
# 					}
# 					if($hash_q{$id}{"freq_mature"} < $hash_q{$id}{"freq_star"}){
# 						## swap read count entries
# 						#$m = $hash_q{$id}{"freq_mature"};
# 						#$hash_q{$id}{"freq_mature"} = $hash_q{$id}{"freq_star"};
# 						#$hash_q{$id}{"freq_star"} = $m;
# #                    $mf = $hash_q{$id}{"freq_star"};
# #                    $sf = $hash_q{$id}{"freq_mature"};
						
# 						## swap sequence entries
						
						
# 						$known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$oid</a></td>";
						
# 						my $a = $hash_q{$id}{'exp'};
# 						$a =~ s/M/Q/g;
# 						$a =~ s/S/M/g;
# 						$a =~ s/Q/S/g;                    
# 						$hash_q{$id}{'exp'} =$a;
						
# 					}else{
# 						if($oid =~ /3p/){
# 							$known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$oid</a></td>";
# 						}
# 					}
					
# 				}
				
# 				## print miRBase id instead of precursor
				
# 				if($options{'d'}){
# 					print HTML "<tr><td nowrap=\"nowrap\">$id</td>\n";
# 				}else{
# 					if($hash_q{$id}{"freq_total"} > 0){
# 						print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$id.pdf\">$id</a></td>\n";
# 					}else{
# 						print HTML "<tr><td nowrap=\"nowrap\">$id</td>\n";
# 					}
# 				}

#             print HTML "<td>";   

# 			for my $k(keys %{$error{$id}}){
# #				die "id: $k\n$error{$id}{$k}\n";
# 				print HTML "id: $k\n$error{$id}{$k}";
# 			}


#                 print HTML <<EOF;
# 				</td>
# 					<td></td>
# 					<td>$hash_q{$id}{'rfam'}</td>
# 					<td>$hash_q{$id}{"freq_total"}</td>
					
# 					<!--			<td>$hash_q{$id}{"freq_mature"}</td> -->
# 					<td>$mf</td>
# <!--            <td>$exprs{$oid}</td>-->
# <td>$hash_q{$id}{"freq_loop"}</td>
# <!--            <td>$hash_q{$id}{"freq_star"}</td> -->
# <td>$sf</td>
# <td></td>
# $known
# <td nowrap="nowrap">$hash_q{$id}{"cons_seed"}</td>
# $blat
# <td><a href=${blast}$hash_q{$id}{'pri_seq'}&JOB_TITLE=$hash_q{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
# <td>$s_mat</td>
# <td>$s_star</td>
# <td>$hash_q{$id}{'pri_seq'}</td>
# </tr>     
# EOF
# if(not $hash_q{$id}{'freq_total'}){$hash_q{$id}{'freq_total'}=0;}
# if(not $hash_q{$id}{'freq_mature'}){$hash_q{$id}{'freq_mature'}=0;}
# if(not $hash_q{$id}{'freq_loop'}){$hash_q{$id}{'freq_loop'}=0;}
# if(not $hash_q{$id}{'freq_star'}){$hash_q{$id}{'freq_star'}=0;}
# if(not $hash_q{$id}{'star_seq'}){$hash_q{$id}{'star_seq'}='-';}
# print CSV "$id\t-\t-\t-\t$hash_q{$id}{'freq_total'}\t$hash_q{$id}{'freq_mature'}\t$hash_q{$id}{'freq_loop'}\t$hash_q{$id}{'freq_star'}\t$oid\t-\t-\t-\t$hash_q{$id}{'mat_seq'}\t$hash_q{$id}{'star_seq'}\t$hash_q{$id}{'pri_seq'}\n" if($csv);

# ##print also star if mapped reads greater than to the mature


# 		 if($hash_q{$id}{'freq_star'} >= $hash_q{$id}{'freq_mature'} and 0){
					
#         ## set blat link
#         if($org eq ""){
#             $blat = '<td> - </td>';
#         }else{
#             $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash_q{$id}{'pri_seq'}\" target=\"_blank\">blat</a></td>";
#         }
#         my $s_star=$hash_q{$id}{'star_seq'};
#         $s_star=$hash_q{$id}{'star_seq_obs'} if($hash_q{$id}{'star_seq_obs'});
        

#         ##here the belonging precursor is shown
#         $known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$_\" target=\"_blank\">$_*</a></td>";
#         $e_count++;

#            if($options{'d'}){
#                 print HTML "<tr><td nowrap=\"nowrap\">$_*</td>\n";
#             }else{
#                 print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$_.pdf\">$_*</a></td>\n";
#             }



#                 print HTML <<EOF;

# 			<td></td>
#             <td></td>
#             <td>$hash_q{$id}{'rfam'}</td>
# 			<td>$hash_q{$id}{"freq_total"}</td>
# <!--			<td>$hash_q{$id}{"freq_mature"}</td> -->
# <td>$mf</td>
#             <td>$hash_q{$id}{"freq_loop"}</td>
# <!--            <td>$hash_q{$id}{"freq_star"}</td> -->
# <td>$sf</td>
#             <td></td>
# 			$known
# 			<td nowrap="nowrap">$hash_q{$id}{"cons_seed"}</td>
#             $blat
#             <td><a href=${blast}$hash_q{$id}{'pri_seq'}&JOB_TITLE=$hash_q{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
#             <td>$s_star</td>
#             <td>$hash_q{$id}{'mat_seq'}</td>
#             <td>$hash_q{$id}{'pri_seq'}</td>
# 			</tr>     
# EOF
# 		 }
 			}

		}
	}
	close CSV;
    return;
    ## read in all miRBase not in data
    open IN,"expression_analyses/expression_analyses_${time}/miRNA_not_expressed.csv" or die "expression_analyses/expression_analyses_${time}/miRNA_not_expressed.csv not found\nplease run quantifier module again and check for errors\n\n";

    my %not_exprs;
    
    while(<IN>){
        my $id;
        my $v;
        if(/(\S+)\s+(\d+)\s*/){
            $id = $1;
            $v = $2;
            if($id !~ /\*/){
                $not_exprs{$id} = $v;
            }
        }
    }
    close IN;

    if($options{'z'}){ ## make 1 if you want to see all the not expressed miRBase miRNAs    
        
        if(scalar(keys %not_exprs) ne 0){
            print HTML "<tr bgcolor=silver><td colspan =16 >miRBase miRNAs not in data</td></tr>\n";
        

    ##print out miRNA not in data
    my $id;
    for(sort keys %not_exprs){
        $id = $mature2hairpin{$_}; ## this is now the name of the precursor
           
        ## set blat link
        if($org eq ""){
            $blat = '<td> - </td>';
        }else{
            $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash_q{$id}{'pri_seq'}\" target=\"_blank\">blat</a></td>";
        }
        my $s_star=$hash_q{$id}{'star_seq'};
        $s_star=$hash_q{$id}{'star_seq_obs'} if($hash_q{$id}{'star_seq_obs'});
        
        ##here the belonging precursor is shown
        $known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$_\" target=\"_blank\">$_</a></td>";

        my $ms =$hash_q{$id}{'mat_seq'};
        my $sse =$s_star;

        if(/3p/){
             $ms = $s_star;
             $sse = $hash_q{$id}{'mat_seq'};
         }

            
            

  if($options{'d'}){
                print HTML "<tr><td nowrap=\"nowrap\">$_</td>\n";
            }else{
                print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$_.pdf\">$_</a></td>\n";
            }
                print HTML <<EOF;

			<td></td>
            <td></td>
            <td>$hash_q{$id}{'rfam'}</td>
			<td>0</td>
			<td></td>
            <td></td>
            <td></td>
            <td></td>
			$known
			<td nowrap="nowrap">$hash_q{$id}{"cons_seed"}</td>
            $blat
            <td><a href=${blast}$hash_q{$id}{'pri_seq'}&JOB_TITLE=$hash_q{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
            <td>$ms</td>
            <td>$sse</td>
            <td>$hash_q{$id}{'pri_seq'}</td>
			</tr>     
EOF
}
        }

}
}## end of function

sub check_Rfam{
    my %hash = @_;
    
    my $err;
    if($options{'q'}){
        open TMP,">expression_analyses/expression_analyses_${time}/identified_precursors.fa" or die "Error: could not create file expression_analyses/expression_analyses_${time}/identified_precursors.fa\n";
    }else{
        open TMP,">mirdeep_runs/run_${time}/identified_precursors.fa" or die "Error: could not create file mirdeep_runs/run_${time}/identified_precursors.fa\n";
    }
    my ($seqm,$seql,$seqs);
    for(keys %hash){
        chomp;
        $seqm = $hash{$_}{'mat_seq'};
        $seqm =~ tr/uUacg/TTACG/;
        
        $seql = $hash{$_}{'loop_seq'};
        $seql =~ tr/uUacg/TTACG/;

        $seqs = $hash{$_}{'star_seq'};
        $seqs =~ tr/uUacg/TTACG/;
        

		print TMP ">$_\n$hash{$_}{'ucsc_seq'}\n";
        print TMP ">M:$_\n$seqm\n";
        if(length($seql) > 15){
            print TMP ">L:$_\n$seql\n";
        }
        print TMP ">S:$_\n$seqs\n";
    }
    close TMP;

    ## get script directory
    my $scripts=`which miRDeep2.pl`;
    $scripts =~ s/miRDeep2.pl//;
    $scripts =~ s/\s+//g;


    ## bowtie index is placed in folder indexes in folder that holds the mirdeep2 scripts
	my $tmp;
	print STDERR "Build bowtie index of Rfam entries\n\n";
    if(not -d "${scripts}indexes"){
        $tmp=`mkdir "${scripts}indexes"`;
		if($tmp){
			print STDERR $tmp,"\n";
		}
    }

    if(not -f "${scripts}indexes/Rfam_index.1.ebwt"){
        $err = `bowtie-build $options{'r'} ${scripts}indexes/Rfam_index`;
		if(not -f "${scripts}indexes/Rfam_index.1.ebwt"){
			print STDERR "$err\n\nRfam index could not be created in ${scripts}indexes/\nPlease check if miRDeep2 is allowed to create the directory ${scripts}indexes/ and write to it
The Rfam analysis will be skipped as long as this is not possible\n\n";
			return 256;
			
		}
	}
	
    #print STDERR "mapping precursor sequences against index\n";
    print STDERR "Mapping mature,star and loop sequences against index\n";
    ## I think 0 MM would be too conservative
    if($options{'q'}){
        $err = `bowtie -f -v 1 -a --best --strata --norc ${scripts}indexes/Rfam_index expression_analyses/expression_analyses_${time}/identified_precursors.fa expression_analyses/expression_analyses_${time}/rfam_vs_precursor.bwt`;
        open IN,"<expression_analyses/expression_analyses_${time}/rfam_vs_precursor.bwt" or die "Error: file expression_analyses/expression_analyses_${time}/rfam_vs_precursor.bwt not found\n";
    }else{
        $err = `bowtie -f -v 1 -a --best --strata --norc ${scripts}indexes/Rfam_index mirdeep_runs/run_${time}/identified_precursors.fa mirdeep_runs/run_${time}/rfam_vs_precursor.bwt`;
        open IN,"<mirdeep_runs/run_${time}/rfam_vs_precursor.bwt" or die "Error: file mirdeep_runs/run_${time}/rfam_vs_precursor.bwt not found\n";
    }
    
    
    my @line;
    my @ids;
    while(<IN>){
        @line = split(/\t/);
        @ids = split(/:/,$line[0]);
	


        if($line[2] =~ /rRNA/i){
            ##     id       M/S     type    counter
            $hash{$ids[1]}{$ids[0]}{'rRNA'}{'c'}++;
        }
        elsif($line[2] =~ /tRNA/i){
            $hash{$ids[1]}{$ids[0]}{'tRNA'}{'c'}++;
        }else{
			$hash{$ids[1]}{$ids[0]}{'rRNA'}{'c'} =0;
			$hash{$ids[1]}{$ids[0]}{'tRNA'}{'c'} =0;
			
            $hash{$ids[1]}{'rfam'} = '';
        }


# 		if(/chr19_7790/){
# 			die "$ids[1]\tfoundX\t$hash{$ids[1]}{$ids[0]}{'rRNA'}{'c'}\t$hash{$ids[1]}{$ids[0]}{'tRNA'}{'c'}\n";
# 		}


    }
    close IN;
    
    my @str = qw(M L S);

    my $count_o;

	


    for my $k(keys %hash){
        $count_o =0;
		
        for(my $i=0; $i < 3; $i++){
            $count_o++ if($hash{$k}{$str[$i]}{'rRNA'}{'c'} > 0);
        }
        if($count_o > 1){
            $hash{$k}{'rfam'} = 'rRNA';
            next;
        }
        
        $count_o = 0;
        for(my $i = 0; $i < 3; $i++){
            $count_o++ if($hash{$k}{$str[$i]}{'rRNA'}{'c'} > 0 or $hash{$k}{$str[$i]}{'tRNA'}{'c'} > 0);
        }
        if($count_o > 1){
            $hash{$k}{'rfam'} = 'rRNA/tRNA';
            next;
        }
        
        $count_o =0;
        for(my $i=0; $i < 3; $i++){
            $count_o++ if($hash{$k}{$str[$i]}{'tRNA'}{'c'} > 0);
        }
        if($count_o > 1){
            $hash{$k}{'rfam'} = 'tRNA';
            next;
        }
    }
# 	for(keys %hash){
# 		print STDERR "$k\t$hash{$k}{'rfam'}\n";
# 	}

    return 0;
}    

sub getEvalue{
    my ($v,$dig) = @_;

  if($v eq 0){ 
        return(0);
    }
    my $sign="";

    if($v < 0){
        $sign = "-";
        $v = abs($v);
    } 

    my $count =0;

    while($v >= 10){
        $count++;
        $v /= 10;
    }

    while($v < 1 and $v < 0.1){
        $count--;
        $v *= 10;
        if($v >= 1){
            last;
        }
    }

    if($v !~ /\./){
        $v .= ".";
    }


    if(length($v) > ($dig+2)){
        $v = substr($v,0,($dig+2));
    }else{
        $v .= "0" x ( ($dig+2)-length($v));
    }
    
    if(not $count){
        return("$sign$v");
    }else{
        if($count > 0){
            return("$sign${v}e+$count\n");
        }else{
            return("$sign${v}e$count\n");
        }
    }
}

sub get_precursor_pos{
	open IN,$options{'p'} or die "no precursor.coords file found with argument -p\n";
	my $id;
	while(<IN>){
		if(/^>*((\S+)_\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s*$/){
			$id = $1;
			$id =~ tr/|/_/;
			$pres_coords{$id}{'chr'} = $2;
			$pres_coords{$id}{'strand'} = $3;
			$pres_coords{$id}{'s'} = $4;
			$pres_coords{$id}{'e'} = $5;
		}
	}
	close IN;
}

sub get_mature_pos{
    open IN,$options{'f'} or die "no input file given\n";

    my %hash;
    my $id;
    my $start = -1;
    my $seq;
    my ($pos,$rpos,$pos_last,$end,$read);

    while(<IN>){
        chomp;
        if(/>(\S+)/){ $id = $1; }
        
        if(/exp\s+(\S+)/) { 
            $seq = $1;
            $start = index($seq,'M'); 
            $end = rindex($seq,'M'); 
            
        }
        
        if(/^(\S+_x\d+)\s+(\S+)/){
            $read = $1;
            $seq = $2;
            $seq =~ tr/acgtuACGTU/AAAAAAAAAA/;
            $pos = index($seq,'A');
            $rpos = rindex($seq,'A');
            if($pos eq $start and $rpos eq $end and not defined $hash{$read}){
                $hash{$id} = $read;
                
                #print "$read\t$id\n";

            }
        }
    }
    close IN;

    open IN,"mirdeep_runs/run_$time/tmp/signature.arf" or die "no signature file given for parsing\n";

    my @line;

    while(<IN>){
        chomp;
        @line = split();
        if($hash{$line[5]} and $hash{$line[5]} eq $line[0]){
            $mature_pos_hash{$line[5]}{'s'} = $line[7];
            $mature_pos_hash{$line[5]}{'e'} = $line[8];
            $mature_pos_hash{$line[5]}{'strand'} = $line[10];
#            die $line[5],"\n";
        }
        
    }
    close IN;
#    exit;
}
        

sub Usage{
	print STDERR "\n\n\n\n[usage]\n\tperl make_html.pl -f miRDeep_outfile [options]\n\n";
	print STDERR "[options]\n-v 2\t only output hairpins with score above 2\n";
    print STDERR "-c\t also create overview in excel format.\n";
    print STDERR "-k file\t supply file with known miRNAs\n";
    print STDERR "-s file\t supply survey file if score cutoff is used to get information about how big is the confidence of resulting reads\n\n\n";
    print STDERR "-e \t report complete survey file\n";
    print STDERR "-g \t report survey for current score cutoff\n";
#    print STDERR "-w project_folder\t automatically used when running webinterface, otherwise don't use it\n";
    print STDERR "-r file\t Rfam file to check for already reported small RNA sequences\n\n";
    print STDERR "-q file\t miRBase.mrd file produced by quantifier module\n";
    print STDERR "-x file\t signature.arf file with mapped reads to precursors\n";
    print STDERR "-t spec\t specify the organism from which your sequencing data was obtained\n";
    print STDERR "-u \t print all available UCSC input organisms\n";
    print STDERR "-y \ttimestamp of this run\n";
    print STDERR "-o \tsort signature by sample in pdf file, default is by beginning position\n";
    print STDERR "-d \t do not generate pdfs\n\n\n";
    print STDERR "-a \tprint genomic coordinates of mature sequence (still testing)\n";
    print STDERR "-b \tsupply confidence file\n";
	print STDERR "-V \tmiRDeep2 version used\n\n";
#	print STDERR "-E \t not fully implemented yet: prints error messages for known precursors that have not been scored by miRDeep2\n";

	exit;
}


__DATA__
Human
Chimp
Orangutan
Rhesus
Marmoset
Mouse
Rat
Guinea Pig
Cat
Dog
Horse
Cow
Opossum
Platypus
Chicken
Zebra finch
Lizard
X. tropicalis
Zebrafish
Tetraodon
Fugu
Stickleback
Medaka
Lamprey
Lancelet
C. intestinalis
S. purpuratus
C. elegans
C. brenneri
C. briggsae
C. remanei
C. japonica
P. pacificus
D. melanogaster
D. simulans
D. sechellia
D. yakuba
D. erecta
D. ananassae
D. pseudoobscura
D. persimilis
D. virilis
D. mojavensis
D. grimshawi
A. gambiae
A. mellifera
S. cerevisiae
