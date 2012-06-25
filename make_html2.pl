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
my %seen;


my $created = 1;

my $in;
my $counter=0;

my %struct; ## counts how often a nt is covered reads

my $i;

my $offset = 0;

my $me=0;     ## mature end coordinate
my @desc;

my %mat_pre_arf =();


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
my %star_exp_hit_pos;
my $blat;
my $spacer;   ## length of longest entry
my $spaces;   ## string of spaces to fill up spacer


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


## pdf histogram colors
my $col_star_exp = 'lightskyblue';
my $col_star_obs = 'darkviolet';
my $col_mature = 'red';
my $col_loop = 'orange';



my %hm;
my %hs;
my %hp;


## options
my %options=();

getopts("ugv:f:ck:os:t:w:er:q:dx:y:ab:i:j:lm:M:PW:",\%options);

my %weighted=();
if(-s $options{'W'}){
	open IN,$options{'W'} or die "Could not open file $options{'W'}\n";
	while(<IN>){
		if(/(\S+)\s+(\d+)/){
			$weighted{$1}=$2;
		}
	}
	close IN;
}



## everything else given to it corresponds to the samples
my @files_mirnaex=split(",",$options{'M'});
foreach(@files_mirnaex){
   print STDERR "$_ file with miRNA expression values\n";
}
 


#die "here in makehtml2\n";

my $time = $options{'y'} or die "no timestamp given with parameter y\n";

if($options{'x'} and not $options{'q'}){
    die "\nError:\n\toption -x can only be used together with option -q\n\n";
}

## determine pdf path when running on a cluster
my $pdf_path;

if($options{'w'}){
    $pdf_path = "http://localhost:8001/links/miRDeep/$options{'w'}";
}


## obtain current working directory
my $cwd = cwd;
if(not $options{'w'}){
    $pdf_path = "file://$cwd";
}

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
    open IN,"<$options{'b'}" or die "file not found\n";
    while(<IN>){
        next if(/\#/);
        chomp;
        my $r = $_;
        $r =~ s/\|/_/g;
        $confident{$r} = 1;
    }
    close IN;
}


PrintQuantifier();

CloseHTML();
system("cp expression_analyses/expression_analyses_${time}/expression_${time}.html expression_${time}.html");


if(not $options{'d'}){
	$mirbase = 1;
	CreateStructurePDFQuantifier(%hash_q);
}
exit;




sub CreateStructurePDFQuantifier{
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
        $sb=$mat_pre_arf{$sid}{'sb'};     ## starting 
        $lmature = 0; ## string of mature 
        $mb= $mat_pre_arf{$sid}{'mb'};     ## mature begin
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

        my @desc2 =split(//,$hash{$sid}{'exp'});
        for (my $i=0;$i < scalar @desc2; $i++){
            if ($desc2[$i] eq "f"){   ## assign_str now starts at 0 not at one
                $assign_str_exp{$i} = "black";
                $assign_str{$i} = "black";
            } elsif ($desc2[$i] eq "l"){
                $assign_str_exp{$i} = $col_loop;
                $assign_str{$i} = $col_loop;

            } elsif ($desc2[$i] eq "S"){
                $assign_str_exp{$i} = $col_star_obs;
                $assign_str{$i} = $col_star_obs;
            } else{
                $assign_str_exp{$i} = $col_mature;
                $assign_str{$i} = $col_mature;
                }   
        }
        
        if(not -d "$cwd/pdfs_$time"){
            mkdir "$cwd/pdfs_$time";
        }
        
        
        open FOLD,">$cwd/pdfs_$time/$filename.tmp" or die "Error: cannot create tmp file$!\n";
        
        $struct = $hash{$sid}{"pri_struct"};
        $lstruct_multi = ((length($struct)+2)*$multiplier);

        
        print FOLD "5${pri_seq}3\n";
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
						if($options{'W'}){
							$dc/=$weighted{$hash{$sid}{"reads"}{$tag}{$read}{"rid"}};
							#die "here  ===  $hash{$sid}{'reads'}{$tag}{$read}{'rid'} $dc\n";
						}

                        $totalreads+= $dc;
                    
                        #$hash2c{$tag}{$v2}+=$1;     ## number of reads with same sequence
						$hash2c{$tag}{$v2}+=$dc;
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
        ## by ref 
        CreatePDFQuantifier(\%hash);

##############################################################################################
##############################################################################################
## insert secondary structure now

        DrawStructure($filename);
        if($totalreads ne '0'){
            CreateHistogramQuantifier();
            #ClosePDF($filename);
            #exit;
}
        ## by ref
        CreateAlignmentQuantifier(\%hash);
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

sub CreateHistogramQuantifier{
#    my ($hash) = @_;
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

    $lastx = $position_hash{0};
    $lasty = (($struct{0}/$totalreads)*100)+$y+50;

    $dline->strokecolor($assign_str{0});
    $dline->move($lastx,$lasty);
    for($i = 0; $i < length($struct); $i++){
        $dline->strokecolor($assign_str{$i});
        $dline->move($lastx,$lasty);
        $dline->line($position_hash{$i+1},(($struct{$i}/$totalreads)*100)+$y+50); ## .25
        $dline->stroke;
        $lastx = $position_hash{$i+1};
        $lasty = (($struct{$i}/$totalreads)*100)+$y+50;
    }

}


sub CreatePDFQuantifier{
    my ($hash) = @_;
    $pdf=PDF::API2->new; 
	
    $spacer = length($sid);
    $pdf->mediabox('A4');
    $page=$pdf->page;
    $gfx=$page->gfx;
    $text=$gfx;
    $trb=$pdf->corefont('Times-Roman', -encode=>'latin1');
 
    
    ## move everything except the structure downwards if $mirbase is set
    my $madd = 60;

    $gfx->textlabel($xposshift+20,$y+300+$downy,$trb,8,"miRBase precursor",-color=>'black');
    $gfx->textlabel($xposshift+110,$y+300+$downy,$trb,8,": $sid",-color=>'black');


    $spaces = " " x ($spacer - length($$hash{$sid}{"freq_total"}));
    $gfx->textlabel($xposshift+20,$y+230+$madd+$downy,$trb,8,"Total read count",-color=>'black');       
    $gfx->textlabel($xposshift+110,$y+230+$madd+$downy,$trb,8,": $$hash{$sid}{'freq_total'}",-color=>'black');

    ## here should be written how many annotated stuff is actually there and how many not
    my $jk =10;
	# old
    #for(sort {$$hash{$sid}{'mapped'}{$b} <=> $$hash{$sid}{'mapped'}{$a}} keys %{$$hash{$sid}{'mapped'}}){
	#new
	my @h2m=split(",",$hairpin2mature{$sid});
	foreach my $h(@h2m){
		next if($h =~ /^\s*$/);

        if($options{'t'}){
            next if($_ !~ $options{'m'});
        }
        $spaces = " " x ($spacer - length($_));
        $gfx->textlabel($xposshift+20,$y+230-$jk+$madd+$downy,$trb,8,"$h read count",-color=>'black');      
        $gfx->textlabel($xposshift+110,$y+230-$jk+$madd+$downy,$trb,8,": $$hash{$sid}{'mapped'}{$h}",-color=>'black');
        $jk+=10;
    }

    $spaces = " " x ($spacer - length("remaining reads"));
    $gfx->textlabel($xposshift+20,$y+230-$jk+$madd+$downy,$trb,8,"remaining reads",-color=>'black');      
    $gfx->textlabel($xposshift+110,$y+230-$jk+$madd+$downy,$trb,8,": $$hash{$sid}{'remaining_rc'}",-color=>'black');
    $jk+=10;
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

sub CreateAlignmentQuantifier{
    my ($hash) = @_;

#    $gfx->textlabel($xposshift+(18+(($mb+1)*$multiplier)),$y+40,$trb,6,$mb-$lflank1+1,-color=>'black');
    $y+=20;
    for my $k1(keys %{$mat_pre_arf{$sid}}){
#        print STDERR "$sid\t$k1\n";
        if($k1 =~ /\*/){
            $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'start'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
        }else{
            $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'start'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_mature);
        }
        $y-=10;
    }

    for(my $i=0; $i < scalar @rna; $i++){
	    $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});   
	}

    $gfx->textlabel($xposshift+25+ $lstruct_multi,$y,$trb,6,'-3\'' ,-color=>'black');
    $gfx->textlabel($xposshift+10,$y,$trb,6,'5\'-' ,-color=>'black');

    if($$hash{$sid}{'obs'}){
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'obs' ,-color=>'black');
    }else{
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'exp' ,-color=>'black');
    }
    
    if($$hash{$sid}{'obs'}){
        $y -= 10;
        for(my $i=0; $i < scalar @rna; $i++){
            $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str_exp{$i}); 
        }
        $gfx->textlabel($xposshift+50+ $lstruct_multi,$y,$trb,6,'exp' ,-color=>'black');
        
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
        for my $tag(keys %{$$hash{$sid}{"reads"}}){
            for my $k(sort { $hash2order{$a} <=> $hash2order{$b} } keys %hash2order){
                next if($hash2sample{$k} ne $tag);
                
                $gfx->textlabel($position_hash{0},$y,$trb,6,$hash2seq{$k},-color=>'black');
                
                ## matches and read numbers
                $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,int($hash2c{$tag}{$hash2key{$k}}) ,-color=>'black');
                $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,$hash2mm{$k} ,-color=>'black');
                $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,$hash2sample{$k} ,-color=>'black');
                
                $y -= 10;
                if($y < 100){
					$page=$pdf->page();
					$pdf->mediabox('A4');
					$text2=$page->text();
					$gfx = $text2;
					$y=800;
					for(my $i=0; $i < scalar @rna; $i++){
						$gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});
					}
					$y+=20;
					for my $k1(keys %{$mat_pre_arf{$sid}}){
						if($k1 =~ /\*/){
							$gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'start'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
						}else{
							$gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'start'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_mature);
						}
						$y-=10;
					}
					$y-=20;
				}
            }
            $y-=10;
        }
    }else{
        for my $k(sort { $hash2order{$a} <=> $hash2order{$b} } keys %hash2order){
            my $tag = $hash2sample{$k};
            $gfx->textlabel($position_hash{0},$y,$trb,6,$hash2seq{$k},-color=>'black');
            
            ## matches and read numbers
            $gfx->textlabel($xposshift+30+ $lstruct_multi,$y,$trb,6,int($hash2c{$tag}{$hash2key{$k}}) ,-color=>'black');
            $gfx->textlabel($xposshift+70+ $lstruct_multi,$y,$trb,6,$hash2mm{$k} ,-color=>'black');
            $gfx->textlabel($xposshift+110+ $lstruct_multi,$y,$trb,6,$hash2sample{$k} ,-color=>'black');
            
            $y -= 10;
            if($y < 100){
                $page=$pdf->page();
                $pdf->mediabox('A4');
                $text2=$page->text();
                $gfx = $text2;
                $y=800;
				$y+=20;
				for my $k1(keys %{$mat_pre_arf{$sid}}){
                    if($k1 =~ /\*/){
                        $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'start'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_star_obs);
                    }else{
                        $gfx->textlabel($xposshift+(20+(($mat_pre_arf{$sid}{$k1}{'start'})*$multiplier)),$y,$trb,6,"$k1",-color=>$col_mature);
                    }
                    $y-=10;
                }
                for(my $i=0; $i < scalar @rna; $i++){
                    $gfx->textlabel($position_hash{$i},$y,$trb,6,$rna[$i],-color=>$assign_str{$i});
                }
                $y-=10;
            }
        }
    }
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
$h{10}{2} = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA map to the miRNA hairpin, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs for the species is displayed.';
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
if($options{'a'}){
    $h{17}{1} = 'genomic position';
}
    

}elsif($hl =~ /miRBase miRNAs in dataset/i){
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
$h{10}{1} = 'miRBase miRNA';
$h{10}{2} = 'this field displays the ids of any reference mature miRNAs for the species that map perfectly (full length, no mismatches) to the reported miRNA hairpin. If this is the case, the reported miRNA hairpin is assigned as a known miRNA. If not, it is assigned as a novel miRNA. If more than one reference mature miRNA map to the miRNA hairpin, then only the id of the miRNA that occurs last in the input file of reference mature miRNAs for the species is displayed. Clicking this field will link to miRBase.';
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
if($options{'a'}){
    $h{17}{1} = 'genomic position';
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

if($options{'P'}){
	$h{5}{2} = 'this is the sum of read counts for the 5p and 3p sequences.';
	$h{6}{1} = '5p read counts';
	$h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 5p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.';
	$h{8}{1} = '3p read counts';
	$h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the 3p sequence, including 2 nts upstream and 5 nts downstream. In parenthesis are normalized read counts shown.';
	$h{9}{1} = 'remaining reads';
	$h{9}{2} = 'this is the number of reads that did not map to any of the 5p or 3p sequences';

}else{
$h{5}{2} = 'this is the sum of read counts for the mature and star miRNAs.';
$h{6}{1} = 'mature read count(s)';
$h{6}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the mature miRNA, including 2 nts upstream and 5 nts downstream. If more than one mature sequence is given this will be a comma separated list. In parenthesis are normalized read counts shown.';
$h{8}{1} = 'star read count';
$h{8}{2} = 'this is the number of reads that map to the miRNA hairpin and are contained in the sequence covered by the star miRNA, including 2 nts upstream and 5 nts downstream. This field is empty unless a reference star miRNA was given as input to quantifier.pl. If more than one mature sequence is given this will be a comman separated list';
$h{9}{1} = 'remaining reads';
$h{9}{2} = 'this is the number of reads that did not map to any of the mature and star sequences';
}

$h{7}{1} = '-    ';
$h{7}{2} = '-    ';



$h{10}{1} = '-'; #'miRBase mature id';
$h{10}{2} = '-'; #'Clicking this field will link to miRBase.';
$h{11}{1} = '-';
$h{11}{2} = '-';
$h{12}{1} = 'UCSC browser';
$h{12}{2} = 'if a species name was input to miRDeep2, then clicking this field will initiate a UCSC blat search of the miRNA precursor sequence against the reference genome.';
$h{13}{1} = 'NCBI blastn';
$h{13}{2} = 'clicking this field will initiate a NCBI blastn search of the miRNA precursor sequence against the nr/nt database (non-redundant collection of all NCBI nucleotide sequences).';
if($options{'P'}){
	$h{14}{1} = 'miRBase 5p sequence(s)';
	$h{14}{2} = 'this is/are the 5p miRNA sequence(s) input to quantifier.pl.';
	$h{15}{1} = 'miRBase 3p sequence(s)';
	$h{15}{2} = 'this is/are the 3p miRNA sequence(s) input to quantifier.pl.';
}else{
	$h{14}{1} = 'miRBase mature sequence(s)';
	$h{14}{2} = 'this is/are the mature miRNA sequence(s) input to quantifier.pl.';
	$h{15}{1} = 'miRBase star sequence(s)';
	$h{15}{2} = 'this is/are the star miRNA sequence(s) input to quantifier.pl. This field is empty unless a reference star miRNA was given as input to quantifier.pl.';
}

$h{16}{1} = 'miRBase precursor sequence';
$h{16}{2} = 'this is the precursor miRNA sequence input to quantifier.pl.';
if($options{'a'}){
    $h{17}{1} = 'genomic position';
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

	next if($_ == 2);
	next if($_ == 3);
	next if($_ == 4);
	next if($_ == 7);
	next if($_ == 10);
	next if($_ == 11);
	

    if($_ ne 16){
        if($_ eq 9 and $hl !~ /not/){

           print HTML "<th><a href=\"http://www.ncbi.nlm.nih.gov/entrez/utils/fref.fcgi?PrId=3051&itool=AbstractPlus-def&uid=15217813&nlmid=9808944&db=pubmed&url=http://bioinformatics.oxfordjournals.org/cgi/pmidlookup?view=long&pmid=15217813\" target=\"blank\"class=\"tooltip3\">$h{$_}{1}$p2$h{$_}{2}$q\n\n";

       }else{ 


        print HTML "$p1$h{$_}{1}$p2$h{$_}{2}$q\n\n";
}

    }else{

        print HTML "$p11$h{$_}{1}$p2$h{$_}{2}$q\n\n";
    }
}
}


sub PrintQuantifier{
    my %not_seen; 
    my %signature;
    my $reads;

    ## create HTML for quantifier module
    open HTML,">expression_analyses/expression_analyses_${time}/expression_${time}.html" or die "cannot create expression_analyses/expression_analyses_${time}/expression_${time}.html\n";
    CreateHTML(); ##
    PrintHtmlTableHeader();
    
    ## now read in the mature_ref_this_species mapped against precursors from the quantifier module 
	## store ids in hashes
	
    my ( $mature, $path0, $extension0 ) = fileparse ( $options{'k'}, '\..*' );
    ##filename
    
    $mature .= '_mapped.arf'; ## change .fa suffix to _mapped.arf
    
    my %exprs;
	open IN,"<expression_analyses/expression_analyses_${time}/miRNA_expressed.csv" or die "Error: File expression_analyses/expression_analyses_${time}/miRNA_expressed.csv not found\n";
    
    ## this is the replacement for exprs hash
    my %exprs2;
    
	while(<IN>){
        chomp;
        next if(/precursor/);
		my @line = split(/\t/);
		
        ## here comes up the trouble when two mature map to same precursor
        $mature2hairpin{$line[0]}=$line[2]; 
		$hairpin2mature{$line[2]}.="$line[0],";
        $exprs{$line[0]} = $line[1];
        $exprs2{$line[2]}{$line[0]} = $line[1];
	}
	close IN;


    my %exprs_sample;
    ## read in sample stuff to expression_value hash
    foreach(@files_mirnaex){
        open IN,"<$_" or die "File $_ not asdfsd found\n";
#        my $sample = $1 if(/_(\S\S\S).csv$/);
		my @samples;
        while(<IN>){
            chomp;
            my @line = split(/\t/);
			if(/precursor/){
				@samples=@line;
			}
            my $sample;
			for(my $index=4;$index <= $#line;$index++){
				$exprs_sample{$samples[$index]}{$line[2]}{$line[0]} = $line[$index];
			}
        }
        close IN;
    }

    my $exdir = "expression_analyses/expression_analyses_${time}";
    
    ## read in mature,precursor and star seuqnces;
    

    open IN,"<$exdir/mature.converted" or die "file mature.converted not found\n";

    my ($seq,$id);
    while(<IN>){
        if(/>(\S+)/){
            $id = $1;
            $seq = <IN>;
            chomp $seq;
            $seq  =lc($seq);
            $seq =~ tr/t/u/;

			if($options{'P'}){
				if($id =~ /-5p/){
					$hm{$id}=$seq;
				}elsif(/-3p/){
					$hs{$id}=$seq;
				}else{
#					print STDERR "option P is used but no 3p or 5p identifier found in $id";
					$hs{$id}=$seq;
					$hm{$id}=$seq;
				}

			}else{
				$hm{$id} = $seq;
			}

#            die "$id\t$seq\n";
        }
    }

	
	my $width = 0;
	for (keys %hm){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
		$width  = length($_) if(length($_) > $width);
	}
	
	$width *=5;

    close IN;
    if(-f "$exdir/star.converted"){
        open IN,"<$exdir/star.converted" or die "file star.converted not found\n";
        
        
        while(<IN>){
            if(/>(\S+)/){
                $id = $1;
                $seq = <IN>;
                chomp $seq;
                $seq  =lc($seq);
                $seq =~ tr/t/u/;
                $hs{$id} = $seq;
            }
        }
        close IN;
    }
    open IN,"<$exdir/precursor.converted" or die "file precursor.converted not found\n";
    while(<IN>){
        if(/>(\S+)/){
            $id = $1;
            $seq = <IN>;
            chomp $seq;
            $hp{$id} = $seq;
        }
    }
    close IN;

    open IN,"<expression_analyses/expression_analyses_${time}/mature2hairpin" or die "Error: File expression_analyses/expression_analyses_${time}/mature2hairpin not found\n";
    my %hairpin2mature2;
	while(<IN>){
        chomp;
		my @line = split(/\t/);
        $hairpin2mature2{$line[0]}=$line[1];
	}
	close IN;
    

    ## read in the expression_analysis signature.arf of mature_sequences mapped to the precursors to get start and end positions
    open IN,"<$options{'i'}" or die "Error:cannot open $options{'i'} file or option -i not given\n\n";
    while(<IN>){
        chomp;
        my @line = split("\t");
        my $id1h = $line[0]; ## this is the mature ID
        my $id2h = $line[5]; ## this is the precursor ID

        ## remove multiple endings if ambigous just for matching with precursor
        $id1h =~ s/\-5p//g;
        $id1h =~ s/\-3p//g;

        ## here is assumed that multiple precursor ids have 3 - in their id, seems to be ok so far
        if($id2h =~/^(\w+\-\w+\-\w+)\-\d+$/){
            $id2h = $1;
        }
        #print STDERR "$id1h\t$line[5]\t$id2h\t$line[0]\n";
        next if($options{'l'} and $id1h !~ /$id2h/i and $id2h !~ /$id1h/i); ## stringent mapping let7a only allowed to map pre-let7a if k is given

        $mat_pre_arf{$line[5]}{$line[0]}{'start'} = $line[7];
        $mat_pre_arf{$line[5]}{$line[0]}{'end'} = $line[8];
        $mat_pre_arf{$line[5]}{$line[0]}{'seq'} = $line[9];
        if($line[0] =~ /\*/){
            $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 1;
        }else{
            $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 0;
        }
    }
    close IN;

    ## think about a stringent option that designates if to show only mappings where the ID of precursor and mature is the same !!!
    ## read in star mapped to mature sequence if a star file was given
    if($options{'j'} and -f $options{'j'}){
        open IN,"<$options{'j'}" or die "Error:cannot open $options{'j'} file or option -j not given\n\n";
        while(<IN>){
            chomp;
            my @line = split("\t");
            my $id1h = $line[0]; ## this is the mature ID
            my $id2h = $line[5]; ## this is the precursor ID

            ## remove multiple endings if ambigous just for matching with precursor
            $id1h =~ s/\-5p//g;
            $id1h =~ s/\-3p//g;
            
            ## here is assumed that multiple precursor ids have 3 - in their id, seems to be ok so far
            #if($id2h =~/^(\w+\-\w+\-\w+)\-\d+$/){
            #    $id2h = $1;
            #}
            next if($options{'l'} and $id1h !~ /$id2h/i and $id2h !~ /$id1h/i); ## stringent mapping let7a only allowed to map pre-let7a if k is given

           # print STDERR "$id2h\t$id1h\t$line[5]\t$line[0]\n";

            $mat_pre_arf{$line[5]}{$line[0]}{'start'} = $line[7];
            $mat_pre_arf{$line[5]}{$line[0]}{'end'} = $line[8];
            $mat_pre_arf{$line[5]}{$line[0]}{'seq'} = $line[9];
            if($line[0] =~ /\*/){
                $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 1;
                if(not $mat_pre_arf{$line[5]}{'sb'}){
                    $mat_pre_arf{$line[5]}{$line[0]}{'sb'} = $line[7];
                }
            }else{
                $mat_pre_arf{$line[5]}{$line[0]}{'star'} = 0;
                if(not $mat_pre_arf{$line[5]}{'mb'}){
                    $mat_pre_arf{$line[5]}{'mb'} = $line[7];
                }
            }
        }
        close IN;
    }


    ## open miRBase.mrd from quantifier module ## precursor ids as hash
	my $oid;
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
        elsif(/^remaining read count\s*(\d+)/){
            $hash_q{$id}{'remaining_rc'}=$1;
            
        }
        
        elsif(/^total read count\s*(\d*)/){
            $hash_q{$id}{"freq_total"}=$1;
        }
        
        ## read in here everything that mapped to the precursor
        elsif(/^(\S+) read count\s*(\d+)/){ ## everything else is just read in as id with real name
            $hash_q{$id}{'mapped'}{$1} = $2;
#			print STDERR "yeah $1\n";
        }
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

            ## this can also be done by just using the information of the mature_signature_arf file
            ## now set labels for mature and star seq in explanation string
            for(my $i=0; $i < length($1); $i++){
                if($d[$i] ne "f"){ 
                    $hash_q{$id}{"ucsc_seq"} .= $s[$i];
                }
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
        }
        elsif(/^(\S\S\S)(\S+)\s+(\S+)\s+(\S+)$/ and $reads){  ## do not read in known miRNAs
        $counter++;
        $hash_q{$id}{"reads"}{$1}{$counter}{"rid"} = "$1$2";
        $hash_q{$id}{"reads"}{$1}{$counter}{"seq"} = $3;
        $hash_q{$id}{"reads"}{$1}{$counter}{"mm"}  = $4;
        }else{}
    }
    close IN;
    print STDERR "parsing miRBase.mrd file finished\n";



    my $sf;
    my $mf;

    for my $id(sort {$hash_q{$b}{"freq_total"} <=> $hash_q{$a}{"freq_total"}} keys %exprs2){
        my %mature = (); ## make hash for mature and star
        my %star = ();
		my $ind;


        my $s_star= $hash_q{$id}{'star_seq'};
            ## now go over everything that is expressed
      
		for my $k2(keys %{$exprs2{$id}}){ ##k1 == precursor, k2 == mirna mapped to it;

			if($options{'P'}){
				if($k2 =~ /-3p/){
					$star{$k2} = $exprs2{$id}{$k2};
				}elsif($k2 =~ /-5p/){
					$mature{$k2} = $exprs2{$id}{$k2};
				}else{
					$ind= index($hash_q{$id}{"pri_seq"},$hash_q{$id}{"mat_seq"});
#					print STDERR "===  $ind  $id $k2  ", length($hash_q{$id}{"pri_seq"})/2,"\n";

					if($ind >  length($hash_q{$id}{"pri_seq"})/2 ){
						$star{$k2} = $exprs2{$id}{$k2};
					}elsif($ind >0){
						$mature{$k2} = $exprs2{$id}{$k2};
					}else{
						print STDERR "Could not determine where $k2 sits in precursor\nPutting it to the 5p species hash\n";
						$mature{$k2} = $exprs2{$id}{$k2};
					}
				}
			}else{
				if($k2 =~ /\*/){
					$star{$k2} = $exprs2{$id}{$k2};
				}else{
					$mature{$k2} = $exprs2{$id}{$k2};
				}
			}

        }




  # for my $k2(keys %{$exprs2{$id}}){ ##k1 == precursor, k2 == mirna mapped to it;
  #           if($k2 =~ /\*/){
  #               $star{$k2} = $exprs2{$id}{$k2};
  #           }else{
  #               $mature{$k2} = $exprs2{$id}{$k2};
  #           }
            
  #       }
        #$id = $mature2hairpin{$oid}; ## this is now the name of the precursor
#            if($id =~ /-1224/){print "$id\n";}
            $hash_q{$id}{"pdf"} = 1;
            $sig++;
            ## set blat link for pre-miRNA sequence
            if($org eq ""){
                $blat = '<td> - </td>';
            }else{
                $blat = "<td><a href=\"http://genome.ucsc.edu/cgi-bin/hgBlat?org=$org&type=BLAT's guess&userSeq=$hash_q{$id}{'pri_seq'}\" target=\"_blank\">blat</a></td>";
            }
            
            $s_star=$hash_q{$id}{'star_seq_obs'} if($hash_q{$id}{'star_seq_obs'});
            my $s_mat = $hash_q{$id}{'mat_seq'};

            ##here the belonging precursor is shown
            $known="<td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$oid</a></td>";
       
            $sf = $hash_q{$id}{"freq_star"};
            $mf = $hash_q{$id}{"freq_mature"};

            if($hash_q{$id}{"freq_mature"} ne $exprs{$oid} ){
                $mf = $exprs{$oid};
                if($sf){
                    $sf = $hash_q{$id}{"freq_mature"};
                }
            }
            

		## do not print precursors with 0 reads mapped to it
		next if(not $hash_q{$id}{"freq_total"});
            ## print miRBase id instead of precursor
            if($options{'d'} or not $hash_q{$id}{"freq_total"}){ ## no link to pdf
                print HTML "<tr><td nowrap=\"nowrap\"><div style=\"visibility: hidden\">$id.pdf </div>$id</a></td>";
            }else{
                print HTML "<tr><td nowrap=\"nowrap\"><a href=\"pdfs_$time/$id.pdf\">$id</a></td>";
            }
            print HTML <<EOF;

            <td>$hash_q{$id}{"freq_total"}</td>            
            <td>

EOF

$mf = "<table>";
$s_mat = "<table>";

#$mf .= "<tr><td> miRNA </td><td>total</td>";
$mf .= "<tr><td WIDTH=$width>id</td>";		
for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
    $mf .= "<td>$sample</td>";
}

if(scalar keys %mature != 0){
for my $mx(keys %mature){ 
#    $mf .= "<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$mx</a></td><td> $mature{$mx} </td>";
$mf .= "\n<tr><td  nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$mx</a></td>";
    for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
        $mf .= "<td><nobr>$exprs_sample{$sample}{$id}{$mx}</nobr></td>";
#		die $exprs_sample{$sample}{$id}{$mx};
    }

    $mf .= "</tr>\n";
    $s_mat .= "<tr><td>$hm{$mx}</td></tr>";
}
}else{
$mf .= "\n<tr><td  nowrap=\"nowrap\"><a>na</a></td>";
for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
        $mf .= "<td><nobr>0</nobr></td>";
#		die $exprs_sample{$sample}{$id}{$mx};
    }

    $mf .= "</tr>\n";
    $s_mat .= "<tr><td>-</td></tr>";
}

$mf .= "</table>\n";
$s_mat .= "</table>";

print HTML "
$mf
</td>
\n";

#		if((scalar keys %star) > 0){
$s_star = "<table>";
$sf = "<table>";
#$sf .= "<tr><td> miRNA </td><td>total    </td>";
$sf .= "<tr><td WIDTH=$width>id</td>";		
for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
    $sf .= "<td>$sample</td>";
}


if((scalar keys %star) != 0){
for my $sx(keys %star){
#    $sf .= "<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$sx</a></td><td>    $star{$sx} </td>";
 
$sf .= "\n<tr><td nowrap=\"nowrap\"><a href=\"http://www.mirbase.org/cgi-bin/query.pl?terms=$id\" target=\"_blank\">$sx</a></td>";

    for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
        $sf .= "<td><nobr>$exprs_sample{$sample}{$id}{$sx}</nobr></td>";
    }
    $sf .="</tr>\n";
    
   $s_star .= "<tr><td>$hs{$sx}</td>\n</tr>";
}
}else{
		$sf .= "\n<tr><td  nowrap=\"nowrap\"><a>na</a></td>";
		for my $sample(sort keys %exprs_sample){                                            #$exprs_sample{$sample}{$line[0]} = $line[1];
			$sf .= "<td><nobr>0</nobr></td>";
#		die $exprs_sample{$sample}{$id}{$mx};
		}

		$sf .= "</tr>\n";
		$s_star .= "<tr><td>-</td></tr>";	
}

$sf .= "</table>\n";
$s_star .= "</table>";


print HTML "
<td>
$sf
</td>
\n";

print HTML <<EOF;

           
			<td nowrap="nowrap">$hash_q{$id}{"remaining_rc"}</td>
            $blat
            <td><a href=${blast}$hash_q{$id}{'pri_seq'}&JOB_TITLE=$hash_q{$id}{"oid"}${blast_query} target="_blank">blast</a></td>
            <td>$s_mat</td>\n<td>$s_star</td>\n<td>$hash_q{$id}{'pri_seq'}</td>
			</tr>     
EOF
        
}

}## end of function
        

sub Usage{
	print STDERR "\n\n\n\n[usage]\n\tperl make_html.pl -f miRDeep_outfile [options]\n\n";
	print STDERR "[options]\n-v 2\t only output hairpins with score above 2\n";
    print STDERR "-c\t also create overview in excel format.\n";
    print STDERR "-k file\t supply file with known miRNAs\n";
    print STDERR "-s file\t supply survey file if score cutoff is used to get information about how big is the confidence of resulting reads\n\n\n";
    print STDERR "-e \t report complete survey file\n";
    print STDERR "-g \t report survey for current score cutoff\n";
    print STDERR "-r file\t Rfam file to check for already reported small RNA sequences\n\n";
    print STDERR "-q file\t miRBase.mrd file produced by quantifier module\n";
    print STDERR "-x file\t signature.arf file with mapped reads to precursors\n";
    print STDERR "-t spec\t specify the organism from which your sequencing data was obtained\n";
    print STDERR "-u \t print all available UCSC input organisms\n";
    print STDERR "-y \ttimestamp of this run\n";
    print STDERR "-o \tsort signature by sample in pdf file, default is by beginning position\n";
    print STDERR "-d \t do not generate pdfs\n\n\n";
    print STDERR "-a \tprint genomic coordinates of mature sequence (still testing)\n";
    print STDERR "-l \tbe stringent when assigning miRNA-precursor connections like mmu-mir only is assigned to mmu-precursor\n";

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
