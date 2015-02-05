#!/usr/bin/perl

use warnings;
use strict;
use IO::File;
use Getopt::Std;
use Cwd;
use File::Copy;
use File::Path;
use File::Basename;
use Term::ANSIColor;


## generate log file for run, all output will be printed to it

my $version="2.0.0.7";

print "

#####################################
#                                   #
# miRDeep$version                    #
#                                   #
# last change: 10/12/2014           #
#                                   #
#####################################

";

my $usage="

miRDeep2.pl reads genome mappings miRNAs_ref/none miRNAs_other/none precursors/none 2>report.log

This script enacts the miRDeep2 pipeline. The input files are:

reads         deep sequences in fasta format. The identifier should contain a prefix, a running
              number and a '_x' to indicate the number of reads that have this sequence.
              There should be no redundancy in the sequences.
genome        genome contigs in fasta format. The identifiers should be unique.
mappings      file_reads mapped against file_genome. The mappings should be in arf format.
              For details on the format see the documentation.
miRNAs_ref    miRBase miRNA sequences in fasta format. These should be the known mature
              sequences for the species being analyzed.
miRNAs_other  miRBase miRNA sequences in fasta format. These should be the pooled known
              mature sequences for 1-5 species closely related to the species being
              analyzed.

precursors    miRBase miRNA precursor sequences in fasta format. These should be the known precursor
              sequences for the species being analyzed.

The output files produced are:

result.html   a html table giving an overview of novel and known miRNAs detected in the
              data. The table is hyperlinked to pdfs showing the signature and structure of
              each hairpin.
result.csv    spread-sheet format of results.html
survey.csv    spread-sheet of prediction accuracy for all score-cutoffs between -10 and 10.
output.mrd    text output of the reported hairpins.

Options:

-a int        minimum read stack height that triggers analysis. Using this option disables
              automatic estimation of the optimal value and all detected precursors are analyzed

-g int        maximum number of precursors to analyze when automatic excision gearing is used.
              default=50.000, if set to -1 all precursors will be analyzed

-b int        minimum score cut-off for predicted novel miRNAs to be displayed in the overview
              table. This score cut-off is by default 0.

-c            disable randfold analysis

-d            disable pdf generation

-t species    species being analyzed - this is used to link to the appropriate UCSC browser entry

-u            output list of UCSC browser species that are supported and exit

-v            remove directory with temporary files

-o            do not sort aligned reads in pdf files by sample, only used if multiple samples were used as input (see Readme for mor information)

-s file       File with known miRBase star sequences

-r string     Prefix for output file names

-z tag        Additional tag appended to current time stamp

-P            use this switch if mature_ref_miRNAs contain miRBase v18 identifiers (5p and 3p) instead of previous ids from v17

Example of use:

miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs -t Mouse 2>report.log

";

## -q file  miRBase.mrd file from quantifier module to show miRBase miRNAs in data that were not scored by miRDeep2

## hardcoded variables
my $read_align_mismatches = 1;


###############################################################################
##
## Don't change anything below this line
##
###############################################################################




## measuring times
my ($sTime,$eTime,$stime,$etime,$sTimeG, $eTimeG, $stimeG, $etimeG);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings);

($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
$second = "0$second" if($second =~ /^\d$/);
$sTimeG = "$hour:$minute:$second";
$stimeG = time;

print "miRDeep2 started at $sTimeG\n\n\n";


my $ctime=time();

my $time=myTime();

my $scripts=`which miRDeep2.pl`;
$scripts =~ s/miRDeep2.pl//;
$scripts =~ s/\s+//g;

############################################# INPUT FILES #########################################################

# foreach my $e(@ARGV){
# 	if($e =~ /(~)/ or $e =~ /(\s)/){
# 		die "*******\nYour call to miRDeep2.pl has an argument $e with allowed characters or whitespace in its name $1 . Please make sure not to have '~' or any whitespace in your arguments\n******\n\n";
# 	}
		
# }

# exit;

print "#Starting miRDeep2\n";
print STDERR "#Starting miRDeep2
$0 @ARGV\n\n";
print STDERR "miRDeep2 started at $sTimeG\n\n\n";

test_first_argument();

my $command_line = "$0 @ARGV\n";

my $file_reads=shift or die "$usage\n\nError: no reads file specified\n";

my $file_genome=shift or die "$usage\n\nError: no genome file specified\n";

my $file_reads_vs_genome=shift or die "$usage\n\nError: no mapping file specified\n";

## remove pathname if any
my $parsed_arf;
if($file_reads_vs_genome =~ /\/*([a-zA-Z_0-9\.]+)$/){
	$parsed_arf = $1;
}else{
	die "could not match arf file\n";
}

my $warning="\n**********\nThe first three arguments to miRDeep2.pl must be files while arguments 4-6 can be files or must be designated as 'none'. Examples:\n
miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf mautre_ref_miRNAs.fa mature_other_miRNAs.fa  hairpin_ref_miRNAs 

or 

miRDeep2.pl reads.fa genome.fa reads_vs_genome.arf none none none


Please check if the supplied files exist and the command call to miRDeep2.pl is correct
*******************\n";

my $file_mature_ref_this_species=shift or die "$usage\n\nError: no file containing miRNAs of investigating species specified\n
Either specify a valid fasta file or say none\n\n$warning";

my $file_mature_ref_other_species=shift or die "$usage\n\nError: no file containing miRNAs of other species specified\n
Either specify a valid fasta file or say none\n\n$$warning";

my $file_precursors=shift or die "$usage\n\nError: no file containing miRNA precursors of this species specified\n
Either specify a valid fasta file or say none\n\n$warning";
if(-f $file_mature_ref_this_species or $file_mature_ref_this_species eq 'none'){}else{die "$usage\n\nError: no file containing miRNAs of investigating species specified\n
Either specify a valid fasta file or say none\n\n$warning";}

if(-f $file_mature_ref_other_species or $file_mature_ref_other_species eq 'none'){}else{die "$usage\n\nError: no file containing miRNAs of other species specified\n
Either specify a valid fasta file or say none\n\n$warning";}

if(-f $file_precursors or $file_precursors eq 'none'){}else{die "$usage\n\nError: no file containing miRNAs precursors of investigating species specified\n
Either specify a valid fasta file or say none\n\n$warning";}


############################################# GLOBAL VARIABLES ####################################################

my %options=();

getopts("a:b:cdt:uvq:s:z:r:p:g:EP",\%options);

my $max_pres=50000;
$max_pres=$options{'g'} if(defined $options{'g'});

## minimal precursor length, used for precheck of precursor file
my $minpreslen=40;
if($options{'p'}){
	$minpreslen=$options{'p'}
};

my $stack_height_min;

if($options{a}){$stack_height_min=$options{a};}

my $dir;
my $dir_tmp;

if(defined $options{'z'}){$time.=$options{'z'}};

################################################## MAIN ############################################################

make_dir_tmp();

test_installed_binaries();

printUsedParameters();

test_input_files();

rna2dna(); ## this makes u to t, lc to uc in sequence and skips everything behind first whitespace in fasta identifier

parse_mappings();

excise_precursors();

prepare_signature();

fold_precursors();

compute_randfold();

miRDeep_core_algorithm();

perform_controls();

make_survey();

output_results();

make_bed();

extract_sequences_from_results();

remove_dir_tmp();

$etimeG = time - $stimeG;
($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
$second = "0$second" if($second =~ /^\d$/);
$eTimeG = "$hour:$minute:$second";

print "

miRDeep runtime: \n
started: $sTimeG
ended: $eTimeG
total:", int($etimeG / 3600),"h:",int(($etimeG % 3600) / 60),"m:",int($etimeG % 60),"s\n\n";

print STDERR "

miRDeep runtime: \n
started: $sTimeG
ended: $eTimeG
total:", int($etimeG / 3600),"h:",int(($etimeG % 3600) / 60),"m:",int($etimeG % 60),"s\n\n";

open OUT,">>$dir/run_${time}_parameters" or die "parameter file not found\n";
print OUT "miRDeep runtime: \n
started: $sTimeG
ended: $eTimeG
total:", int($etimeG / 3600),"h:",int(($etimeG % 3600) / 60),"m:",int($etimeG % 60),"s\n\n";
close OUT;


exit;


################################################## SUBS ##############################################################


sub test_input_files{

################################
## precheck for each file
################################
    open IN,"<$file_reads";
	my $line=<IN>;
	chomp $line;
    if($line !~ /^>\S+/){
        printErr();
        die "The first line of file $file_reads does not start with '>identifier'
Reads file $file_reads is not a valid fasta file\n\n";
    }
	if($line =~ /\s/){
		printErr();
        die "Reads file $file_reads has not allowed whitespaces in its first identifier\n\n";
	}

    if(<IN> !~ /^[ACGTUNacgtun]*$/){
        printErr();
        die "File $file_reads contains not allowed characters in sequences
Allowed characters are ACGTUN
Reads file $file_reads is not a fasta file\n\n";
    }
    close IN;

    open IN,"<$file_genome";
	$line=<IN>;
	chomp $line;
    if($line !~ />\S+/){
        printErr();
        die "The first line of file $file_genome does not start with '>identifier'
Genome file $file_genome is not a fasta file\n\n";
    }
	if($line =~ /\s/){
		printErr();
        die "Genome file $file_genome has not allowed whitespaces in its first identifier\n\n";
	}

	## get genome ids
	my $tmps=`grep ">" $file_genome`;
	my %genomeids = map { $_ => 1 } split("\n",$tmps);



    if(<IN> !~ /^[ACGTUNacgtun]*$/){
        printErr();
        die "File $file_genome contains not allowed characters in sequences
Allowed characters are ACGTUN
Genome file $file_genome is not a fasta file\n\n";
    }
    close IN;

    open IN,"<$file_reads_vs_genome";
    if(<IN> !~ /^(\S+_x\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+(\d+)\s*([mDIM]*)$/){
        printErr();
        die "Mapping file $file_reads_vs_genome is not in arf format\n
Each line of the mapping file must consist of the following fields
readID_wo_whitespaces  length  start  end read_sequence genomicID_wo_whitspaces length  start   end     genomic_sequence  strand  #mismatches editstring
The editstring is optional and must not be contained
The readID must end with _xNumber and is not allowed to contain whitespaces.
The genomeID is not allowed to contain whitespaces.
";
    }
    close IN;

	## get ids from arf file and compare them with ids from the genome file
	$tmps = `cut -f6 $file_reads_vs_genome|sort -u`;
	foreach my $s(split("\n",$tmps)){
		if(not $genomeids{">$s"}){ die "The mapped reference id $s from file $file_reads_vs_genome is not an id of the genome file $file_genome\n\n";}
	}


    if($file_mature_ref_this_species !~ /none/){
        open IN,"<$file_mature_ref_this_species";
		my $line=<IN>;
		chomp $line;
        if($line !~ />\S+/){
            printErr();
            die "The first line of file  $file_mature_ref_this_species does not start with '>identifier'
miRNA reference this species file $file_mature_ref_this_species is not a fasta file\n\n";
        }
		if($line =~ /\s/){
			printErr();
			die "miRNA reference this species file $file_mature_ref_this_species has not allowed whitespaces in its first identifier\n\n";
		}
		$line=<IN>;
        if($line !~ /^[ACGTUNacgtun]*$/){
			printErr();
            die "File $file_mature_ref_other_species contains not allowed characters in sequences
Allowed characters are ACGTUN
miRNA reference this species file $file_mature_ref_this_species is not a fasta file\n\n";
        }

        close IN;
    }

    if($file_mature_ref_other_species !~ /none/){
        open IN,"<$file_mature_ref_other_species";
		my $line=<IN>;
		chomp $line;
        if($line !~ />\S+/){
            printErr();
            die "The first line of file  $file_mature_ref_other_species does not start with '>identifier'
miRNA reference other species file $file_mature_ref_other_species is not a fasta file\n\n";
        }
		if($line =~ /\s/){
			printErr();
			die "miRNA reference other species file $file_mature_ref_other_species has not allowed whitespaces in its first identifier\n\n";
		}
        if(<IN> !~ /^[ACGTUNacgtun]*$/){
            printErr();
            die "File $file_mature_ref_other_species contains not allowed characters in sequences
Allowed characters are ACGTUN
miRNA reference other species file $file_mature_ref_other_species is not a fasta file\n\n";
        }
        close IN;
    }

    if($file_precursors !~ /none/){
        open IN,"<$file_precursors";
		my $line=<IN>;
		chomp $line;
        if($line !~ />\S+/){
            printErr();
            die "The first line of file $file_precursors does not start with '>identifier'
precursor file $file_precursors is not a fasta file\n\n";
        }
		if($line =~ /\s/){
			printErr();
			die "precursor file $file_precursors has not allowed whitespaces in its first identifier\n\n";
		}
		$line=<IN>;
        if($line !~ /^[ACGTUNacgtun]*$/){
            printErr();
            die "File $file_precursors contains not allowed characters in sequences
Allowed characters are ACGTUN
precursor file $file_precursors is not a fasta file\n\n";
        }
		if(length($line) < $minpreslen){
            printErr();
            die "The precursor file $file_precursors does not contain sequences of at least $minpreslen nt
Please make sure that you provided the correct file and the correct parameter ordering when calling miRDeep2.pl
If you have precursors with less than $minpreslen please use option -p <int> to specify this length
\n";

}
        close IN;
    }


##################################################
# precheck finished
##################################################


    #do stringent testing of all input files

    print "#testing input files\n";
    print STDERR "#testing input files\n";

    if($file_mature_ref_this_species !~ /none/){
        start();
        print STDERR "sanity_check_mature_ref.pl $file_mature_ref_this_species\n\n";
        my $ret_file_mature_ref_this_species=`sanity_check_mature_ref.pl $file_mature_ref_this_species 2>&1`;

        if($ret_file_mature_ref_this_species){
            printErr();
            die "problem with $file_mature_ref_this_species\n".$ret_file_mature_ref_this_species;
        }
        end();
    }

    if($file_mature_ref_other_species !~ /none/){
        start();

        print STDERR "sanity_check_mature_ref.pl $file_mature_ref_other_species\n\n";
        my $ret_file_mature_ref_other_species=`sanity_check_mature_ref.pl $file_mature_ref_other_species 2>&1`;

        if($ret_file_mature_ref_other_species){
            printErr();
            die "problem with $file_mature_ref_other_species\n".$ret_file_mature_ref_other_species;
        }
        end();
    }


    print STDERR "sanity_check_reads_ready_file.pl $file_reads\n\n";
    start();
    my $ret_test_file_reads=`sanity_check_reads_ready_file.pl $file_reads 2>&1`;

    if($ret_test_file_reads){
        printErr();
        die "problem with $file_reads\n".$ret_test_file_reads;
    }
    end();
    start();
    print STDERR "sanity_check_genome.pl $file_genome\n\n";
    my $ret_test_file_genome=`sanity_check_genome.pl $file_genome 2>&1`;

    if($ret_test_file_genome){
        printErr();
        die "problem with $file_genome\n".$ret_test_file_genome;
    }
    end();
    start();

    print STDERR "sanity_check_mapping_file.pl $file_reads_vs_genome\n\n";
    my $ret_test_file_reads_genome=`sanity_check_mapping_file.pl $file_reads_vs_genome 2>&1`;

    if($ret_test_file_reads_genome){
        printErr();
        die "problem with $file_reads_vs_genome\n". $ret_test_file_reads_genome;
    }
    end();

    if($file_precursors !~ /none/){
        start();

        print STDERR "sanity_check_mature_ref.pl $file_precursors\n\n";
        my $ret_file_precursors=`sanity_check_mature_ref.pl $file_precursors 2>&1`;

        if($ret_file_precursors){
            printErr();
            die "problem with $file_precursors\n".$ret_file_precursors;
        }
        end();

        start();
        if($file_mature_ref_this_species !~ /none/i){
            print STDERR "Quantitation of expressed miRNAs in data\n\n\n";

            my $species='';
            if($options{'t'}){
                $species="-t $options{'t'}";
			}

            my $file_star='';
            if($options{'s'}){
				if(-s $options{'s'}){
					$file_star = "-s $options{'s'}";
				}else{
					print STDERR "File $options{'s'} specified by option -s is empty or not found\n";
					$options{'s'}=0;
				}
			}

            print "#Quantitation of known miRNAs in data\n";
	        my $dopt="";
			my $Popt="";
	        if($options{'d'}){$dopt="-d";}
			if($options{'P'}){$Popt="-P";}

            my $quant = "quantifier.pl -p $file_precursors -m $file_mature_ref_this_species  -r $file_reads $file_star $species -y $time -k $dopt $Popt";
            print STDERR $quant,"\n";
            `$quant`;
            $options{'q'} = "expression_analyses/expression_analyses_$time/miRBase.mrd";

            end();
        }else{
            print STDERR "Pre-quantitation is skipped caused by missing file with known miRNAs\n\n\n";
        }

    }else{
        print STDERR "Pre-quantitation is skipped caused by missing file with known precursor miRNAs\n\n\n";
    }

    return;
}




sub make_dir_tmp{

    #make temporary directory
    if(not -d "mirdeep_runs"){
        mkdir("mirdeep_runs");
    }


    $dir="mirdeep_runs/run_$time";

    print STDERR "mkdir $dir\n\n";
    mkdir("$dir");

    $dir_tmp = "$dir/tmp";

    mkdir("$dir_tmp");

    return;
}




sub rna2dna{

##process_input mirna files
    if($file_mature_ref_this_species !~ /none/i){
        start();
        ## copy file
        my ( $file_mature_ref_this_species_tmp, $path0, $extension0 ) = fileparse (  $file_mature_ref_this_species, '\..*' );
        print STDERR "rna2dna.pl $file_mature_ref_this_species > $dir_tmp/$file_mature_ref_this_species_tmp$extension0\n\n";
        my $ret_parse_mature_ref_this_species=`rna2dna.pl $file_mature_ref_this_species > $dir_tmp/$file_mature_ref_this_species_tmp$extension0`;
        ##rename orig file
        $file_mature_ref_this_species = $file_mature_ref_this_species_tmp.$extension0;
    }

    if($file_mature_ref_other_species !~ /none/i){
        ##copy file
        my ( $file_mature_ref_other_species_tmp, $path0, $extension0 ) = fileparse (  $file_mature_ref_other_species, '\..*' );
        print STDERR "rna2dna.pl $file_mature_ref_other_species > $dir_tmp/$file_mature_ref_other_species_tmp$extension0\n\n";
        ##here give file name
        my $ret_parse_mature_ref_other_species=`rna2dna.pl $file_mature_ref_other_species > $dir_tmp/$file_mature_ref_other_species_tmp$extension0`;
        ##rename orig file
        $file_mature_ref_other_species =  $file_mature_ref_other_species_tmp.$extension0;
        end();
    }

	if($file_precursors !~ /none/i){
        ##copy file
        my ( $file_precursors_tmp, $path0, $extension0 ) = fileparse (  $file_precursors, '\..*' );
        print STDERR "rna2dna.pl $file_precursors > $dir_tmp/$file_precursors_tmp$extension0\n\n";
        ##here give file name
        my $ret_parse_precursors=`rna2dna.pl $file_precursors > $dir_tmp/$file_precursors_tmp$extension0`;
        ##rename orig file
        $file_precursors =  $file_precursors_tmp.$extension0;
        end();
    }

    return 0;
}


sub parse_mappings{

    #parse mappings to retain only perfect mappings of reads 18 nt <= length <= 25 nt that map perfectly to five loci or less
    print "#parsing genome mappings\n";
    print STDERR "#parsing genome mappings\n";



	print STDERR "parse_mappings.pl $file_reads_vs_genome -a 0 -b 18 -c 25 -i 5 > $dir_tmp/${parsed_arf}_parsed.arf\n\n";
    start();


    my $ret_parse_mappings=`parse_mappings.pl $file_reads_vs_genome -a 0 -b 18 -c 25 -i 5 > $dir_tmp/${parsed_arf}_parsed.arf`;
    end();

    return 0;
}


sub excise_precursors{

    #excise precursors from the genome

    print "#excising precursors\n";
    print STDERR "#excising precursors\n";

    start();
    my $ret_excise_precursors;

    if($options{a}){

        print STDERR "excise_precursors.pl $file_genome $dir_tmp/${parsed_arf}_parsed.arf $dir_tmp/precursors.coords -a $stack_height_min > $dir_tmp/precursors.fa\n\n";
        $ret_excise_precursors=`excise_precursors.pl $file_genome $dir_tmp/${parsed_arf}_parsed.arf $dir_tmp/precursors.coords -a $stack_height_min > $dir_tmp/precursors.fa`;}


    else{
        print STDERR "excise_precursors_iterative_final.pl $file_genome $dir_tmp/${parsed_arf}_parsed.arf $dir_tmp/precursors.fa $dir_tmp/precursors.coords $max_pres\n";
        $ret_excise_precursors=`excise_precursors_iterative_final.pl $file_genome $dir_tmp/${parsed_arf}_parsed.arf $dir_tmp/precursors.fa $dir_tmp/precursors.coords $max_pres`;

		open OSS,"<$dir_tmp/precursors.fa_stack" or die "No file $dir_tmp/precursors.fa_stack found\n";


		$stack_height_min=<OSS>;
		chomp $stack_height_min;
		close OSS;
    }

    end();

	die "No precursors excised\n" if(-z "$dir_tmp/precursors.fa" or not -f "$dir_tmp/precursors.fa");

    return 0;
}



sub prepare_signature{


    #prepare signature file:

    print "#preparing signature\n";
    print STDERR "#preparing signature\n";

    if($file_mature_ref_this_species !~ /none/i){

        print STDERR "prepare_signature.pl $file_reads $dir_tmp/precursors.fa $read_align_mismatches -a $dir_tmp/$file_mature_ref_this_species -o $dir_tmp/signature.arf 2>>error_${time}.log\n\n";
        start();
        my $ret_prepare_signature=`prepare_signature.pl $file_reads $dir_tmp/precursors.fa $read_align_mismatches -a $dir_tmp/$file_mature_ref_this_species -o $dir_tmp/signature.arf 2>>error_${time}.log`;
        end();

    }else{

        print STDERR "prepare_signature.pl $file_reads $dir_tmp/precursors.fa $read_align_mismatches -o $dir_tmp/signature.arf 2>>error_${time}.log\n\n";
        start();
        my $ret_prepare_signature=`prepare_signature.pl $file_reads $dir_tmp/precursors.fa $read_align_mismatches -o $dir_tmp/signature.arf 2>>error_${time}.log`;
        end();

    }
    return 0;
}


sub fold_precursors{

    #predicting RNA secondary structures with RNAfold

    print "#folding precursors\n";
    print STDERR "#folding precursors\n";
    print STDERR "RNAfold < $dir_tmp/precursors.fa -noPS > $dir_tmp/precursors.str\n\n";
    start();
	my $ret_fold_precursors=system("RNAfold < $dir_tmp/precursors.fa -noPS > $dir_tmp/precursors.str 2>>error_${time}.log");
	if($ret_fold_precursors){
		$ret_fold_precursors=system("RNAfold < $dir_tmp/precursors.fa --noPS > $dir_tmp/precursors.str");
		if($ret_fold_precursors){
			die "Some RNAfold error occurred. Error $ret_fold_precursors\n";
		}
	}

    end();
    return;
}



sub compute_randfold{

    if($options{c}){return;}

    #compute randfold p-values for the subset of precursors which are plausible Dicer substrates

    print "#computing randfold p-values\n";
    print STDERR "#computing randfold p-values\n";
    print STDERR "select_for_randfold.pl $dir_tmp/signature.arf $dir_tmp/precursors.str > $dir_tmp/precursors_for_randfold.ids\n\n";
    start();
    my $ret_select_for_randfold=`select_for_randfold.pl $dir_tmp/signature.arf $dir_tmp/precursors.str > $dir_tmp/precursors_for_randfold.ids`;
    end();
    start();
    print STDERR "fastaselect.pl $dir_tmp/precursors.fa $dir_tmp/precursors_for_randfold.ids > $dir_tmp/precursors_for_randfold.fa\n\n";
    my $ret_fasta_select=`fastaselect.pl $dir_tmp/precursors.fa $dir_tmp/precursors_for_randfold.ids > $dir_tmp/precursors_for_randfold.fa`;
    end();
    start();
    print STDERR "randfold -s $dir_tmp/precursors_for_randfold.fa 99 > $dir_tmp/precursors_for_randfold.rand\n\n";
    my $ret_randfold=`randfold -s $dir_tmp/precursors_for_randfold.fa 99 > $dir_tmp/precursors_for_randfold.rand`;
    end();
    return;
}



sub miRDeep_core_algorithm{

    #run miRDeep core algorithm

    print "#running miRDeep core algorithm\n";
    print STDERR "#running miRDeep core algorithm\n";
    my $line;

	my $longest_id=40;
	if($file_mature_ref_this_species !~ /none/i){
		$longest_id= get_longest_id("$dir_tmp/$file_mature_ref_this_species");
	}


    start();

    if($file_mature_ref_other_species !~ /none/i){

        $line="miRDeep2_core_algorithm.pl $dir_tmp/signature.arf $dir_tmp/precursors.str -s $dir_tmp/$file_mature_ref_other_species -v -50 -l $longest_id";

    }else{

        $line="miRDeep2_core_algorithm.pl $dir_tmp/signature.arf $dir_tmp/precursors.str -v -50 -l $longest_id";
    }

    unless($options{c}){$line.=" -y $dir_tmp/precursors_for_randfold.rand";}

    print STDERR "$line > $dir/output.mrd\n";
    my $ret_miRDeep_core=`$line > $dir/output.mrd`;
	if($options{'E'}){
		$ret_miRDeep_core=`$line -t > $dir/error.output.mrd`;
	}

    end();

	#check if file is empty
	if(-z "$dir/output.mrd"){
		print STDERR "Error:\n\tFile $dir/output.mrd is empty\n\n";
		print STDERR "Now running miRDeep2_core_algorithm.pl with option -t to see why all precursors were discarded\n";
		my $ret_miRDeep_core=`$line -t > error.output.mrd_$time`;
		print STDERR "The debug file is called error.output.mrd_$time\n";
		die "\nExiting now\n\n";
	}

    return;
}


sub perform_controls{


    #run permuted controls:
    print "#running permuted controls\n";
    print STDERR "#running permuted controls\n";
    start();
    my $line;
    if($file_mature_ref_other_species !~ /none/i){

        $line="miRDeep2_core_algorithm.pl $dir_tmp/signature.arf $dir_tmp/precursors.str -s $dir_tmp/$file_mature_ref_other_species -v -50";

    }else{

        $line="miRDeep2_core_algorithm.pl $dir_tmp/signature.arf $dir_tmp/precursors.str -v -50";

    }

    unless($options{c}){$line.=" -y $dir_tmp/precursors_for_randfold.rand";}

    print STDERR "echo '$line > $dir/output.mrd' > $dir_tmp/command_line\n\n";
    my $ret_command_line=`echo '$line > $dir/output.mrd' > $dir_tmp/command_line`;
    print STDERR "perform_controls.pl $dir_tmp/command_line $dir_tmp/precursors.str 100 -a > $dir_tmp/output_permuted.mrd 2>>error_${time}.log\n\n";
    my $ret_perform_controls=`perform_controls.pl $dir_tmp/command_line $dir_tmp/precursors.str 100 -a > $dir_tmp/output_permuted.mrd 2>>error_${time}.log`;
    end();
    return;
}



sub make_survey{

    #get overview of the output:

    print "#doing survey of accuracy\n";
    print STDERR "#doing survey of accuracy\n";

    if($file_mature_ref_this_species !~ /none/i){

        print STDERR "survey.pl $dir/output.mrd -a $dir_tmp/output_permuted.mrd -b $dir_tmp/$file_mature_ref_this_species -c $dir_tmp/signature.arf -d $stack_height_min > $dir/survey.csv\n\n";
        start();
        my $ret_survey=`survey.pl $dir/output.mrd -a $dir_tmp/output_permuted.mrd -b $dir_tmp/$file_mature_ref_this_species -c $dir_tmp/signature.arf -d $stack_height_min > $dir/survey.csv`;
        end();

    }else{

        print STDERR "survey.pl $dir/output.mrd -a $dir_tmp/output_permuted.mrd -d $stack_height_min > $dir/survey.csv\n\n";
        start();
        my $ret_survey=`survey.pl $dir/output.mrd -a $dir_tmp/output_permuted.mrd -d $stack_height_min > $dir/survey.csv`;
        end();
    }

    return;
}



sub output_results{

    #making final results html file:

    print "#producing graphic results\n";
    print STDERR "#producing graphic results\n";
    start();

    ## sort aligned reads in pdf not by sample if option -o is given
    my $sort_by_sample = '-o';
    $sort_by_sample = '' if($options{'o'});

    my $line;

    ## choose file to use for counting miRNAs in data
    my $xopt = "$dir_tmp/signature.arf";
    if(-f "expression_analyses/expression_analyses_$time/miRNA_expressed.csv"){
        $xopt = "expression_analyses/expression_analyses_$time/miRNA_expressed.csv";
    }
	my $sc=0;
	$sc=$options{'b'} if($options{'b'});


	my $OE="";
	if($options{'E'}){ $OE=" -E";}


    if($file_mature_ref_this_species !~ /none/i){

        if($options{'q'}){
            $line="make_html.pl -f $dir/output.mrd -k $dir_tmp/$file_mature_ref_this_species -p $dir_tmp/precursors.coords -s $dir/survey.csv -c -e -q $options{'q'} -x $xopt -r ${scripts}Rfam_for_miRDeep.fa -v $sc -y $time $sort_by_sample $OE";
        }else{
            $line="make_html.pl -f $dir/output.mrd -k $dir_tmp/$file_mature_ref_this_species -p $dir_tmp/precursors.coords -s $dir/survey.csv -c -e -r ${scripts}Rfam_for_miRDeep.fa -v $sc -y $time  $sort_by_sample $OE";
        }
    }else{
        if($options{'q'}){
            $line="make_html.pl -f $dir/output.mrd -p $dir_tmp/precursors.coords -s $dir/survey.csv -c -e -q $options{'q'}  -x $xopt -r ${scripts}Rfam_for_miRDeep.fa -v $sc -y $time $sort_by_sample $OE";
        }else{
            $line="make_html.pl -f $dir/output.mrd -p $dir_tmp/precursors.coords -v $sc -s $dir/survey.csv -c -e -r ${scripts}Rfam_for_miRDeep.fa -y $time $sort_by_sample $OE";
        }
    }

    if($options{t}){$line.=" -t $options{t}";}

	my $dopt="";
	if($options{'d'}){$dopt="-d";}
	



    print STDERR "$line -V $version $dopt\n\n";
    my $ret_make_html=`$line -V $version $dopt`;

    end();
    return;
}



##remove temporary directory
sub remove_dir_tmp{
    if($options{v}){

        print STDERR "rmtree($dir_tmp)\n\n";
        rmtree($dir_tmp);
    }
    return;
}




sub test_first_argument{

    if(not $ARGV[0]){die $usage;}

    if($ARGV[0] eq '-u'){system("make_html.pl -u -y 1"); exit;}

    if($ARGV[0] eq '-h' or $ARGV[0] eq '--help'){die $usage;}

    return;
}



## measuring times
sub start{
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $second = "0$second" if($second =~ /^\d$/);
    $sTime = "$hour:$minute:$second";
    $stime = time;
    print STDERR "started: $sTime\n"
}

sub end{
    $etime = time - $stime;
    ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $second = "0$second" if($second =~ /^\d$/);
    $eTime = "$hour:$minute:$second";

    print STDERR "
ended: $eTime
total:", int($etime / 3600),"h:",int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";
}


sub printUsedParameters{
    open OUT,">$dir/run_${time}_parameters" or die "\nError:\n\tcannot create runtime file to save parameters\n";
	print OUT "Start: $time\n";
    print OUT "Script\t$0\n";
    print OUT "args $command_line\n";

    print OUT "dir_with_tmp_files\tdir_miRDeep2_$time\n";

    my $d=cwd();
    print OUT "dir\t$d\n";
    print OUT "file_reads\t$file_reads\n";
    print OUT "file_genome\t$file_genome\n";
    print OUT "file_reads_vs_genome\t$file_reads_vs_genome\n";
    print OUT "file_mature_ref_this_species\t$file_mature_ref_this_species\n";
    print OUT "file_mature_ref_other_species\t$file_mature_ref_other_species\n";

    if($options{'a'}){print OUT "option{a} =\t$options{'a'}\n";}
    if($options{'b'}){print OUT "option{b} =\t$options{'b'}\n";}
    if($options{'c'}){print OUT "option{c} =\t$options{'c'}\n";}
    if($options{'t'}){print OUT "option{t} =\t$options{'t'}\n";}
    if($options{'v'}){print OUT "option{v} =\tused\n";}
#    if($options{'q'}){print OUT "option{q} =\t$options{'q'}\n";}

    close OUT;
}

sub printErr{
    print STDERR color 'bold red';
    print STDERR "Error: ";
    print STDERR color 'reset';
}


sub myTime{
    my ($sec,$min,$hour,$day,$month,$year) = localtime($ctime);
    $year+=1900;
    $month++;
    my $ret=sprintf "%02d_%02d_%02d_t_%02d_%02d_%02d", $day, $month, $year, $hour, $min, $sec;
    return($ret);
}

sub test_installed_binaries{
	my $stdm = "
If you used the install.pl script make sure that you started a complete new shell window after installation.
If this did not help please restart youer workstation.\n\n";

	my $ret;

	$ret = checkBIN("bowtie --version","version");
	die "Error: \tbowtie not found\nCheck if bowtie is correctly installed and all Pathes were set correctly.\n$stdm" if($ret);

	$ret = checkBIN("RNAfold -h","gamma");
	die "Error: \tRNAfold not found\nCheck if RNAfold is correctly installed and all Pathes were set correctly.\n$stdm" if($ret);

	$ret = checkBIN("randfold","let7");
	die "Error: \trandfold not found\nCheck if randfold is correctly installed and all Pathes were set correctly.\n$stdm" if($ret);

	$ret = checkBIN("perl -e \'use PDF::API2; print \"installed\";\'","installed");
	die "Error: \tPerl PDF::API2 package not found\nCheck if the perl PDF::API2 package is correctly installed and all Pathes were set correctly.\n$stdm" if($ret);

	if(not -f "$scripts/Rfam_for_miRDeep.fa"){
		die "Error:\t Rfam_for_miRDeep.fa not found in your miRDeep2 scripts directory\nPlease copy this file from the miRDeep2 archive to your miRDeep2 scripts directory\n\n";
	}

	return 0;
}

sub make_bed{
	my $res=`mirdeep2bed.pl result_${time}.csv > result_${time}.bed`;
	if(!$res){
		return 0;
	}else{
		print STDERR $res,"\n";
		return 1;
	}
}

sub extract_sequences_from_results{
	my $od="mirna_results_${time}";
	
    my $res=`get_mirdeep2_precursors.pl -r result_${time}.csv -p -d -o $od`;
	my $res1=`get_mirdeep2_precursors.pl -r result_${time}.csv -m -p -d -o $od`;
	my $res2=`get_mirdeep2_precursors.pl -r result_${time}.csv -k -p -d -o $od`;
	if(!$res and !$res1 and !$res2){
		print STDERR "fasta and bed files have been created in subfolder $od\n";
		return 0;
	}else{
		print STDERR $res,"\n";
		print STDERR $res1,"\n";
		print STDERR $res1,"\n";
		return 1;
	}
}

sub checkBIN{
    my ($a,$b) = @_;
	my $e = system("$a 1>$dir/tmp/binaries 2>$dir/tmp/binaries2");

    open IN,"<$dir/tmp/binaries";
    my $found = 1;
    while(<IN>){
		if(/$b/){
			$found =0;
		}
	}
	close IN;
	if($found){
		open IN,"<$dir/tmp/binaries2";
		while(<IN>){
			if(/$b/){
				$found =0;
			}
		}
	}
    close IN;
    return $found;
}

sub get_longest_id{
	my ($f) = @_;
	my $l = 0;
	open IN,$f or die "No file given for checking\n";
	while(<IN>){
		if(/>(\S+)/){
			$l = length($1) if(length($1) > $l);
		}
	}
	close IN;
	return $l;
}
