#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use File::Copy;
use File::Path;
use Term::ANSIColor;


####################################### USAGE ####################################################

my $usage =
"$0 input_file_reads

This script takes as input a file with deep sequencing reads (these can be in
different formats, see the options below). The script then processes the reads
and/or maps them to the reference genome, as designated by the options given.
Options:

Read input file:
-a              input file is seq.txt format
-b              input file is qseq.txt format
-c              input file is fasta format
-e              input file is fastq format
-d              input file is a config file (see miRDeep2 documentation).
                options -a, -b or -c must be given with option -d.

Preprocessing/mapping:
-g              three-letter prefix for reads (by default 'seq')
-h              parse to fasta format
-i              convert rna to dna alphabet (to map against genome)
-j              remove all entries that have a sequence that contains letters
                other than a,c,g,t,u,n,A,C,G,T,U,N
-k seq          clip 3' adapter sequence
-l int          discard reads shorter than int nts, default = 18
-m              collapse reads

-p genome       map to genome (must be indexed by bowtie-build). The 'genome'
                string must be the prefix of the bowtie index. For instance, if
                the first indexed file is called 'h_sapiens_37_asm.1.ebwt' then
                the prefix is 'h_sapiens_37_asm'.
-q              map with one mismatch in the seed (mapping takes longer)

-r int          a read is allowed to map up to this number of positions in the genome
                default is 5 

Output files:
-s file         print processed reads to this file
-t file         print read mappings to this file 

Other:
-u              do not remove directory with temporary files
-v              outputs progress report

-n              overwrite existing files

-o              number of threads to use for bowtie 

Example of use:

$0 reads_seq.txt -a -h -i -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p h_sapiens_37_asm -s reads.fa -t reads_vs_genome.arf -v
";


###################################### INPUT #######################################################

## create a log file for the mapper.pl
## the latest run of mapper will be on top of the log file
##
##

if(-f "mapper.log"){
	`mv mapper.log mapper.log_bak`;
}else{
	`touch mapper.log_bak`;
}

open MAP,">mapper.log_tmp" or die "could not create mapper.log_tmp\n";

my $cdir = `pwd`;

print MAP "current dir:\t$cdir";
print MAP "mapper command:\t$0 @ARGV\n"; 

my $file_reads=shift or die $usage;

if(not -f "$file_reads"){
	die "No config or reads file could be found\n$usage";
}

my $line=$file_reads;

foreach(@ARGV){
    $line.="\t$_";
}

check_line($line);

my %options=();
getopts("abcdeg:hijk:l:mp:qs:t:uvnr:o:",\%options);


`rm $options{'s'}` if(defined $options{'s'} and -f $options{'s'} and $options{'n'});
`rm $options{'t'}` if(defined $options{'t'} and -f $options{'t'} and $options{'n'});

if(not $options{'l'}){$options{'l'} = 18; }

check_options();




#################################### GLOBAL VARIABLES ################################################
my $threads=1;
$threads=$options{'o'} if(exists $options{'o'});


## check number of cores on the system and threads to be used
my $cores=`grep -ic ^processor /proc/cpuinfo`;
if($cores !~ /^\d+$/){
	$cores=`sysctl -n hw.physicalcpu`;
	if($cores !~ /^\d+$/){
		$cores=`sysctl -n hw.logicalcpu`;
	}
}
if($cores !~ /^\d+$/){
	$cores=1;
}

if($threads > $cores){ print STDERR "More threads specified than cores on the system. Reducing the number of threads to $cores\n"; $threads=$cores;}

my $orig_file_reads;

my $mismatches_seed=0;

if($options{q}){$mismatches_seed=1;}

my $prefix_global="seq";

if($options{g}){$prefix_global=$options{g};}

####################################### MAIN ########################################################


my $dir;#=make_dir_tmp();

if($options{d}){

    handle_config_file($file_reads);

}else{

    handle_one_file($file_reads,$prefix_global);
}

remove_dir_tmp();

print MAP "#"x60,"\n\n";

close MAP;

`cat mapper.log_tmp mapper.log_bak > mapper.log`;
`rm mapper.log_tmp mapper.log_bak`;


## get some statistics about mapped reads if options{'s'} and options{'t'} are supplied
if($options{'s'} and $options{'t'}){
	read_stats();
}
######################################### SUBS #####################################################




sub handle_config_file{

    my $file=shift;

    open (FILE, "$file") or die "can not open $file\n";
    while (<FILE>){

	if(/(^\S+)\s+(\S+)\s*.*$/){

	    my $file=$1;
	    my $prefix=$2;

		if(length($file) < length($prefix)){
			$file=$2;
			$prefix=$1;
		}

	    test_prefix($prefix);

		print MAP "\nhandling file \'$file\' with prefix \'$prefix\'\n";

	    if($options{v}){print STDERR "\nhandling file \'$file\' with prefix \'$prefix\'\n";}

	    handle_one_file($file,$prefix);
	}
    }
    close FILE;
    return;
}


sub make_dir_tmp{
    my ($pref)=@_;
    my $ctime=time();
    my ($sec,$min,$hour,$day,$month,$year) = localtime($ctime);
    $year+=1900;
    $month++;
    my $time=sprintf "%02d_%02d_%02d_t_%02d_%02d_%02d", $day, $month, $year, $hour, $min, $sec;

	print MAP "\ntimestamp:\t$time\n\n"; 
	
	my $num=rand(1);
	my $chance=substr($num,2,10);
    #make temporary directory
    my $dir="dir_mapper${pref}_${chance}_$time";
    mkdir $dir;
    return $dir;
}



sub handle_one_file{

    my($file_reads,$prefix)=@_;
 
    my $file_reads_latest=process_reads($file_reads,$prefix);

    if($options{p}){
        my $file_mapping_latest=map_reads($file_reads_latest);
    }
    return;
}



sub process_reads{
   
    my($file_reads_latest,$prefix)=@_;
    $orig_file_reads=$file_reads_latest;
    $dir=make_dir_tmp("_${prefix}_$orig_file_reads"); 
    #parse Solexa to fasta
    if($options{h}){

        ## parse fastq to fasta
        if($options{e}){

			print MAP "parsing fastq to fasta format\n";

            if($options{v}){print STDERR "parsing fastq to fasta format\n";}

			print MAP "fastq2fasta.pl $file_reads_latest > $dir/reads.fa\n";

            my $ret_format=`fastq2fasta.pl $file_reads_latest > $dir/reads.fa`;
            $file_reads_latest="$dir/reads.fa";
        }else{
			
			print MAP "parsing Solexa / Illumina output to fasta format\n";

            if($options{v}){print STDERR "parsing Solexa / Illumina output to fasta format\n";}
            
            my $line="illumina_to_fasta.pl $file_reads_latest";
    
            if($options{b}){$line.=" -a";}
            
			print MAP "$line > $dir/reads.fa\n";

            my $ret_format=`$line > $dir/reads.fa`;
            
            $file_reads_latest="$dir/reads.fa";
        }
    }

   #rna2dna
    if($options{i}){

		print MAP "converting rna to dna alphabet\n";

	if($options{v}){print STDERR "converting rna to dna alphabet\n";}

		print MAP "rna2dna.pl $file_reads_latest > $dir/reads_dna.fa\n";
		
	my $ret_rna2dna=`rna2dna.pl $file_reads_latest > $dir/reads_dna.fa`;
    
	$file_reads_latest="$dir/reads_dna.fa";
    }


    #discard entries that contain non-canonical letters
    if($options{j}){

		print MAP "discarding sequences with non-canonical letters\n";

	if($options{v}){print STDERR "discarding sequences with non-canonical letters\n";}

		print MAP "fastaparse.pl $file_reads_latest -b > $dir/reads_letters.fa 2>$dir/reads_discarded.fa\n";

	my $ret_clip=`fastaparse.pl $file_reads_latest -b > $dir/reads_letters.fa 2>$dir/reads_discarded.fa`;

	$file_reads_latest="$dir/reads_letters.fa";
    }


    #clip 3' adapters
    if($options{k}){

		print MAP "clipping 3' adapters\n";

	if($options{v}){print STDERR "clipping 3' adapters\n";}

		print MAP "clip_adapters.pl $file_reads_latest $options{k} > $dir/reads_clip.fa\n";

	my $ret_clip=`clip_adapters.pl $file_reads_latest $options{k} > $dir/reads_clip.fa`;

	$file_reads_latest="$dir/reads_clip.fa";
    }


    #discard short reads
    if($options{l}){

		print MAP "discarding short reads\n";

	if($options{v}){print STDERR "discarding short reads\n";}

		print MAP "fastaparse.pl $file_reads_latest -a $options{l} > $dir/reads_no_short.fa 2>$dir/reads_too_short\n";

	my $ret_rem_short=`fastaparse.pl $file_reads_latest -a $options{l} > $dir/reads_no_short.fa 2>$dir/reads_too_short.fa`;

	$file_reads_latest="$dir/reads_no_short.fa";
    }


    #collapse reads
    if($options{m}){

		print MAP "collapsing reads\n";

	if($options{v}){print STDERR "collapsing reads\n";}

		print MAP "collapse_reads_md.pl $file_reads_latest $prefix > $dir/reads_nr.fa\n";

	my $ret_collapse=`collapse_reads_md.pl $file_reads_latest $prefix > $dir/reads_nr.fa`;

	$file_reads_latest="$dir/reads_nr.fa";
    }

    #printing reads
    if($options{s}){

	cat_to($file_reads_latest,$options{s});

#	my $ret=`cat $file_reads_latest >> $options{s}`;	
    }

    return($file_reads_latest);
}

sub map_reads{

    my $file_reads_latest=shift;

    #map reads to genome

	print MAP "mapping reads to genome index\n";

    if($options{v}){print STDERR "mapping reads to genome index\n";}
    
    my $file_genome_latest=$options{p};
    
	my $mapping_loc=5;
	if(defined $options{'r'}){
		$mapping_loc=$options{'r'};
	}

	print MAP "bowtie -p $threads -f -n $mismatches_seed -e 80 -l 18 -a -m $mapping_loc --best --strata $file_genome_latest  --al $dir/${orig_file_reads}_mapped --un $dir/${orig_file_reads}_not_mapped  $file_reads_latest $dir/mappings.bwt\n\n";
#bowtie -f -n $mismatches_seed -e 80 -l 18 -a -m $mapping_loc --best --strata $file_genome_latest $file_reads_latest $dir/mappings.bwt\n\n";

    my $ret_mapping=`bowtie -p $threads -f -n $mismatches_seed -e 80 -l 18 -a -m $mapping_loc --best --strata $file_genome_latest  --al $dir/${orig_file_reads}_mapped --un $dir/${orig_file_reads}_not_mapped  $file_reads_latest $dir/mappings.bwt`;
    
    my $file_mapping_latest="$dir/mappings.bwt";
    
	print MAP "convert_bowtie_output.pl $file_mapping_latest > $dir/mappings.arf\n";

    my $ret_parse_to_arf=`convert_bowtie_output.pl $file_mapping_latest > $dir/mappings.arf`;
    
    $file_mapping_latest="$dir/mappings.arf";
    
    #trim unmapped nts in the 3' end

	print MAP "trimming unmapped nts in the 3' ends\n";
	
    if($options{v}){print STDERR "trimming unmapped nts in the 3' ends\n";}
    
	print MAP "parse_mappings.pl $file_mapping_latest -j > $dir/mappings_trim.arf\n\n";

    my $ret_trim=`parse_mappings.pl $file_mapping_latest -j > $dir/mappings_trim.arf`;
    
    $file_mapping_latest="$dir/mappings_trim.arf";
        
    #printing mappings
    if($options{t}){
		cat_to($file_mapping_latest,$options{t});
#	my $ret=`cat $file_mapping_latest >> $options{t}`;
    }

    return($file_mapping_latest);
}


sub remove_dir_tmp{

    #remove temporary directory
    
    unless($options{u}){

		print MAP "remove tmp dir\nrmtree($dir)\n\n";
		
		rmtree($dir);
    }
    return;
}




sub check_options{

 
    my $formats=0;
    
    if($options{a}){$formats++;}
    
    if($options{b}){$formats++;}
    
    if($options{c}){$formats++;}
    
    if($options{e}){$formats++;}
    
    unless($formats==1){die "exactly one input format (-a, -b , -e or -c) must be designated\n";}
    
    
    my $processing_steps=0;
    
    if($options{h}){$processing_steps++;}
    
    if($options{i}){$processing_steps++;}

    if($options{j}){$processing_steps++;}

    if($options{k}){$processing_steps++;}
    
    if($options{l}){$processing_steps++;}
    
    if($options{m}){$processing_steps++;}

    if($options{p}){$processing_steps++;}
    
    unless($processing_steps>0){die "at least one processing/mapping step (-h, -i, -j, -k, -l, -m or -p) must be designated\n";}
    
    
    my $files_output=0;
    
	if(exists $options{'o'}){
		if($options{'o'} =~ /\d+/ and $options{'o'} > 0){}else{

		die "options{'o'} must be a positive integer\n";}
	}

    if($options{s}){$files_output++;}
    
    if($options{t}){$files_output++;}
    
    unless($files_output>0){die "at least one output file (-s or -t) must be designated\n";}
    
    if($options{s} and -f $options{s} and not $options{'n'}){die "file $options{s} already exists\n";}

    if($options{t} and -f $options{t} and not $options{'n'}){die "file $options{t} already exists\n";}
    
    if(($options{a} or $options{b} or $options{e}) and not($options{h})){die "raw illumina output must be parsed to fasta format with options -h\n";}

    if($options{c} and $options{h}){die "input file is already designated as a fasta file, so option -h should not be used\n";}

    if($options{c} and not($options{i} or $options{j} or $options{k} or $options{l} or $options{m} or $options{p})){die "at least one processing/mapping step (-i, -j, -k, -l, -m or -p) must be designated\n";}

    if($options{d} and not($options{a} or $options{b} or $options{c} or $options{e})){die "option -d must be given with option -a, -b, -c or -e \n";}

    if($options{d} and $options{g}){die "option -d and -g are mutually exclusive. If -d is given, the prefixes must be contained in the config file\n";}

    if($options{g}){test_prefix($options{g});}

    if($options{i} and not($options{c} or $options{h})){die "option -i must be used on reads in fasta format or with option{h}\n";}

    if($options{j} and not($options{c} or $options{h})){die "option -j must be used on reads in fasta format or with option{h}\n";}

    if($options{k} and not($options{c} or $options{h})){die "option -k must be used on reads in fasta format or with option{h}\n";}

    if($options{l} and not($options{c} or $options{h})){die "option -l must be used on reads in fasta format or with option{h}\n";}

    if($options{m} and not($options{c} or $options{h})){die "option -m must be used on reads in fasta format or with option{h}\n";}

    if($options{p} and not($options{c} or $options{h})){die "option -p must be used on reads in fasta format or with option{h}\n";}

    if($options{q} and not($options{p})){die "option -q must be given with option -p\n";}

    if($options{s} and not($options{h} or $options{i} or $options{j} or $options{k} or $options{l} or $options{m} or $options{p})){die "at least one processing step (-h, -i, -j, -k, -l, -m or -p) must be designated if processed file should be output (-s)\n";}

    if($options{t} and not($options{p})){die "reads must be mapped (-p) if mappings are to be output (-t)\n";}
   
    if($options{k} and $options{k}=~/^-/){die "please make sure that the adapter sequence designated with the -k option is correct\n";}

    if($options{l} and $options{l}=~/^-/){die "please make sure that the int given with the -l option is correct\n";}

    if($options{p} and $options{p}=~/ebwt$/){die "please make sure that you are using the -p option correctly.\nThe argument given after -p must be the _prefix_ of the bowtie\nindexed files and should not contain 'ebwt'. For instance,\nif the first indexed file is called 'h_sapiens_37_asm.1.ebwt'\nthen the prefix is 'h_sapiens_37_asm'.\n";}

    if($options{p} and $options{p}=~/^-/){die "please make sure that the genome index designated with the -p option is correct\n";}

## added by SM to check if bowtie is installed when reads should be mapped to genome
    if($options{p}){
        my $binst=`bowtie --version 2>&1`;
        if(not $binst){
            printErr();
            die "Bowtie mapping tool not installed.\n 
Please download from http://downloads.sourceforge.net/project/bowtie-bio/bowtie/ the latest version and install it.\n\n"}
    }

    if($options{s} and $options{s}=~/^-/){die "please make sure that the output file designated with the -s option is correct\n";}

    if($options{t} and $options{t}=~/^-/){die "please make sure that the output file designated with the -t option is correct\n";}

    return;
}


sub check_line{

    my $line=shift;

    if($line=~/-h\s+\d/ or $line=~/-h\s+\w/){die "option -h should not be given with an integer or string\n";}

    if($line=~/-i\s+\d/ or $line=~/-i\s+\w/){die "option -i should not be given with an integer or string\n";}
    
    if($line=~/-j\s+\d/ or $line=~/-j\s+\w/){die "option -j should not be given with an integer or string\n";}
    
    if($line=~/-m\s+\d/ or $line=~/-m\s+\w/){die "option -m should not be given with an integer or string\n";}
    
    if($line=~/-q\s+\d/ or $line=~/-q\s+\w/){die "option -q should not be given with an integer or string\n";}

    return;
}


sub test_prefix{

    my $prefix=shift;

    unless($prefix=~/^\w\w\w$/ and not($prefix=~/_/)){

	die "prefix $prefix does not contain exactly three alphabet letters\n";
    }
    return;
}




sub cat_to{
    
    my($file_1,$file_2)=@_;


    open OUT, ">>$file_2" or die "cannot print to $file_2\n";

    open IN, "<$file_1" or die "cannot read from $file_1\n";
    
    while(my $line = <IN>){

	print OUT "$line";
    }

    close IN;

    close OUT;

    return;
}

sub printErr{
    print STDERR color 'bold red';
    print STDERR "\nError: ";
    print STDERR color 'reset';
}

sub read_stats{
	my %hash;
	my $count;
	my %k2;
	my $total;
	
	open IN,"$options{'s'}" or die "No reads file in fasta format given\n";
	while(<IN>){
		if(/^>*((\S\S\S)\S+_x(\d+))/){
			next if($hash{$1});
			$hash{$1} = 1;
			$count+=$3;
			$k2{$2}+=$3;
		}
	}
	close IN;
	my %hash2;
	my $count2;
	my %k22;
	
	print STDERR "Mapping statistics\n";
	open IN, "$options{'t'}" or die "No mapping file given\n";
	while(<IN>){
		if(/^>*((\S\S\S)\S+_x(\d+))/){
			next if($hash2{$1});
			$hash2{$1} = 1;
			$count2+=$3;
			$k22{$2}+=$3;
		}
	}
	
	print STDERR "\n#desc\ttotal\tmapped\tunmapped\t%mapped\t%unmapped\n";
	print STDERR "total: ",$count,"\t",$count2,"\t",$count-$count2,"\t";
	printf STDERR "%.3f\t%.3f\n",$count2/$count,1-($count2/$count);
	foreach(sort keys %k2){
		print STDERR "$_: ",$k2{$_},"\t",$k22{$_},"\t",$k2{$_}-$k22{$_},"\t";
		printf STDERR "%.3f\t%.3f\n",$k22{$_}/$k2{$_},1-($k22{$_}/$k2{$_});
	}
}
