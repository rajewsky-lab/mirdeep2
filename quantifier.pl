#!/usr/bin/perl

######################
## Author: SM
## Date: 07/12/2011
######################

use File::Path;
use strict;
use File::Basename;
use Getopt::Std;

my %hash;
my %hash_star;
my %hash_sample;
my %hash_star_sample;
my %total;
my $total_t;

my %seen;
my $species = "none";
my $time= time();


my %organisms;
my %rorganisms;
while(<DATA>){
    chomp;
    if(/^(\S+)\s+(\S+)$/){
        $organisms{$1}=$2;
        $rorganisms{$2}=$1;
    }
}


## options
my %options=();
getopts("p:m:r:s:t:y:dokuc:nxg:e:f:vjwT:",\%options);

## number of mismatches when mapping reads to precursors, default one
my $mismatches = 1;
$mismatches = $options{'g'} if(defined $options{'g'});


my $threads=1;
$threads = $options{'T'} if($options{'T'});

my $upstream = 2;
my $downstream = 5;

$upstream = $options{'e'} if(defined $options{'e'});
$downstream = $options{'f'} if(defined $options{'f'});


if($options{'u'}){
    print STDERR "\n\nAllowed species arguments that have an entry at UCSC\n\n";
    for(keys %organisms){
        print STDERR "$_\t$organisms{$_}\n";
    }
    
    die "\n";
}

my $usage="usage:
  \tperl quantifier.pl [options] -p precursor.fa -m mature.fa -r reads.fa -s star.fa -t species -y [timestamp] -d [pdfs] -o [sort] -k [stringent] -c config.txt -g [number of mismatches in reads vs precursor mappings]

[options]

[mandatory parameters]
  \t-u\tlist all values allowed for the species parameter that have an entry at UCSC 

  \t-p precursor.fa\t miRNA precursor sequences from miRBase
  \t-m mature.fa\t miRNA sequences from miRBase
  \t-r reads.fa\t your read sequences

[optional parameters]
  \t-c [file]    config.txt file with different sample ids... or just the one sample id 
  \t-s [star.fa] optional star sequences from miRBase     
  \t-t [species] e.g. Mouse or mmu
  \t             if not searching in a specific species all species in your files will be analyzed
  \t             else only the species in your dataset is considered
  \t-y [time]    optional otherwise its generating a new one
  \t-d           if parameter given pdfs will not be generated, otherwise pdfs will be generated
  \t-o           if parameter is given reads were not sorted by sample in pdf file, default is sorting
  \t-k           also considers precursor-mature mappings that have different ids, eg let7c
  \t             would be allowed to map to pre-let7a
  \t-n           do not do file conversion again
  \t-x           do not do mapping against precursor again
  \t-g [int]     number of allowed mismatches when mapping reads to precursors, default 1
  \t-e [int]     number of nucleotides upstream of the mature sequence to consider, default 2
  \t-f [int]     number of nucleotides downstream of the mature sequence to consider, default 5
  \t-j           do not create an output.mrd file and pdfs if specified\n
  \t-w           considers the whole precursor as the 'mature sequence'
\n";

if(not $options{'p'} or not $options{'r'}){
    die $usage;
}

if(not $options{'w'} and not $options{'m'}){
    die $usage;
}
if($options{'w'}){
    $options{'m'}="$options{'p'}.dummy";
    open IN,"$options{'p'}";
    open OUT,">$options{'p'}.dummy";
    while(<IN>){
      if(/>/){print OUT;
      }else{
        print OUT substr($_,0,18),"\n";
      }
    }
    close IN;
    close OUT;

    $options{'e'}=0;
    $options{'f'}=0;
}




if($options{'t'}){
  $species = $options{'t'};
}

if($rorganisms{$species}){
  $species = $rorganisms{$species};
}

if($options{'y'}){
   $time = $options{'y'}
}

my $opt_d ="";
if($options{'d'}){
   $opt_d = "-d";
}

## sort pdf reads by sample
my $opt_o ='';
if(not $options{'o'}){
 $opt_o = "-o";
}



my ( $name0, $path0, $extension0 ) = fileparse ( $options{'p'}, '\..*' );
my ( $name1, $path1, $extension1 );
my ( $name1, $path1, $extension1 ) = fileparse ( $options{'m'}, '\..*' );# if(not defined $options{'w'});
my ( $name2, $path2, $extension2 ) = fileparse ( $options{'r'}, '\..*' );
my ( $name3, $path3, $extension3 );
if($options{'s'}){
	( $name3, $path3, $extension3 ) = fileparse ( $options{'s'}, '\..*' );
}

my $dir="expression_analyses";

if(not -d $dir){
  mkdir($dir);
}


my $outdir="${dir}/${dir}_${time}";
if(not -d $outdir){
  mkdir($outdir);
}

## check if reads file has correct format by quickly checking the first line
open IN,"<$options{'r'}" or die "File $options{'r'} not found\n";
my $line = <IN>;
if($line !~ /^>\S\S\S_\d+_x\d+/){
         die "\n$options{'r'} ids do not have the correct format

it must have the id line >SSS_INT_xINT\n
SSS is a three letter code indicating the sample origin
INT is just a running number
xINT is the number of read occurrences\n\n

You can use the mapper.pl module to create such a file from a fasta file with 
mapper.pl $options{'r'} -c -m -s $options{'r'}.collapsed

See also the mapper.pl help for more information on preprocessing input files.

";
}
close IN;

my %samples;

print STDERR "getting samples and corresponding read numbers\n\n";

if(not $options{'c'}){
	open IN,"<$options{'r'}";
    while(<IN>){
       next if($_ !~ /^>/);
       if(/^>(\S\S\S)_\S+_x(\d+)$/){
		$samples{$1} += $2;
       }
    }   
    close IN;
}else{
	if(-f $options{'c'}){
		open IN,"<$options{'c'}";
		while(<IN>){
			chomp;
			if(/^(\S+)\s+(\S+)$/){
             open FIN,"<$1" or die "file $1 not found\n";
             while(<FIN>){
                next if($_ !~ /^>/);
                if(/^>(\S\S\S)_\S+_x(\d+)$/){
		           $samples{$1} += $2;
                }
             }   
             close FIN;
			}
		}
		close IN;
	}else{
		$samples{$options{'c'}}++;
	}
}

for(keys %samples){
   print STDERR "$_\t$samples{$_} reads\n";
}
print STDERR "\n\n";



##convert input files to bowtie accepting format
if(not $options{'n'}){
print STDERR "Converting input files\n";
ConvertFastaFile($options{'p'},$name0,'precursor',$species);
ConvertFastaFile($options{'m'},$name1,'mature',$species);
ConvertFastaFile($options{'r'},$name2,"","");
if($options{'s'}){
	ConvertFastaFile($options{'s'},$name3,'star',$species);
}
}
if(not $options{'x'}){
chdir($outdir);

Mapping();
}else{
chdir($outdir);
}
##now analyze expression file
print STDERR "analyzing data\n";
ReadinPrecursorFile(); 

ReadinMatureMappingFile();

if($options{'s'}){
	ReadinStarMappingFile(); 
}
ReadinReadsMappingFile();
chdir("../../"); 

PrintExpressionValues();
PrintExpressionValuesSamples();

print STDERR "\nCreating miRBase.mrd file\n\n";
die "exiting here and not creating mirdeep.mrd file\nif you want this created do not specify option -j\n" if($options{'j'});
CreateOutputMRD();
#CreateOutputMRD_orig();

my $opt_l ='-l';

if($options{'k'}){
$opt_l = '';
}

my $opt_m='';
if($options{'t'}){
    if($organisms{$options{'t'}}){
        $opt_m = "-m $options{'t'}";
    }elsif($rorganisms{$options{'t'}}){
        $opt_m = "-m $rorganisms{$options{'t'}}";
    }else{$opt_m = "-m $options{'t'}";}
}

my $t;
my $command='';

my $starf='';
if($options{'s'}){$starf ="-j $outdir/${name3}_mapped.arf";}

if( $organisms{$species} ){
   $command = "make_html2.pl -q $outdir/miRBase.mrd -k $name1$extension1 -z -t $organisms{$species} -y $time $opt_d $opt_o -i $outdir/${name1}_mapped.arf $starf $opt_l $opt_m miRNAs_expressed_all_samples_$time.csv";


#  for my $sample(keys %hash_sample){
#         $command .= " miRNAs_expressed_${time}_${sample}.csv";
#  }
        


    print STDERR "$command\n";

    $t=`$command`;

}elsif($rorganisms{$species}){

  
     $command="make_html2.pl -q $outdir/miRBase.mrd -k $name1$extension1 -z -t $species -y $time $opt_d $opt_o -i $outdir/${name1}_mapped.arf $starf $opt_l $opt_m miRNAs_expressed_all_samples_$time.csv";



#  for my $sample(keys %hash_sample){
#         $command .= " miRNAs_expressed_${time}_${sample}.csv";
#  }
        


    print STDERR "$command\n";

    $t=`$command`;



}else{
    $command = "make_html2.pl -q $outdir/miRBase.mrd -k $name1$extension1 -z -y $time $opt_d $opt_o -i $outdir/${name1}_mapped.arf $starf $opt_l $opt_m miRNAs_expressed_all_samples_$time.csv";


# for my $sample(keys %hash_sample){
#        $command .= " miRNAs_expressed_${time}_${sample}.csv";
# }

 print STDERR "$command\n";

 $t=`$command`;

}
exit;


######################################
#                                    #
# subroutines                        #
#                                    #
######################################

sub Mapping{
    my $err;
## build bowtie index
    print STDERR "building bowtie index\n";
    $err = `bowtie-build precursor.converted miRNA_precursor`;
    
## map mature sequences against precursors
    print STDERR "mapping mature sequences against index\n";
#    print STDERR "\nbowtie -f -v 0 -a --best --strata --norc miRNA_precursor $name1.converted ${name1}_mapped.bwt\n\n";
    $err = `bowtie -p $threads -f -v 0 -a --best --norc miRNA_precursor mature.converted ${name1}_mapped.bwt`;
    
## map reads against precursors
    print STDERR "mapping read sequences against index\n"; 
    $err=`bowtie -p $threads -f -v $mismatches -a --best --norc miRNA_precursor $name2.converted ${name2}_mapped.bwt`;
    
    if($options{'s'}){
        print STDERR "mapping star sequences against index\n";
        $err = `bowtie -p $threads -f -v 0 -a --best --norc miRNA_precursor star.converted ${name3}_mapped.bwt`;
    }
} 


sub ConvertFastaFile{ 
    my $file = shift; 
    my $ofile= shift;
	my $des = shift;
    my $sp = shift;
    
    open IN,"$file" or die "File $file not found\n"; 
    if($des eq ""){
        open OUT,">$outdir/$ofile.converted" or die "file $outdir/$ofile.converted could not be created\n";
    }else{
        open OUT,">$outdir/$des.converted" or die "file $outdir/$des.converted could not be created\n";
    }
    my $line;
    my $id;
	my $tmpid;
	my $seq;
	my $first = 1;

    my $sp_hits=0;

    while($line = <IN>){
		chomp $line;
		if($line =~ /^(>\S+)\s*(\S*)/){
			$tmpid = $1;

			if(not $first){
				if($sp eq 'none'){

					if($des eq ""){
						if($seq !~ /N/){   ## skip reads that contain an N in the sequence
							print OUT "$id\n$seq\n";
							$sp_hits++;
						}

					}else{
						if($seq !~ /N/){
							print OUT "$id\n$seq\n";
							$sp_hits++;
						}
					}
				}elsif($id =~ /$sp/i){
					if($des eq ""){
						if($seq !~ /N/){   ## skip reads that contain an N in the sequence
							print OUT "$id\n$seq\n";
						$sp_hits++;
						}
					}else{
						if($seq !~ /N/){
							print OUT "$id\n$seq\n";
							$sp_hits++;
						}
					}
				}else{}

			}else{
				$first = 0;
			}
			$seq="";
			$id = $tmpid;
		}else{
            $line = uc($line);
			$line =~ s/U/T/g;
			$seq .= $line;
		}
	}

 	if($sp eq 'none'){

		if($des eq ""){
			if($seq !~ /N/){   ## skip reads that contain an N in the sequence
				print OUT "$id\n$seq\n";
				$sp_hits++;
			}
		}else{
			if($seq !~ /N/){
				print OUT "$id\n$seq\n";
			}
			$seq="";
			$sp_hits++;
		}
 	}elsif($id =~ /$sp/i){
		if($des eq ""){
			if($seq !~ /N/){   ## skip reads that contain an N in the sequence
				print OUT "$id\n$seq\n";
				$sp_hits++;
			}
		}else{
			if($seq !~ /N/){
				print OUT "$id\n$seq\n";
				$seq="";
				$sp_hits++;
			}
		}
 	}else{}

    close IN;
    close OUT;
    
    if(not $sp_hits){
		die "\nError: No entrys for species \"$options{'t'} or $species\" found in file $file
Please make sure that the given species argument matches the species id in your file $file or say none\n\n\n";
	}
}



sub ReadinPrecursorFile{
    my $id;
    open IN,"precursor.converted" or die "Precursor file precursor.converted not found\n";
    while(<IN>){
        chomp;

        if(/^>(\S+)/){
            $id = $1;
            $hash{$id}{'seq'} = "";
            $hash_star{$id}{'seq'} = "";
            ## make it for different samples now
             for my $sample (keys %samples){
                 $hash_sample{$sample}{$id}{'seq'} = "";
                 $hash_star_sample{$sample}{$id}{'seq'} = "";
                 $hash_sample{$sample}{$id}{'c'} = 0;
                 $hash_sample{$sample}{$id}{'end'} = $hash{$id}{'end'};
             }
        }else{
            $hash{$id}{'seq'} = "$hash{$id}{'seq'}$_"; ## get complete precursor in one line
            $hash_star{$id}{'seq'} = "$hash{$id}{'seq'}$_";
        }
        $hash{$id}{'c'} = 0;
        
		$hash{$id}{'end'} = length($hash{$id}{'seq'});
        
        
        $hash_star{$id}{'c'} = 0;
		$hash_star{$id}{'end'} = length($hash_star{$id}{'seq'});

    }
    close IN;
}


sub ReadinMatureMappingFile{
    my @line;
    my $matches;
    open OUT,">mature2hairpin" or die "cannot create file mature2hairpin\n";
    open IN,"${name1}_mapped.bwt" or die "Mature mapping file ${name1}_mapped.bwt not found \n";
	my $cx;
    my $id1 ='';
    my $id2='';


    while(<IN>){
        $id1= '';
        $id2='';
        @line = split(/\t/);

        $id1 = $line[0]; ## this is the mature ID
        $id2 = $line[2]; ## this is the precursor ID

        ## remove multiple endings if ambigous just for matching with precursor
        $id1 =~ s/\-5p//g;
        $id1 =~ s/\-3p//g;

        ## here is assumed that multiple precursor ids have 3 - in their id, seems to be ok so far
        if($id2 =~/^(\w+\-\w+\-\w+)\-\d+$/){
            $id2 = $1;
        }
        next if(not $options{'k'} and $id1 !~ /$id2/i and $id2 !~ /$id1/i); ## stringent mapping let7a only allowed to map pre-let7a if k is given

		
		$cx++;
		$hash{$line[2]}{'c'}++;           ## how many mature mapped to this precursor
        for my $sample(keys %samples){
            $hash_sample{$sample}{$line[2]}{'c'}++; 
        }
        $matches = $hash{$line[2]}{'c'};

        ## there is a problem, Hash id is from precursor sequence, mature 7a and 7b map to same precursor
        $hash{$line[2]}{$matches}{'beg'} = $line[3]-$upstream;
        $hash{$line[2]}{$matches}{'beg'} = 0 if($hash{$line[2]}{$matches}{'beg'} < 0);
        $hash{$line[2]}{$matches}{'end'} = $line[3]+length($line[4])-1+$downstream;
        $hash{$line[2]}{$matches}{'score'} = 0;
        $hash{$line[2]}{$matches}{'mature'} = $line[0]; ## assign unique mature sequence to precursor

        for my $sample(keys %samples){
            $hash_sample{$sample}{$line[2]}{$matches}{'beg'} = $hash{$line[2]}{$matches}{'beg'};
            $hash_sample{$sample}{$line[2]}{$matches}{'end'} = $hash{$line[2]}{$matches}{'end'};
            $hash_sample{$sample}{$line[2]}{$matches}{'score'} = $hash{$line[2]}{$matches}{'score'};
            $hash_sample{$sample}{$line[2]}{$matches}{'mature'} = $hash{$line[2]}{$matches}{'mature'};
        }


	}
    
	print "\n$cx mature mappings to precursors\n\n";
    close OUT;
    close IN;
}

sub ReadinStarMappingFile{
    my @line;
    my $matches;

    open IN,"${name3}_mapped.bwt" or die "Mature mapping file ${name3}_mapped.bwt not found \n";
	my $cx;
	my $ltmp = "qwertyuiop";
    my $id1 ='';
    my $id2='';


    while(<IN>){
        $id1= '';
        $id2='';
        @line = split(/\t/);

        $id1 = $line[0]; ## this is the mature ID
        $id2 = $line[2]; ## this is the precursor ID

        ## remove multiple endings if ambigous just for matching with precursor
        $id1 =~ s/\*//g;
        $id1 =~ s/\-5p//g;
        $id1 =~ s/\-3p//g;
        if($id1 =~/^(\w+\-\w+\-\w+)\-\d+$/){
            $id1 = $1;
        }
        if($id2 =~/^(\w+\-\w+\-\w+)\-\d+$/){
            $id2 = $1;
        }
        next if(not $options{'k'} and $id1 !~ /$id2/i and $id2 !~ /$id1/i);## maybe this can be removed
		
		$cx++;
		$hash_star{$line[2]}{'c'}++;
        for my $sample(keys %samples){
            $hash_star_sample{$sample}{$line[2]}{'c'}++; 
        }


        #print "$line[2]\t$hash{$line[2]}{'c'}\n";
        $matches = $hash_star{$line[2]}{'c'};
        ## there is a problem, Hash id is from precursor sequence, mature 7a and 7b map to same precursor 
        $hash_star{$line[2]}{$matches}{'beg'} = $line[3]-$upstream;
        $hash_star{$line[2]}{$matches}{'beg'} = 0 if($hash_star{$line[2]}{$matches}{'beg'} < 0);
        $hash_star{$line[2]}{$matches}{'end'} = $line[3]+length($line[4])-1+$downstream;
        $hash_star{$line[2]}{$matches}{'score'} = 0;
        $hash_star{$line[2]}{$matches}{'mature'} = $line[0];


        for my $sample(keys %samples){
            $hash_star_sample{$sample}{$line[2]}{$matches}{'beg'} = $hash_star{$line[2]}{$matches}{'beg'};
            $hash_star_sample{$sample}{$line[2]}{$matches}{'end'} = $hash_star{$line[2]}{$matches}{'end'};
            $hash_star_sample{$sample}{$line[2]}{$matches}{'score'} = $hash_star{$line[2]}{$matches}{'score'};
            $hash_star_sample{$sample}{$line[2]}{$matches}{'mature'} = $hash_star{$line[2]}{$matches}{'mature'};
        }



	}
	print "\n$cx star mappings to precursors\n\n";
}




sub ReadinReadsMappingFile{
    my @line;
    my $rb;
    my $re;
    my @scores;
    my $len_sc;
    open IN,"${name2}_mapped.bwt" or die "Reads mapping File ${name2}_mapped.bwt not found \n";


    
    my $matched = 0;
    my $sample;

    while(<IN>){
        $matched = 0;
        @line = split(/\t/);
        if($species ne "none"){
            next if($line[2] !~ /$species/);
        }

        $rb = $line[3]; 
        $re = ($line[3]+length($line[4])-1);
        
        next if(not $hash{$line[2]}{'c'});

        for(my $i = 1; $i <= $hash{$line[2]}{'c'}; $i++){
			if($options{'w'}){## if consider complete precursor as mature seq
				@scores = split(/x/,$line[0]);
				
				$sample = $1 if($scores[0] =~ /^(\S\S\S)_/); ## get sample id here 
				$len_sc = $scores[$#scores];
				$hash{$line[2]}{$i}{'score'}+= $len_sc; ## hash of pre -> mature -> score
				$hash_sample{$sample}{$line[2]}{$i}{'score'}+= $len_sc;
				$total{$sample}+=$len_sc;
				$total_t+=$len_sc;
#                print "$line[2] ==== $line[0]\t$len_sc\n";
				$matched = 1;

			}else{
				if($rb >= $hash{$line[2]}{$i}{'beg'} and $re <= $hash{$line[2]}{$i}{'end'}){
					@scores = split(/x/,$line[0]);
					
					$sample = $1 if($scores[0] =~ /^(\S\S\S)_/); ## get sample id here 
					$len_sc = $scores[$#scores];
					$hash{$line[2]}{$i}{'score'}+= $len_sc; ## hash of pre -> mature -> score
					$hash_sample{$sample}{$line[2]}{$i}{'score'}+= $len_sc;
					$total{$sample}+=$len_sc;
					$total_t+=$len_sc;
#                print "$line[2] ==== $line[0]\t$len_sc\n";
					$matched = 1;
					
				}
			}
        }


       for(my $i = 1; $i <= $hash_star{$line[2]}{'c'}; $i++){

		   if($options{'w'}){
			   @scores = split(/x/,$line[0]);
			   $sample = $1 if($scores[0] =~ /^(\S\S\S)_/); ## get sample id here 
			   $len_sc = $scores[$#scores];
			   $hash_star{$line[2]}{$i}{'score'}+= $len_sc;          
			   $hash_star_sample{$sample}{$line[2]}{$i}{'score'}+= $len_sc;
			   $matched = 1;
		   }else{

			   if($rb >= $hash_star{$line[2]}{$i}{'beg'} and $re <= $hash_star{$line[2]}{$i}{'end'}){
				   
				   @scores = split(/x/,$line[0]);
				   $sample = $1 if($scores[0] =~ /^(\S\S\S)_/); ## get sample id here 
				   $len_sc = $scores[$#scores];
				   $hash_star{$line[2]}{$i}{'score'}+= $len_sc;          
				   $hash_star_sample{$sample}{$line[2]}{$i}{'score'}+= $len_sc;
				   $total{$sample}+=$len_sc;
				   $total_t+=$len_sc;
				   $matched = 1;
				   
			   }
		   }
        }
        if(not $matched){
            @scores = split(/x/,$line[0]);
            $sample = $1 if($scores[0] =~ /^(\S\S\S)_/); ## get sample id here 
            $hash{$line[2]}{'r'} += $scores[$#scores];
            $hash_sample{$sample}{$line[2]}{'r'}+= $scores[$#scores];
        }
    }
}

sub PrintExpressionValues{
    my $mat;

    open OUT1,">$outdir/miRNA_expressed.csv";
#    open OUT1B,">miRNAs_expressed_$time.csv";
    open OUT2,">$outdir/miRNA_not_expressed.csv";
    print OUT1 "#miRNA\tread count\tprecursor\n";
#    print OUT1B "#miRNA\tread count\tprecursor\n";
    print OUT2 "#miRNA\tread count\n";

	my %seen;
	my %not_seen;


  
    ## check which mature sequences have a mapped read and which not;
    for(sort keys %hash){
        if($species ne "none"){
            next if($_ !~ /$species/);
        }
        for(my $i = 1; $i <= $hash{$_}{'c'}; $i++){
            if($hash{$_}{$i}{'score'} > 0){
                print OUT1 "$hash{$_}{$i}{'mature'}\t$hash{$_}{$i}{'score'}\t$_\n";
#                print OUT1B "$hash{$_}{$i}{'mature'}\t$hash{$_}{$i}{'score'}\t$_\n";
            }else{
                print OUT2 "$hash{$_}{$i}{'mature'}\t0\n";
            }
            if($hash{$_}{$i}{'score'} == 0){
                print OUT1 "$hash{$_}{$i}{'mature'}\t$hash{$_}{$i}{'score'}\t$_\n";
#                print OUT1B "$hash{$_}{$i}{'mature'}\t$hash{$_}{$i}{'score'}\t$_\n";
            }
        }

    }

    ## now for the star sequences
    
    for(sort keys %hash_star){
        if($species ne "none"){
            next if($_ !~ /$species/);
        }
        for(my $i = 1; $i <= $hash_star{$_}{'c'}; $i++){
            if($hash_star{$_}{$i}{'score'} > 0){
                print OUT1 "$hash_star{$_}{$i}{'mature'}\t$hash_star{$_}{$i}{'score'}\t$_\n";
#                print OUT1B "$hash_star{$_}{$i}{'mature'}\t$hash_star{$_}{$i}{'score'}\t$_\n";
            }else{
                print OUT2 "$hash_star{$_}{$i}{'mature'}\t0\n";
            }
        }
    }
    print STDERR "Expressed miRNAs are written to $outdir/miRNA_expressed.csv
not expressed miRNAs are written to $outdir/miRNA_not_expressed.csv\n";
    close OUT1;
    close OUT1B;
    close OUT2;
}

sub PrintExpressionValuesSamples{
    

    open OUTG,">miRNAs_expressed_all_samples_$time.csv";
    print OUTG "#miRNA\tread count\tprecursor\ttotal";
    for my $sample(sort keys %hash_sample){
		next if($sample =~ /config/);
        print OUTG "\t$sample";
    }

	for my $sample(sort keys %hash_sample){
		next if($sample =~ /config/);
        print OUTG "\t$sample (norm)";
    }

    print OUTG "\n";

    for(sort keys %hash){
        if($species ne "none"){
            next if($_ !~ /$species/);
        }

		if($options{'w'}){
			my $i = 1;
			print OUTG "$hash{$_}{$i}{'mature'}\t$hash{$_}{$i}{'score'}\t$_\t$hash{$_}{$i}{'score'}";
			for my $sample(sort keys %hash_sample){
				next if($sample =~ /config/);
				if($hash_sample{$sample}{$_}{$i}{'score'} > 0){
					#print OUTG "\t$hash_sample{$sample}{$_}{$i}{'score'}";
					print OUTG "\t$hash_sample{$sample}{$_}{$i}{'score'}";
				}else{
					print OUTG "\t0";
				}
			}
			
			for my $sample(sort keys %hash_sample){
				next if($sample =~ /config/);
				if($hash_sample{$sample}{$_}{$i}{'score'} > 0){
					#print OUTG "\t$hash_sample{$sample}{$_}{$i}{'score'}";
					print OUTG "\t",sprintf("%.2f",$total_t*$hash_sample{$sample}{$_}{$i}{'score'}/$total{$sample});
				}else{
					print OUTG "\t0";
				}
			}




print OUTG "\n";


		}else{
			for(my $i = 1; $i <= $hash{$_}{'c'}; $i++){
				print OUTG "$hash{$_}{$i}{'mature'}\t$hash{$_}{$i}{'score'}\t$_\t$hash{$_}{$i}{'score'}";
				for my $sample(sort keys %hash_sample){
					next if($sample =~ /config/);
					if($hash_sample{$sample}{$_}{$i}{'score'} > 0){
#						print OUTG "\t$hash_sample{$sample}{$_}{$i}{'score'}";
						print OUTG "\t$hash_sample{$sample}{$_}{$i}{'score'}";
					}else{
						print OUTG "\t0";
					}
				}

				
				for my $sample(sort keys %hash_sample){
					next if($sample =~ /config/);
					if($hash_sample{$sample}{$_}{$i}{'score'} > 0){
						#print OUTG "\t$hash_sample{$sample}{$_}{$i}{'score'}";
						print OUTG "\t",sprintf("%.2f",$total_t*$hash_sample{$sample}{$_}{$i}{'score'}/$total{$sample});
					}else{
						print OUTG "\t0";
					}
				}

				print OUTG "\n";
			}
		}

    }
    if($options{'s'}){
		for(sort keys %hash_star){
        if($species ne "none"){
            next if($_ !~ /$species/);
        }
		my $star_keys = scalar keys %hash_star_sample;
		if($star_keys > 0){
		if($options{'w'}){
			my $i =1;
			
			
			print OUTG "$hash_star{$_}{$i}{'mature'}\t$hash_star{$_}{$i}{'score'}\t$_\t$hash_star{$_}{$i}{'score'}";
			for my $sample(sort keys %hash_star_sample){
				next if($sample =~ /config/);
				if($hash_star_sample{$sample}{$_}{$i}{'score'} > 0){
#					print OUTG "\t$hash_star_sample{$sample}{$_}{$i}{'score'}";
					print OUTG "\t$hash_star_sample{$sample}{$_}{$i}{'score'}";
				}else{
					print OUTG "\t0";
				}
				
			}

			for my $sample(sort keys %hash_star_sample){
				next if($sample =~ /config/);
				if($hash_star_sample{$sample}{$_}{$i}{'score'} > 0){
#					print OUTG "\t$hash_star_sample{$sample}{$_}{$i}{'score'}";
					print OUTG "\t",sprintf("%.2f",$total_t*$hash_star_sample{$sample}{$_}{$i}{'score'}/$total{$sample});
				}else{
					print OUTG "\t0";
				}
				
			}



			print OUTG "\n";
		}else{
			for(my $i = 1; $i <= $hash_star{$_}{'c'}; $i++){
				print OUTG "$hash_star{$_}{$i}{'mature'}\t$hash_star{$_}{$i}{'score'}\t$_\t$hash_star{$_}{$i}{'score'}";
				for my $sample(sort keys %hash_star_sample){
					next if($sample =~ /config/);
					if($hash_star_sample{$sample}{$_}{$i}{'score'} > 0){
			#			print OUTG "\t$hash_star_sample{$sample}{$_}{$i}{'score'}";
						print OUTG "\t$hash_star_sample{$sample}{$_}{$i}{'score'}";
					}else{
						print OUTG "\t0";
					}
				}

				for my $sample(sort keys %hash_star_sample){
					next if($sample =~ /config/);
					if($hash_star_sample{$sample}{$_}{$i}{'score'} > 0){
#					print OUTG "\t$hash_star_sample{$sample}{$_}{$i}{'score'}";
						print OUTG "\t",sprintf("%.2f",$total_t*$hash_star_sample{$sample}{$_}{$i}{'score'}/$total{$sample});
					}else{
					print OUTG "\t0";
					}
					
				}


				print OUTG "\n";
			}
		}
	}

    }
}
    close OUTG;

#     for my $sample(sort keys %hash_sample){
# 		next if($sample =~ /config/);
#         my $mat;
#         open OUT,">miRNAs_expressed_${time}_${sample}.csv";
#         print OUT "#miRNA\tread count\tprecursor\n";

#         my %seen;
#         my %not_seen;
        

  
#         ## check which mature sequences have a mapped read and which not;
#         for(sort keys %{$hash_sample{$sample}}){
#             if($species ne "none"){
#                 next if($_ !~ /$species/);
#             }
#             for(my $i = 1; $i <= $hash_sample{$sample}{$_}{'c'}; $i++){
#                 print OUT "$hash_sample{$sample}{$_}{$i}{'mature'}\t$hash_sample{$sample}{$_}{$i}{'score'}\t$_\n";
#             }
#         }
        
#         ## now for the star sequences
#         if($options{'s'}){
#         for(sort keys %{$hash_star_sample{$sample}}){
#             if($species ne "none"){
#                 next if($_ !~ /$species/);
#             }
#             for(my $i = 1; $i <= $hash_star{$_}{'c'}; $i++){
#                 print OUT "$hash_star_sample{$sample}{$_}{$i}{'mature'}\t$hash_star_sample{$sample}{$_}{$i}{'score'}\t$_\n";
#             }
#         }
# 	}

#         close OUT;
#     }
#     print STDERR "Expressed miRNAs are written to miRNA_expressed_sample.csv\n";
}



sub CreateOutputMRD{
    my %exprs;
    my @ex;
    ## everything goes in here stars and other stuff
    open IN,"<$outdir/miRNA_expressed.csv" or die "File $outdir/miRNA_expressed.csv not found";
    while(<IN>){
        chomp;
        next if(/precursorID/);
        @ex = split("\t");
        $exprs{$ex[2]}{$ex[0]} = $ex[1]; ## precursor-entitiy = count
    }
    close IN;

	chdir($outdir);
	open OUT,">miRBase.mrd" or die "could not create file $outdir/miRBase.mrd\n";
	## get mature mappings
	my ($mature,$star,$reads);
	$mature=`convert_bowtie_output.pl ${name1}_mapped.bwt > ${name1}_mapped.arf`;


	my @line;
	my @tmp;

	my $ltmp;

    my $id1= '';
    my $id2='';
	open IN,"<${name1}_mapped.arf";
	while(<IN>){        
        $id1= '';
        $id2='';
        @line = split(/\t/);

        if($species ne "none"){
            next if($line[5] !~ /$species/);
        }


        $id1 = $line[0]; ## this is the mature ID
        $id2 = $line[5]; ## this is the precursor ID

        ## remove multiple endings if ambigous just for matching with precursor
        $id1 =~ s/\-5p//g;
        $id1 =~ s/\-3p//g;

        ## here is assumed that multiple precursor ids have 3 - in their id, seems to be ok so far
        if($id2 =~/^(\w+\-\w+\-\w+)\-\d+$/){
            $id2 = $1;
        }
        next if(not $options{'k'} and $id1 !~ /$id2/i and $id2 !~ /$id1/i); ## if k then stringent checking is on ;may this can be removed

        if(not $hash{$line[5]}{'struct'}){
            @tmp  = split(//,"f" x length($hash{$line[5]}{'seq'}));
        }else{
            @tmp  = split(//,$hash{$line[5]}{'struct'});
        }
## this here gets complicated when two times is a 5p and 3p there 

        if($line[0] =~ /5p/){
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = '5';
            }
        }elsif($line[0] =~ /3p/){
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = '3';
            }
        }else{
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = 'M';
            }   
        }
		$hash{$line[5]}{'struct'} = join('',@tmp);
        
	}
	close IN;

    
    ## now treat the star sequences if given

    if($options{'s'}){
        $star=`convert_bowtie_output.pl ${name3}_mapped.bwt > ${name3}_mapped.arf`;
        open IN,"<${name3}_mapped.arf";
        while(<IN>){
            @line = split(/\t/);
            if($species ne "none"){
                next if($line[5] !~ /$species/);
            }

            $id1= '';
            $id2='';
            @line = split(/\t/);
            
            $id1 = $line[0]; ## this is the mature ID
            $id2 = $line[5]; ## this is the precursor ID

            ## remove multiple endings if ambigous just for matching with precursor
            $id1 =~ s/\*//g;
            $id1 =~ s/\-5p//g;
            $id1 =~ s/\-3p//g;
            if($id1 =~/^(\w+\-\w+\-\w+)\-\d+$/){
                $id1 = $1;
            }
            if($id2 =~/^(\w+\-\w+\-\w+)\-\d+$/){
                $id2 = $1;
            }

            next if(not $options{'k'} and $id1 !~ /$id2/i and $id2 !~ /$id1/i); ## maybe this can be removed
            @tmp  = split(//,$hash{$line[5]}{'struct'});
            
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = 'S';
            }
            $hash{$line[5]}{'struct'} = join('',@tmp);
        }
        close IN;
    }
    
 #    print "after star thingy\n\n";
#     foreach(keys %hash){
#         print ">$_\n$hash{$_}{'struct'}\n";
#     }
#     print "\n\n\n";


    ## create loop entrys in struct if there is a mature and star sequence
    my @loop;
    my $found = 0;
    my $loop =0;
	if($options{'s'}){
		for my $k(keys %hash){
			@loop = split(//,$hash{$k}{'struct'});
			for(my $i=0; $i< scalar @loop; $i++){
				if($loop[$i] ne 'f' and $found == 0){
					$found = 1;
				}elsif($loop[$i] eq 'f' and $found == 1){
					$loop = 1;
				}elsif($loop[$i] ne 'f' and $found == 1 and $loop ==1){
					$found = 2;
				}else{}
			}
        ## set loop right now if $found = 2;
			if($found == 2 ){
				
				$found = 0;
				$loop =0;
				for(my $i=0; $i< scalar @loop; $i++){
					if($loop[$i] ne 'f' and $found == 0){
						$found = 1;
					}elsif($loop[$i] eq 'f' and $found == 1 and $loop <2 ){
						$loop[$i] = 'l';
						$loop = 1;
					}elsif($loop[$i] ne 'f' and $found == 1 and $loop > 0){
						$found = 2;
					}else{}
				}
			}
			$hash{$k}{'struct'} = join('',@loop);
		}
	}


	$reads=`convert_bowtie_output.pl ${name2}_mapped.bwt > ${name2}_mapped.arf`;
	open IN,"<${name2}_mapped.arf";
	my (@rseq,@qseq);
	my $counter; 
	my $col1_width = 40;
	my $spacer;

 
#     print "after loop thingy\n\n";
#     foreach(keys %hash){
#         print ">$_\n$hash{$_}{'struct'}\n";
#     }
#     print "\n\n\n";


    ## process the reads now ## stopped here , check if something else needs to be done here or not
#    die "${name2}_mapped.arf";
    my @rc;
	while(<IN>){
		chomp;
		next if(/^\s*$/);

		@line = split(/\t/);
		if($species ne "none"){
            next if($line[5] !~ /$species/);
        }
		
        @rc = split("_x",$line[0]);
        
        $hash{$line[5]}{'reads'}{'c'}++;
		$hash{$line[5]}{'reads'}{'tot'}+= $rc[1]; ## sum up total read count for this precursor
		$counter = $hash{$line[5]}{'reads'}{'c'};
		
		$hash{$line[5]}{'reads'}{$counter} = $_; ## all lines of mapping in hash now
	}	
    print "after READS READ IN thing\n\n";



	for my $k1(sort keys %hash){
		if($species ne "none"){
			next if($k1 !~ /$species/);
		}
	
	
		$hash{$k1}{'seq'} =~ tr/ACGTU/acguu/;
	
		my @str=`echo $hash{$k1}{'seq'} |RNAfold`;
		
		my @str1 = split(/\s+/,$str[1]);
  
        ## not yet exactly correct because miRNAs with same precursor have different stars
		print OUT ">$k1\n";
        $spacer = " " x ($col1_width - length('total read count'));
        
        ## print total read count to precursor
        print OUT "total read count$spacer",$hash{$k1}{'reads'}{'tot'},"\n"; ### mmm error is here but not seen yet

        ## print all mature ids given in mature file for this precuror
        for my $k2 (keys %{$hash{$k1}}){
            next if($k2 !~ /^\d+$/);    ## skip all keys that are not a number
            my $mat = $hash{$k1}{$k2}{'mature'};
            $spacer = " " x ($col1_width - length("$mat read count"));
            print OUT "$mat read count$spacer$hash{$k1}{$k2}{'score'}\n";
        }
        
         for my $k2 (keys %{$hash_star{$k1}}){
            next if($k2 !~ /^\d+$/);    ## skip all keys that are not a number
            my $mat = $hash_star{$k1}{$k2}{'mature'};
            $spacer = " " x ($col1_width - length("$mat read count"));
            print OUT "$mat read count$spacer$hash_star{$k1}{$k2}{'score'}\n";
        }
        $spacer = " " x ($col1_width - length('remaining read count'));
        print OUT "remaining read count$spacer",$hash{$k1}{'reads'}{'tot'}-$hash{$k1}{1}{'score'}-$hash_star{$k1}{1}{'score'},"\n";


#        $spacer = " " x ($col1_width - length('star read count'));
#        print OUT "star read count$spacer$hash{$k1}{'s'}\n";

		$spacer = " " x ($col1_width - length('exp')); 
		print OUT "exp$spacer$hash{$k1}{'struct'}\n";
		$spacer = " " x ($col1_width - length('pri_seq')); 
		print OUT "pri_seq$spacer$hash{$k1}{'seq'}\n";
		$spacer =   " " x ($col1_width - length('pri_struct'));
		print OUT "pri_struct$spacer$str1[0]\t#MM\n";	
		
        

		my @reads_arr;
		my %reads_hash;
	
		my @pseq = split(//,lc $hash{$k1}{'seq'});
		
		## put all reads of key in an array then use mirdeep2 routine
		for my $k2(keys %{$hash{$k1}{'reads'}}){
			
			next if($k2 eq 'c');
			push(@reads_arr,$hash{$k1}{'reads'}{$k2});
		}
		
		## now use miRDeep2_core routine
		my $lr = scalar @reads_arr;
		my $rrc = 0;
		
		foreach(@reads_arr){
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
			my $rseq  = lc $reads_hash{$_}{'seq'};
			$rseq =~ tr/t/u/;
			my $bef="." x ($reads_hash{$_}{'beg'}-1);
			my $after = "." x ($hash{$k1}{'end'} - $reads_hash{$_}{"end"});
			my $spacer = " " x ($col1_width - length($reads_hash{$_}{'id'})); 
			my @sread = split(//,$rseq);
			
			
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
			print OUT "$reads_hash{$_}{'id'}$spacer$bef$rseq$after\t$reads_hash{$_}{'mm'}\n";
		}
		print OUT "\n\n\n";
		
	} ## close $for my $k1
	close OUT;
	chdir("../../");
}## close sub

sub CreateOutputMRD_orig{
    my %exprs;
    my @ex;
    open IN,"<$outdir/miRNA_expressed.csv" or die "File $outdir/miRNA_expressed.csv not found";
    while(<IN>){
        chomp;
        next if(/precursorID/);
        @ex = split("\t");
        $exprs{$ex[2]} = $ex[1];
    }
    close IN;

	chdir($outdir);
	open OUT,">miRBase.ord" or die "could not create file $outdir/miRBase.mrd\n"; 
	my ($mature,$star,$reads);
	$mature=`convert_bowtie_output.pl ${name1}_mapped.bwt > ${name1}_mapped.arf`;

	my @line;
	my @tmp;

	my $ltmp;


	open IN,"<${name1}_mapped.arf";
	while(<IN>){
#        next if(not $exprs{$line[5]}); ## skipp mapping if not expressed precursor
		@line = split(/\t/);
		if($species ne "none"){
            next if($line[5] !~ /$species/);
        }
        
        $ltmp = "qwertyuuiioo";
		## get all precursor ids # some are of the form let7a-1 let7a-2 etc. 
		if($line[5] =~ /^(\S+)-\d+$/){ ## this is necessary because sometimes exist more than one precursor sequence
			$ltmp = $1;
		}
        next if($line[0] !~ /$line[5]/i and $line[5] !~ /$line[0]/i and $line[0] !~ /$ltmp/i);
        
        if(not $hash{$line[5]}{'struct'}){
            @tmp  = split(//,"f" x length($hash{$line[5]}{'seq'}));
        }else{
            @tmp  = split(//,$hash{$line[5]}{'struct'});
        }


        if($line[0] =~ /5p/){
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = 'M';
            }
        }elsif($line[0] =~ /3p/){
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = 'S';
            }
        }else{
            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = 'M';
            }   
        }
		$hash{$line[5]}{'struct'} = join('',@tmp);
	   
	}
	close IN;

    if($options{'s'}){
        $star=`convert_bowtie_output.pl ${name3}_mapped.bwt > ${name3}_mapped.arf`;
        open IN,"<${name3}_mapped.arf";
        while(<IN>){
            @line = split(/\t/);
            if($species ne "none"){
                next if($line[5] !~ /$species/);
            }
            $ltmp = "qwertyuiop";
            ## get all precursor ids # some are of the form let7a-1 let7a-2 etc. 
            if($line[5] =~ /^(\S+)-\d+$/){ ## this is necessary because sometimes exist more than one precursor sequence
                $ltmp = $1;
            }
            next if($line[0] !~ /$line[5]/i and $line[5] !~ /$line[0]/i and $line[0] !~ /$ltmp/i);


            @tmp  = split(//,$hash{$line[5]}{'struct'});

            for(my $i = $line[7]-1; $i <= $line[8]-1; $i++){
                $tmp[$i] = 'S';
            }

            $hash{$line[5]}{'struct'} = join('',@tmp);

        }
        close IN;
    }
    
    ## create loop entrys in struct if there is a mature and star sequence
    my @loop;
    for my $k(keys %hash){
        my ($m,$s);
        @loop = split(//,$hash{$k}{'struct'});
        for(my $i=0; $i< scalar @loop; $i++){
            if($loop[$i] eq 'M' and not $s){
                $m =1;
            }elsif($loop[$i] eq 'M' and $s){
                $m =2;
            }elsif($loop[$i] eq 'S' and not $m){
                $s = 1;
            }
            elsif($loop[$i] eq 'S' and $m){
                $s =2;
            }else{}
        }
        
        if($m == 1 and $s == 2){
            $m =0;
            $s =0;
            for(my $i=0; $i< scalar @loop; $i++){
                if($loop[$i] eq 'M'){
                    $m = 1;
                }elsif($m and $loop[$i] eq 'f' and not $s){
                    $loop[$i]  ='l';
                }elsif($loop[$i] eq 'S'){
                    $s = 1;
                }else{}
            }
        }elsif($s == 1 and $m == 2){
            $m =0;
            $s =0;
            for(my $i=0; $i< scalar @loop; $i++){
                if($loop[$i] eq 'S'){
                    $s = 1;
                }elsif($s and $loop[$i] eq 'f' and not $m){
                    $loop[$i]  ='l';
                }elsif($loop[$i] eq 'M'){
                    $m = 1;
                }else{}
            }
        }else{}
        $hash{$k}{'struct'} = join('',@loop);
    }

	$reads=`convert_bowtie_output.pl ${name2}_mapped.bwt > ${name2}_mapped.arf`;
	open IN,"<${name2}_mapped.arf";
	my (@rseq,@qseq);
	my $counter; 
	my $col1_width = 40;
	my $spacer;

	while(<IN>){
		chomp;
		next if(/^\s*$/);

		@line = split(/\t/);
		if($species ne "none"){
            next if($line[5] !~ /$species/);
        }
		
		$hash{$line[5]}{'reads'}{'c'}++;
		$counter = $hash{$line[5]}{'reads'}{'c'};
		
		$hash{$line[5]}{'reads'}{$counter} = $_; ## all lines of mapping in hash now
	}	
	   
		
	for my $k1(sort keys %hash){
		if($species ne "none"){
			next if($k1 !~ /$species/);
		}
	
	
		$hash{$k1}{'seq'} =~ tr/ACGTU/acguu/;
	
		my @str=`echo $hash{$k1}{'seq'} |RNAfold`;
		
		my @str1 = split(/\s+/,$str[1]);
  
        ## not yet exactly correct because miRNAs with same precursor have different stars
		print OUT ">$k1\n";
        $spacer = " " x ($col1_width - length('total read count'));
        print OUT "total read count$spacer",($hash{$k1}{'m'}+$hash{$k1}{'s'} + $hash{$k1}{'r'}),"\n"; ### mmm error is here but not seen yet
        $spacer = " " x ($col1_width - length('mature read count'));
        print OUT "mature read count$spacer$hash{$k1}{'m'}\n";
        $spacer = " " x ($col1_width - length('loop read count'));
        print OUT "loop read count$spacer\n";
        $spacer = " " x ($col1_width - length('star read count'));
        print OUT "star read count$spacer$hash{$k1}{'s'}\n";
		$spacer = " " x ($col1_width - length('exp')); 
		print OUT "exp$spacer$hash{$k1}{'struct'}\n";
		$spacer = " " x ($col1_width - length('pri_seq')); 
		print OUT "pri_seq$spacer$hash{$k1}{'seq'}\n";
		$spacer =   " " x ($col1_width - length('pri_struct'));
		print OUT "pri_struct$spacer$str1[0]\t#MM\n";	
		

		my @reads_arr;
		my %reads_hash;
	
		my @pseq = split(//,lc $hash{$k1}{'seq'});
		
		## put all reads of key in an array then use mirdeep2 routine
		for my $k2(keys %{$hash{$k1}{'reads'}}){
			
			next if($k2 eq 'c');
			push(@reads_arr,$hash{$k1}{'reads'}{$k2});
		}
		
		## now use parts of the miRDeep2_core routine code
		my $lr = scalar @reads_arr;
		my $rrc = 0;
		
		foreach(@reads_arr){
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
			my $rseq  = lc $reads_hash{$_}{'seq'};
			$rseq =~ tr/t/u/;
			my $bef="." x ($reads_hash{$_}{'beg'}-1);
			my $after = "." x ($hash{$k1}{'end'} - $reads_hash{$_}{"end"});
			my $spacer = " " x ($col1_width - length($reads_hash{$_}{'id'})); 
			my @sread = split(//,$rseq);
			
			
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
			print OUT "$reads_hash{$_}{'id'}$spacer$bef$rseq$after\t$reads_hash{$_}{'mm'}\n";
		}
		print OUT "\n\n\n";
		
	} ## close $for my $k1
	close OUT;
	chdir("../../");
}## close sub




__DATA__
hsa                 Human
ptr                 Chimp
na                  Orangutan
na                  Rhesus
na                  Marmoset
mmu                 Mouse
rno                 Rat
na                  Guinea Pig
lca                 Cat
cfa                 Dog
eca                 Horse
bta                 Cow
na                  Opossum
na                  Platypus
gga                 Chicken
na                  Zebra finch
na                  Lizard
xtr                 X. tropicalis
tni                 Zebrafish
tni                 Tetraodon
fru                 Fugu
na                  Stickleback
na                  Medaka
na                  Lamprey
bfl                 Lancelet
cin                 C.intestinalis
spu                 S.purpuratus
cel                 C.elegans
na                  C.brenneri
cbr                 C.briggsae
na                  C.remanei
sja                 C.japonica
na                  P.pacificus
dme                 D.melanogaster
dsi                 D.simulans
dse                 D.sechellia
dya                 D.yakuba
der                 D.erecta
dan                 D.ananassae
dps                 D.pseudoobscura
dpe                 D.persimilis
dvi                 D.virilis
dmo                 D.mojavensis
dgr                 D.grimshawi
aga                 A.gambiae
ame                 A.mellifera
na                  S.cerevisiae
