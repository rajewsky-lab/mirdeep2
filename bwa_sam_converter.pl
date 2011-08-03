#!/usr/bin/perl

use strict;

my $usage = "\nError:\nperl bwa_sam_converter.pl mapped.sam reads.fa reads_vs_genome.arf multiple\n
If output_reads.fa is not specified a file called reads_ready.fa is automatically created\n
If mapped.arf is not specified a file called signature.arf is automatically created\n\n";


my $sam = shift or die $usage;
my $rout = shift;
my $arf = shift;
my $multi =shift;

open IN,"<$sam" or die "$usage";

my @line;
my @edit_string;
my @ref_seq;
my $num;

my @processed_seq;
my @processed_edit;

my $rev=0;

my $offset;

my @cigar;          ## the cigar string in sam file refers just to the read sequence
my $print_read;     


if(not $rout){
    open OUT,">reads.fa" or die "cannot create file reads.fa\n";
}else{
    open OUT,">$rout" or die "cannot create file $rout.fa\n";
}

if(not $arf){
    open ARF,">reads_vs_genome.arf" or die "cannot create file reads_vs_genome.arf\n";
}else{
    open ARF,">$arf" or die "cannot create file $arf\n";
}

my $count=0;

my $edit;
my $edit_str;
my $strand;
my $rev = 0;

my $genome_seq;
my $edit_s;

my $cline;

my %reads;

while(<IN>){
	next if(/^\@/);
    $cline = $_;
    $strand = "+";
    $edit = 0;
    $edit_str="";
    @line = split(/\t/);

    ## next if read is not aligned
    next if($line[1] eq 4);

    ## determine if read is coming from minus strand

#    if($line[1] ne 0){
    $rev=FLAGinfo($line[1]);
#    }else{
#        $rev =0;
#    }
    
#    print "$strand\n";
    ## set strand
    $strand = "-" if ($rev);
 #   print "$strand\n";
    ## print ID of read to reads_ready.fa
    #print OUT ">$line[0]\n";
    
    ## grep edit string, this one corresponds to the genome sequence 
    if($cline =~ m/MD:Z:(\S+)\s+/){
        @edit_string = split(//,$1);
    }
    
    @ref_seq = split(//,$line[9]);


    $print_read="";
    
    $offset = 0;
    
    @cigar = split(//,$line[5]);

    $num = "";

    #print "i\tnum\toffset\n";
    for(my $i=0; $i < scalar (@cigar) ; $i++){

        if($cigar[$i] =~ m/\d/){
            $num .= $cigar[$i];

        }elsif($cigar[$i] =~ m/M/){
            $edit_str.= 'm' x $num;
            $print_read .= join(/ /,@ref_seq[$offset..($num-1+$offset)]);
           
            $offset += $num;

            $num="";

        }elsif($cigar[$i] =~ m/I/){
            $offset += $num;
            $edit_str.= 'I';
            $num="";
            $edit++;

        }elsif($cigar[$i] =~ m/D/){
            $edit_str.= 'D';
            $print_read .= ('N' x $num);
            $num="";

		## process N's in Cigar string
        }elsif($cigar[$i] =~ m/N/){
			$num="";
        }else{
		}
    }

    $offset = 0;

    $num="";

    @processed_seq = split(//,$print_read); ## right genome sequence length with N for deleted chars in read sequence
    @processed_edit = split(//,$edit_str);

    ## now process edit string
    for(my $i=0; $i < scalar (@edit_string) ; $i++){
        
        if($edit_string[$i] =~ m/\d/){
            $num .= $edit_string[$i];

        }elsif($edit_string[$i] =~ m/\^/){ ## get the deleted nt in read sequence
            $i++;
            $processed_seq[$num+$offset] = $edit_string[$i]; 
            $edit++;

            $offset+= $num+1;
            $num="";

        }elsif($edit_string[$i] =~ m/\w/){
            $edit++;
            $offset+=$num;

            $processed_seq[$offset] = $edit_string[$i];
            $processed_edit[$offset] = 'M';

            $offset++;
            $num="";
        }else{}
    }
    
    $genome_seq = join("",@processed_seq);
    $edit_s = join("",@processed_edit);


### here reverse if not multi
    if(not $multi){

        if($strand eq "-"){
            $genome_seq = reverse($genome_seq);
            $genome_seq =~ tr/ACGT/TGCA/;
            $line[9] = reverse($line[9]);
            $line[9] =~ tr/ACGT/TGCA/;
            $edit_s = reverse($edit_s);
            
        }
    }

    if($reads{$line[0]} and $reads{$line[0]}{'mm'} > $edit ){
        $reads{$line[0]}{'mm'} = $edit;
        $reads{$line[0]}{'seq'} = $line[9];
    }else{
        $reads{$line[0]}{'mm'} = $edit;
        $reads{$line[0]}{'seq'} = $line[9];
    }

    #print OUT $genome_seq;
    #print OUT "\n";

    if($line[0] !~ /_x\d+/){
        $line[0].="_x1";
    }

    print ARF "$line[0]\t",length($line[9]),"\t1\t",length($line[9]),"\t",lc $line[9],"\t$line[2]\t",length($genome_seq),"\t$line[3]\t",($line[3] -1 + (length($genome_seq))),"\t",lc $genome_seq,"\t$strand\t$edit\t$edit_s\n";
}  
close IN;

for(keys %reads){
    print OUT ">$_\n$reads{$_}{'seq'}\n";
}

close OUT;


sub FLAGinfo{
    my $in=shift;
    my %bwa_codes;

    ## read in bwa hex codes
    while(<DATA>){
        chomp;
        if(/^(\d+)\s+(.+)$/){
            $bwa_codes{$1} = $2;
        }
    }
    
    my $rest;
    
    $rest = $in;
    
    my @arr;
    
    ##modulo operations to determine binary number of decimal number
    while($rest ne 0){
        push(@arr,$rest%2);
    $rest = int($rest/2);
    }
    
    
    
    my $hex;
    my $bin;
    
    my $rev = 0;
    
    ## translate binary to hexadecimal number and check if read is on minus strand
    for(my $i=0; $i < scalar @arr; $i++){
        $bin = $arr[$i] * 2**$i;
        $hex = sprintf("%x", $bin);
        if($arr[$i] ne 0){
            $rev = 1 if($hex eq 10);
        }
    }
    return($rev);
}

__DATA__
0   .
1	the read is paired in sequencing
2	the read is mapped in a proper pair
4	the read sequence is unmapped
8	the mate is unmapped
10	read is mapped to minus strand (given seq in col 10 is therefore the reverse complement of the plus strand)
20	strand of the mate
40	the read is the first read in a pair
80	the read is the second read in a pair
100	the alignment is not primary
200 the read fails plattform/vendor quality checks
400 the read is either a PCR duplicate or an optical duplicate
