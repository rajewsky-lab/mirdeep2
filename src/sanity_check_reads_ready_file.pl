#!/usr/bin/perl -W

use strict;

my $counter=0;
my %hash_seq;
my %hash_num;

my $id;
my $tag;

my $hint="Please check your file for the following issues:\n
I.  Sequences are allowed only to comprise characters [ACGTNacgtn].
II. Identifiers are not allowed to have withespaces.\n";

## check number of lines in input file file
my $lines = 0;
my $buffer;
open(FILE, $ARGV[0]) or die "Can't open `$ARGV[0]': $!";
while (sysread FILE, $buffer, 4096) {
    $lines += ($buffer =~ tr/\n//);
}
close FILE;
#print Nicenumber($lines/2), " reads in file\n";




if($lines/2 > 5000000){
    my %rhash;
    my $rn = 0;
    my $seq;
    for(my $i=0; $i < 5000000; $i++){
        while($rhash{$rn}){
            $rn = (2*int(rand($lines/2)));
        }
        $rhash{$rn} = 1;
    }
    open IN,"<$ARGV[0]";
    while(<IN>){
        chomp;
        if(/^\>(.+)$/){
            $counter++;
            $id=$1;
            if($id =~ /\s+/){
                die "Error in line ",Nicenumber($counter),": The identifier\n 
$_\n 
contains white spaces\n

Please make sure that none of the identifiers contain whitepaces.
You could run remove_white_space_in_id.pl $ARGV[0] > newfile
This will remove everything from the id line after the first whitespace
";
            }elsif($id !~ /^(\S\S\S)_(\d+)_(x\d+)$/){
                die "Error in line ",Nicenumber($counter),": The identifier\n
$id\n
has to have the format\nname_uniqueNumber_xnumber\n

Please make sure that all identifiers are unique and have the format described above.
";


            }elsif($id =~ /^(\S\S\S)_(\d+)_(x\d+)$/ and $rhash{$counter}){
                $hash_num{$2}++;
                $tag = $1;
            }else{}
        }elsif(/^([A|C|G|T|N|a|c|g|t|n]{17,})$/){
            $seq = $1;
            if($rhash{$counter}){
                if($hash_seq{$tag}{$seq}){
                    die "Error in line ",Nicenumber($counter),": The sequence\n
$1\n
occures at least twice in your reads file.\n
At first it occured at line ",Nicenumber($hash_seq{$seq}),"

Please make sure that your reads file only contains unique sequences.\n
";
                }else{
                    $hash_seq{$tag}{$seq} = $counter;
                    
                }
            }
        }else{
            die "Error in line ",Nicenumber($counter),": Either the sequence\n

$_\n
contains less than 17 characters or contains characters others than [acgtunACGTUN]\n
Please make sure that your file only comprises sequences that have at least 17 characters\n
containing letters [acgtunACGTUN]\n
";
        }
    }
#   print "read in complete file once\nchecking now for duplicate lines\n";
    close IN;
    $counter=0;

    open IN,"<$ARGV[0]";
    while(<IN>){
        
        if(/^\>(\S\S\S)_.+$/){
            $tag = $1;
            $counter++;
        }
        else{
            if($rhash{$counter}){
                
            }elsif(/^([A|C|G|T|N|a|c|g|t|n]{17,})$/){
                $seq = $1;
                if($hash_seq{$tag}{$seq}){
                    die "Error in line ",Nicenumber($counter),": The sequence\n
$1\n
occures at least twice in sample $tag in your reads file.\n
At first it occured at line ",Nicenumber($hash_seq{$tag}{$seq}),"

Please make sure that your reads file only contains unique sequences within each sample.\n
";
                }
            }
        }
    }
    close IN;

## check now if not in random mode
}else{
    while(<>){
        $counter++;
        if(/^\>(\S\S\S)(.+)$/){
            $id="$1$2";
            $tag =$1;
            
            if($id =~ /\s+/){
die "Error in line ",Nicenumber($counter),": The identifier\n 
$_\n 
contains white spaces\n

Please make sure that none of the identifiers contain whitepaces.
You could run remove_white_space_in_id.pl $ARGV[0] > newfile
This will remove everything from the id line after the first whitespace
";
            }elsif($id !~ /^(\S\S\S)_(\d+)_(x\d+)$/){
die "Error in line ",Nicenumber($counter),": The identifier\n
$id\n
has to have the format\nname_uniqueNumber_xnumber\n

Please make sure that all identifiers are unique and have the format described above.
";
            }elsif($id =~ /^(\S\S\S)_(\d+)_(x\d+)$/){
                $hash_num{$2}++;
            }else{}#print "$id\n";}
            
        }elsif(/^([A|C|G|T|U|N|a|c|g|t|u|n]{17,})$/){
            if($hash_seq{$tag}{$1}){
                die "Error in line ",Nicenumber($counter),": The sequence\n
$1\n
occures at least twice in sample $tag in your reads file.\n
At first it occured at line ",Nicenumber($hash_seq{$tag}{$1}),"

Please make sure that your reads file only contains unique sequences within each sample.\n
";
            }else{
                $hash_seq{$tag}{$1} = $counter;   
            }
        }
        else{
            die "Error in line ",Nicenumber($counter),": Either the sequence\n
$_\n
contains less than 17 characters or contains characters others than [acgtunACGTUN]\n
Please make sure that your file only comprises sequences that have at least 17 characters\n
containing letters [acgtunACGTUN]\n
";
        }      
    }

}
exit;

## subroutine to insert each 3 digits a dot 
sub Nicenumber{
    my @numarr=split(/[.,]/,shift);
    
    my $number = $numarr[0];

    my @n = split(//,reverse($number));
    my $res="";
    for(my $j=0; $j < length($number); $j++){
        if($j%3 eq 0 and $j ne 0){
            $res.=".$n[$j]";
        }else{
            $res.="$n[$j]";
        }
    }
    $res=reverse($res);
    $res.=",$numarr[1]" if($numarr[1]);
    return($res);
}
