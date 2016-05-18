#!/usr/bin/perl -W

use strict;

my $counter=0;
my %hash_num;

my $id;

my $hint="Please check your file for the following issues:\n
I.  Sequences are allowed only to comprise characters [ACGTNacgtn].
II. Identifiers are not allowed to have withespaces.\n";

while(<>){
    $counter++;
    if(/^\>(.+)$/){
        $id=$1;        
        if($id =~ /\s+/){
die "Error in line ",Nicenumber($counter),": The identifier\n 
$id\n 
contains white spaces\n

$hint

You could run remove_white_space_in_id.pl inputfile > newfile
This will remove everything from the id line after the first whitespace

";
        }else{
            $hash_num{$id}++;
        }
    }elsif(not /^([A|C|G|T|U|N|a|c|g|t|u|n]+)$/){
        die "Error in line ",Nicenumber($counter),": The sequence\n
$_\n
contains characters others than [acgtunACGTUN]\n
$hint
";
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
