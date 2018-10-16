#!/usr/bin/perl -W

# miRDeep2 sanity-check-genome perl script
# Copyright (C) 2008 - 2011  Marc Friedländer
# Copyright (C) 2009 - 2011  Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
            chomp;
            die "Error in line ",Nicenumber($counter),": The identifier\n
$_\n
contains white spaces\n
$hint
You could run remove_white_space_in_id.pl inputfile > newfile
This will remove everything from the id line after the first whitespace
";
        }else{
            $hash_num{$id}++;
        }
    }elsif(not (/^([A|C|G|T|N|a|c|g|t|n]+)$/ or /^\s$/)){
        chomp;
        die "Error in line ",Nicenumber($counter),": The sequence\n
$_\n
contains characters others than [acgtnACGTN]\n
$hint
";
    }
}
$counter = 1;
for(values %hash_num){
    if($_ > 1){
        $counter++;
    }
}
if($counter > 1 ){print "\nError: $counter read identifiers are not unique\n";}


exit;

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
