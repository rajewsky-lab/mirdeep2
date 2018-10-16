#!/usr/bin/perl -W

# miRDeep2 sanity-check-mapping-file perl script
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
my $rid;
my $rl;
my $rs;
my $re;
my $rseq;
my $gid;
my $gl;
my $gs;
my $ge;
my $gseq;
my $strand;
my $mm;
my $edit;

my $counter = 0;

while(<>){
    chomp;
    $counter++;
    #    read-id   l       s       e      seq   gid     l        s       e      seq
    if(/^(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+([+-])\s+(\d+)\s*([mDIM]*)$/){
     
    }else{
        die "\nWrong format in line ",Nicenumber($counter),": The row\n
$_\n
does not correspond to the format\n
read_id_wo_whitespaces\tlength\tstart\tend\tread_sequence         \tgenomicID_wo_whitspaces\tlength\tstart\tend\tgenomic_sequence       \tstrand\t#mismatches\teditstring\n
e.g. read_22_x10000 \t22    \t1    \t22 \tagtcgtgactgactgactgacg\tchromosomeIII_x12312312\t22    \t1001 \t1022\tagtcgtgactgactgactgacg\t+-    \t0          \tmmmmmmmmmmmmmmmmmmmmm\n
Please make sure that all lines have the above described format.
\n
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
