#!/usr/bin/perl

# miRDeep2 extract-miRNAs perl script
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

my $id;
my $seq;
my @species = split(",",$ARGV[1]);
die "Error: Specify species or comma separated list of species\n\nUsage:\n\tperl $0 miRBase_fasta_file species [mature/star]\n\te.g. $0 mature.fa cel mature\n\n" if(not $ARGV[1]);
my $in = 0;

open IN,"<$ARGV[0]" or die "$ARGV[0] not found $!";

my $first =1;

my $str;

my $sp;

my $type='mature';
if($ARGV[2]){
	if(lc($ARGV[2]) eq 'mature' or lc($ARGV[2]) eq 'star'){
		$type = $ARGV[2];
	}else{
		die 'argv 2 must either be mature or star\n';
	}
}

foreach my $sp(@species){

$in = 0;

open IN,"<$ARGV[0]" or die "$ARGV[0] not found $!";

$first =1;

$str;

my $id;

#die $type;

while(<IN>){
    chomp;
    if(/(>$sp\S+)/){
		$id=$1;
		if($type eq 'mature' and $id =~/\*$/){ $in= 0; next;}
		if($type eq 'star' and $id !~ /\*/){$in=0; next;}
	#print "=== $id\n";
       if($first){
            print "$id\n";
            $first = 0;
        }else{
            print "\n$id\n";
        }
        $in =1;
        next;
     }elsif(/>/){
         $in = 0;
     }elsif($in == 1){
	#print "\n==== $_ ===\n";
       $str = uc($_);
       $str =~ tr/U/T/;
       print $str;
     }else{
	$id='';
	$str='';
	$in=0;
     }
}

if(not $first){
    print "\n";
}
close IN;
}
