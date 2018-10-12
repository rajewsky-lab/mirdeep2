#!/usr/bin/perl

# miRDeep2 permute-structure perl script
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

use warnings;
use strict;
use Getopt::Std;


my $usage=
"$0 file_struct

Permutes the order of id's in a output file from RNAfold
";

my $file_struct=shift or die $usage;


my %hash_desc;
my %hash_seq;
my %hash_struct;
my %hash_mfe;


parse_file_struct($file_struct);
permute_order();

exit;



sub permute_order{
    
    my @ids=sort keys %hash_seq;
    my $number_ids=scalar @ids;

    my @ids_cur=@ids;
    my $number_ids_cur=$number_ids;

    for(my $i=0; $i<$number_ids; $i++){

	my $id=$ids[$i];

	my $rand=int(rand($number_ids_cur));
	my $id_cur=$ids_cur[$rand];
	
	my $seq=$hash_seq{$id_cur};
	my $struct=$hash_struct{$id_cur};
	my $mfe=$hash_mfe{$id_cur};

	splice(@ids_cur,$rand,1);
#	delete($hash_seq{$id_cur});

	$number_ids_cur--;

	print ">$id\n$seq\n$struct ($mfe)\n";
    }
}








sub parse_file_struct{
    #parses the output from RNAfoldand reads it into hashes

    my($file) = @_;
    my($id,$desc,$seq,$struct,$mfe) = ();

    open (FILE_STRUCT, "<$file") or die "can not open $file\n";
    while (<FILE_STRUCT>)
    {
        chomp;
        if (/^>(\S+)\s*(.*)/)
	{
	    $id          = $1;
	    $desc        = $2;
	    $seq         = "";
	    $struct      = "";
	    $mfe         = "";
	    while (<FILE_STRUCT>){
                chomp;
                if (/^>(\S+)\s*(.*)/){
		    $hash_desc{$id}   = $desc;
		    $hash_seq{$id}    = $seq;
		    $hash_struct{$id} = $struct;
		    $hash_mfe{$id}    = $mfe;

		    $id          = $1;
		    $desc        = $2;
		    $seq         = "";
		    $struct      = "";
		    $mfe         = "";

		    next;
                }
		if(/^\w/){
#		    tr/uU/tT/;
		    $seq .= $_;
		}if(/((\.|\(|\))+)/){
		    $struct .=$1;
		}
		if(/\((\s*-\d+\.\d+)\)/){
		    $mfe = $1;
		}
	    
	    }
        }
    }

    $hash_desc{$id}        = $desc;
    $hash_seq{$id}         = $seq;
    $hash_struct{$id}      = $struct;
    $hash_mfe{$id}         = $mfe;

    close FILE_STRUCT;
    return;
}




sub round {
    
    my($number) = shift;
    return int($number + .5);
    
}
