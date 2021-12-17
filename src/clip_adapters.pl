#!/usr/bin/env perl

# miRDeep2 clip-adapters perl script
# Copyright (C) 2008 - 2012  Marc Friedländer
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

my $usage =
"$0 file_fasta seq_adapter

Removes 3' adapters from deep sequencing reads. Works in sequence space.
";

my $file_fasta=shift or die $usage;
my $seq_adapter=shift or die $usage;

my $seq_test="TCGTATGCCGTCTTCTGCTTGT";

my $prefix=substr($seq_adapter,0,6);

$prefix=~tr/[acgtun\.]/[ACGTTNN]/;

remove_adapters($file_fasta,$prefix);

exit;


sub remove_adapters{

    my ($file_fasta,$prefix) = @_;
    my ($id,$seq) = ();

    open (FASTA, "<$file_fasta") or die "can not open $file_fasta\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)/)
        {
            $id  = $1;
	    $seq = "";
            while (<FASTA>){
                chomp;
                if (/^>(\S+)/){
		    remove_adapter($id,$seq,$prefix);
                    $id    = $1;
                    $seq   = "";
                    next;
                }
                $seq .= $_;
            }
        }
    }
    remove_adapter($id,$seq,$prefix);
    close FASTA;
    return;
}





sub remove_adapter{

    my($id,$seq,$prefix)=@_;

    $seq=~tr/[acgtun\.]/[ACGTTNN]/;

    my $seq_clipped;

    if($seq=~/(\w+)$prefix/){

	$seq_clipped=$1;

    }elsif(substr($seq,0,6) eq $prefix){
	$seq_clipped=$prefix;

    }else{

	my $finish=0;

	while(not $finish and (length($prefix)>0)){

	    chop $prefix;

	    if($seq=~/(\w+)$prefix$/){

		$seq_clipped=$1;
		$finish=1;
	    }
	}

    }

    if(not $seq_clipped){


	$seq_clipped=$seq;
    }

    print ">$id\n$seq_clipped\n";
}
