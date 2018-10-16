#!/usr/bin/perl

# miRDeep2 FASTA-parse perl script
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
"$0 fastafile

This script parses a fastafile.

-a int    only output entries where the sequence is minimum int nts long
-b        remove all entries that have a sequence that contains letters
          other than a,c,g,t,u,n,A,C,G,T,U,N.
-s        output progress
";

my $file_fasta=shift or die $usage;

my %options=();
getopts("a:bs",\%options);


my $running=0;


parse_fasta($file_fasta);


sub parse_fasta{

    my ($file_fasta) = @_;
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
		    resolve($id,$seq);
                    $id    = $1;
                    $seq   = "";
                    next;
                }
                $seq .= $_;
            }
        }
    }
    resolve($id,$seq);
    close FASTA;
    return;
}



sub resolve{

    my($id,$seq)=@_;

    $running++;

    if($options{s}){print STDERR "$running\r";}

    my $lng=length $seq;

    if($options{a} and $lng<$options{a}){
        print STDERR ">$id\n$seq\n";
	return;
    }

    if($options{b} and not($seq=~/^(a|c|g|t|u|n)+$/i)){

        print STDERR ">$id\n$seq\n";
	return;
    }

    print ">$id\n";

    print "$seq\n";

    return;
}

















