#!/usr/bin/perl

# miRDeep2 FASTQ-to-FASTA perl script
# Copyright (C) 2008 - 2011  Marc Friedländer
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

my $file = $ARGV[0] or die "Usage: \n\n $0 fastq file > outputfile\n";;

open IN,"<$file" or die "File $file not found error here\n";
my $c = 0;
my @line = ();
my $id ='';
my $seq = '';

my $processed = 0;

while(<IN>){
    chomp;
    $c++;
    if($c == 1){
        $processed++;
#        print STDERR "$processed reads processed\r";
        @line = split();
        $id = $line[0];
        $id =~ s/\@//;
        
    }elsif($c == 2){
        $seq = $_;
    }elsif($c == 4){
        $c = 0;
        print ">seq_$processed\n$seq\n";
        $id ='';
        @line =();
        $seq ='';
    }else{}
}
close IN;
exit;
#print STDERR "$processed reads processed\n\n";
