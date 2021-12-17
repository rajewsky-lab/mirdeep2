#!/usr/bin/env perl

# miRDeep2 find-read-count perl script
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

use warnings;
use strict;
use Getopt::Std;

my $usage =
"$0 file_fasta

finds the total read count if '_x\\d' format is used.
";

my $file_fasta=shift or die $usage;

my %options=();
getopts("",\%options);


my $count_reads=0;

open (FILE, "<$file_fasta") or die "can not open $file_fasta\n";

while (my $line=<FILE>){

    if($line=~/_x(\d+)/){
	$count_reads+=$1;
    }
}

print "$count_reads\n";
