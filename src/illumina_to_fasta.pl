#!/usr/bin/perl

# miRDeep2 Illumina-to-FASTA perl script
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
"$0 _seq.txt

This script parses _seq.txt to fasta format. Options:
-a    format is qseq.txt

Example of use:
$0 mouse_seq.txt > mouse.fa
";

my $file_solexa=shift or die $usage;

my %options=();
getopts("a",\%options);

my $running=0;



parse_file_solexa($file_solexa);

exit;




sub parse_file_solexa{

    my ($file) = @_;

    open (FILE, "<$file") or die "can not open $file\n";
    while (<FILE>){
		next if(/^\s*$/);
	my @fields=split(/\s+/);

	my $seq;

	if($options{a}){

	    $seq=$fields[8];

	}else{	

	    $seq=$fields[4];
	}

	chomp $seq;

	$seq=~s/\./N/g;

	print ">seq\_$running\n";
	
	print "$seq\n";

	$running++;
    }
    return;
}

