#!/usr/bin/perl 

# miRDeep2 mirdeep-to-BED perl script
# Copyright (C) 2009 - 2012  Sebastian Mackowiak
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

open IN,"$ARGV[0]" or die "No csv file given for bed conversion\n";

my ($known,$novel, $not, $line,$thres,$score,$line,@l,$strand,$label,$end);

while(<IN>){
	if(/novel miRNAs predicted by miRDeep2/){
		$novel=1;
		$known=0;
		$not=0;
	}elsif(/mature miRBase miRNAs detected/){
		$novel=0;
		$known=1;
		$not=0;
	}elsif(/miRBase miRNAs not detected/){
		last;
	}else{
		@l=split();
		if($l[$#l] =~ /(\S+):(\d+)\.\.(\d+):(\S)/){
			$end=$3;
			$strand='255,0,0' if($4 eq '+');
			$strand='0,0,255' if($4 eq '-');
			$label="known" if($known);
			$label="novel" if($novel);
			print "$1\t$2\t$end\t${label}:$l[0]\t$l[1]\t$4\t$2\t$end\t$strand\n";
		}
	}
}
