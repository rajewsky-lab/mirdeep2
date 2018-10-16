#!/usr/bin/perl

# miRDeep2 remove-white-space-in-ID perl script
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

my $usage =
"$0 file_fasta

Removes whitespaces from id line in a fasta file\n
";

my $file_fasta=shift or die $usage;


parse_fasta($file_fasta);

exit;




sub parse_fasta{

    my ($file) = @_;

    open (FASTA, "<$file") or die "can not open $file\n";
    while (<FASTA>){

		if (/^(>\S+)/){

			print "$1\n";

		}else{
			print uc($_);
		}
	}

    close FASTA;
    return 0;
}
