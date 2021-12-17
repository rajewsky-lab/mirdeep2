#!/usr/bin/env perl

# miRDeep2 GEO-to-FASTA perl script
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


my $usage ="

perl geo2fasta.pl GSM > outfile

reads in all files starting with GEO and produces one big output files containing only valid sequences

modified by Sebastian M

";

my $running=0;

if(not $ARGV[0]){ die $usage;}

parse_file_geo($ARGV[0]);

exit;


sub parse_file_geo{
    my @files=<@_*.txt>;

    foreach(@files){
        my $file = $_;

        print STDERR "processing file $file\n";

        open(FILE, $file) or die "Could not open file $file";

        while (my $line = <FILE>){

            if($line=~/^(\S+)\s+(\d+)/){

                my $seq=$1;

                my $cnt=$2;

                while($cnt>0){

                    if($seq =~ /^[ACGTNUacgtnu]+$/){

                        print ">deflated\_$running\n";

                        print "$seq\n";

                        $running++;

                        $cnt--;
                    }else{
                        print STDERR "$running\n$seq\n";
                        $cnt = 0;
                    }
                }
            }
        }
    }
    return;
}
