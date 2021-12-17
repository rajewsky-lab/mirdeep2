#!/usr/bin/env perl

# miRDeep2 FASTA-select perl script
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
"$0 file_fasta idfile (-a)

This script only prints out the fastaentries that has an id
that is present in the idfile. Using the option -a only prints
out entries that has an id that is not present in the id file.
";

my $file_fasta=shift or die $usage;
my $file_ids=shift or die $usage;

my %options=();
getopts("a",\%options);

my %hash;



parse_file_ids(\$file_ids,\%hash);

parse_fasta(\$file_fasta,\%hash);

exit;




sub parse_fasta{

    my ($file,$hash) = @_;
    my ($id, $desc, $sequence) = ();

    open (FASTA, "<$$file") or die "can not open $$file\n";
    while (<FASTA>)
    {
        chomp;
        if (/^>(\S+)(.*)/)
        {
            $id       = $1;
            $desc     = $2;
            $sequence = "";
            while (<FASTA>){
                chomp;
                if (/^>(\S+)(.*)/){

		    if((defined $$hash{$id} and not $options{a}) or (not defined $$hash{$id} and $options{a})){

			print ">$id$desc\n$sequence\n";
		    }

                    $id         = $1;
                    $desc       = $2;
                    $sequence   = "";
                    next;
                }
                $sequence .= $_;
            }
        }
    }
    if((defined $$hash{$id} and not $options{a}) or (not defined $$hash{$id} and $options{a})){

    	print ">$id$desc\n$sequence\n";
    }

    close FASTA;
    return;
}






sub parse_file_ids{

    my ($file,$hash) = @_;

    open (FILE, "<$$file") or die "can not open $$file\n";

    while (my $line=<FILE>){

	if($line=~/^(\S+)/){

	    my $id=$1;

	    $$hash{$id}=1;
	}
    }
    return;
}
