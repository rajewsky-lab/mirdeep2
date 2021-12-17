#!/usr/bin/env perl

# miRDeep2 perform-controls perl script
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
use File::Copy;
use File::Path;



my $usage=
"$0 file_command_line file_structure rounds_controls

-a   Output progress to screen
";

my $file_command_line=shift or die $usage;
my $file_structure=shift or die $usage;
my $rounds=shift or die $usage;

#options
my %options=();
getopts("a",\%options);


my $ltime=time();
my $dir="dir_perform_controls$ltime";


my $command_line=parse_file_command_line($file_command_line);

if($options{a}){print STDERR "total number of rounds controls=$rounds\n";}

perform_controls();


sub perform_controls{

    mkdir $dir;

    my $round=1;

    while($round<=$rounds){

	if($options{a}){print STDERR "$round\r";}

	system("permute_structure.pl $file_structure > $dir/precursors_permuted.str 2> /dev/null");

	my $ret=`$command_line 2> /dev/null`;

	print "permutation $round\n\n";

	print "$ret";

	$round++;
    }

    rmtree($dir);

    if($options{a}){print STDERR "controls performed\n\n";}
}


sub parse_file_command_line{

    my ($file) = @_;

    open (FILE, "<$file") or die "can not open $file\n";
    while (my $line=<FILE>){

	if($line=~/(\S+)/){

	    chomp $line;

	    $line=~s/$file_structure/$dir\/precursors_permuted.str/;

	    $line=~s/>.+//;

	    return $line;

	}
    }
    die "$file is empty\n";
}




