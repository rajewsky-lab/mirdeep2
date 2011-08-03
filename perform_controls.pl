#!/usr/bin/perl

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




