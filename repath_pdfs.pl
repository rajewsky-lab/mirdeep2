#!/usr/bin/perl

use strict;
use File::Basename;
use Cwd;

my $cwd = cwd;
$cwd .="/";

my $pdfsdir = $ARGV[1] or die "Usage: $0 html_file pdfs_dir > outfile.html\n";

open IN,"$ARGV[0]" or die "No html input file given\n";
my ( $name0, $path0, $extension0 ) = fileparse ( $ARGV[0], '\..*' );


open OUT,">${name0}_repath.html" or die "Could not create file ${name0}_repath.html\n";

while(<IN>){
	if(/a href="file:\/\/(\/\S+)$pdfsdir/){
		$_ =~ s/$1/$cwd/;
	}
	print OUT;
}
close IN;
close OUT;
