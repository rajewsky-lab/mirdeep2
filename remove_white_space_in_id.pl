#!/usr/bin/perl

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
