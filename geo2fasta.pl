#!/usr/bin/perl

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
