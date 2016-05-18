#!/bin/bash

bowtie-build cel_cluster.fa cel_cluster
ec=`echo $?`
if [ $ec != 0 ];then
	echo An error occured, exit code $ec
fi

mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p cel_cluster -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v -n 
ec=`echo $?`
if [ $ec != 0 ];then
	echo An error occured, exit code $ec
fi

miRDeep2.pl reads_collapsed.fa cel_cluster.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa -t C.elegans

ec=`echo $?`
if [ $ec != 0 ];then
	echo An error occured, exit code $ec
fi

