package Bio::miRDeep2;

use 5.020002;
use strict;
use warnings;
use Data::Dumper;
use Carp qw(verbose confess carp croak);
require Exporter;
our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use pipeline-util ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
#our %EXPORT_TAGS = ( 'all' => [ qw(
#	'TargetVcfReAnnotator'
#) ] );

#our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(
	LinePosLastHandleLogger
	TargetVcfReAnnotator
	WalkToTarget
	NewWalk
	WalkToNext
	AnnotateTargetRecords
	ADFilterTargetRecords
	FormatWalkTargetLineAsVcfLine
	FormatWalkTargetLineAsVcfHeader
	_formatwalkasvcflineswithfile
);

our $VERSION = "0.1.3";

# Preloaded methods go here.
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Bio::miRDeep2 - miRDeep2 cmdline tools for analysis of small rna hts sequencing datasets 

=head1 SYNOPSIS

miRDeep2 is a software package for identification of novel and known miRNAs in
 deep sequencing data. Furthermore, it can be used for miRNA expression
 profiling across samples. Last, a new script for preprocessing of raw Illumina
 sequencing data produces files for downstream analysis with the miRDeep2 or
 quantifier script. Colorspace sequencing data is currently not supported by the
 preprocessing module but it is planed to be implemented. Preprocessing is
 performed with the mapper.pl script. Quantification and expression profiling is
 done by the quantifier.pl script. miRNA identification is done by the miRDeep2.pl script.

=head1 DESCRIPTION

read README.md on github

=head2 EXPORT

None by default.

=head1 SEE ALSO

	https://github.com/rajewsky-lab/mirdeep2

=head1 AUTHOR

Sebastian Mackowiak & Marc Friedländer


=head1 COPYRIGHT AND LICENSE

Copyright Sebastian Mackowiak & Marc Friedländer

GPLv3 licence

=cut

