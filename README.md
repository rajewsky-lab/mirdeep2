[![Build Status](https://travis-ci.org/rajewsky-lab/mirdeep2.svg?branch=master)](https://travis-ci.org/rajewsky-lab/mirdeep2)

# miRDeep2 `README`

## About

Authors: Sebastian Mackowiak & Marc Friedländer

This is miRDeep2 developed by Sebastian Mackowiak & Marc Friedländer.
miRDeep2 discovers active known or novel miRNAs from deep sequencing data
(Solexa/Illumina, 454, ...).

(minor edits to `README`, `TUTORIAL`, `CHANGELOG`, and `FAQ`, convertion to
Markdown, trailing whitespace removal & CI setup by Marcel Schilling)


## Requirements

Linux system, 2GB Ram, enough disk space dependent on your deep sequencing data


## Testing version

MacOSX with Xcode and gcc compiler installed. (This can be obtained from the
appstore, if there are any issues with installing it please look for help
online).

To compile the Vienna package it may be necessary to have GNU grep installed
since the MacOSX grep is BSD based and sometimes not accepted by the installer.
To get a GNU grep you could for example install homebrew by typing

```sh
ruby -e "$(curl -fsSL \
  https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

(the link could be out of date, in that case look up online what to do)

After that typing

```sh
brew tap homebrew/dupes; brew install grep
```

will install GNU grep as `ggrep` in `/usr/local/bin/`


## Installation

### Option 1: with the provided `install.pl` script

Type

```sh
perl install.pl
```

### Option 2. without the install mirdeep script

Follow the instructions given below

#### Dependencies

First download all necessary packages listed here

1. [bowtie short read aligner][bowtie-source]
2. [Vienna package with RNAfold][vienna-source]
3. [SQUID library][squid-source] goto Squid and download it
4. [randfold][randfold-source]
5. [Perl package PDF::API2][pdf-api2-source]

[bowtie-source]: http://bowtie-bio.sourceforge.net/index.shtml
[vienna-source]: http://www.tbi.univie.ac.at/~ivo/RNA/
[squid-source]: http://eddylab.org/software.html
[randfold-source]: http://bioinformatics.psb.ugent.be/software/details/Randfold
[pdf-api2-source]: http://search.cpan.org/search?query=PDF%3A%3AAPI2&mode=all

#### Manual installation

When packages are downloaded

1. attach the miRDeep2 executable path to your PATH

```sh
echo 'export PATH=$PATH:your_path_to_mirdeep2/src' >> ~/.bashrc
```

2. `unzip bowtie-0.11.3-bin-linux-x86_64.zip`

3. put the bowtie directory into your `PATH` variable, *e.g.*

```sh
echo 'export PATH=$PATH:your_path_tobowtie' >> ~/.bashrc
```

4. `tar xvvzf ViennaRNA-1.8.4.tar.gz`

5. `cd` to the Vienna dir

6. type

```sh
./configure --prefix=your_path_to_Vienna/install_dir
make
make install
```

7. add Vienna binaries to your `PATH` variable, *e.g.*

```sh
echo 'export PATH=$PATH:your_path_to_Vienna/install_dir/bin' >> ~/.bashrc
```

8. `tar xxvzf squid-1.9g.tar.gz`

9. `tar xvvzf randfold-2.0.tar.gz`

10. `cd randfold2.0`

11. edit Makefile, *e.g.* `emacs Makefile`:

change line with `INCLUDE=-I.` to
`INCLUDE=-I. -I<your_path_to_squid-1.9g> -L<your_path_to_squid-1.9g>`,
*e.g.* `INCLUDE=-I. -I/home/Pattern/squid-1.9g/ -L/home/Pattern/squid-1.9g/`

12. `make`

13. add randfold to your `PATH` variable, *e.g.*

```sh
echo 'export PATH=$PATH:your_path_to_randfold' >> ~/.bashrc
```

14. `tar xvvzf PDF-API2-0.73.tar.gz`

15. `cd` to your PDF_API2 directory

16. then type in

```sh
perl Makefile.PL INSTALL_BASE=your_path_to_miRDeep2 LIB=your_path_to_miRDeep2/lib
make
make test
make install
```

17. add your library to the `PERL5LIB`, *e.g.*

```sh
echo \
  'export PERL5LIB=PERL5LIB:your_path_to_miRDeep2/lib/perl5' \
  >> ~/.bashrc
```


18. `cd` to your mirdeep2 directory (the one containing `install.pl`)

19. `touch install_successful`

20. start a new shell session to apply the changes to environment variables

#### Test installation

To test if everything is installed properly type in

1. `bowtie`
2. `RNAfold -h`
3. `randfold`
4. `make_html.pl`

You should not get any error messages. Otherwise something is not correctly
installed.


### Install Paths

Everything that is download by the installer will be in a directory called
`<your_path_to_mirdeep2>/essentials`


## Script Reference

miRDeep2 analyses can be performed using the three scripts `miRDeep2.pl`,
`mapper.pl` and `quantifier.pl`.


### `miRDeep2.pl`

#### Description

Wrapper function for the `miRDeep2.pl` program package. The script runs all
necessary scripts of the miRDeep2 package to perform a microRNA detection deep
sequencing data anlysis.

#### Input

* A FASTA file with deep sequencing reads,
* a FASTA file of the corresponding genome,
* a file of mapped reads to the genome in miRDeep2 ARF format,
* an optional FASTA file with known miRNAs of the analysed species, and
* an optional FASTA file of known miRNAs of related species.

#### Output

* A spreadsheet and
* an HTML file

with an overview of all detected miRNAs in the deep sequencing input data.

#### Options

| option         | description                                                |
|----------------|------------------------------------------------------------|
| `‑a <int>`     | minimum read stack height that triggers analysis. Using this option disables automatic estimation of the optimal value. |
| `‑b <int>`     | minimum score cut-off for predicted novel miRNAs to be displayed in the overview table. This score cut-off is by default 0. |
| `‑c`           | disable randfold analysis                                  |
| `‑t <species>` | species being analyzed - this is used to link to the appropriate UCSC browser |
| `‑u`           | output list of UCSC browser species that are supported and exit |
| `‑v`           | remove directory with temporary files                      |
| `‑q <file>`    | `miRBase.mrd` file from quantifier module to show miRBase miRNAs in data that were not scored by miRDeep2 |

#### Examples:

The miRDeep2 module identifies known and novel miRNAs in deep sequencing data.
The output of the mapper module can be directly plugged into the miRDeep2
module.

##### Example use 1

The user wishes to identify miRNAs in mouse deep sequencing data, using default
options.
The `miRBase_mmu_v14.fa` file contains all miRBase mature mouse miRNAs, while
the `miRBase_rno_v14.fa` file contains all the miRBase mature rat miRNAs.
The `2>` will pipe all progress output to the `report.log` file.

```sh
miRDeep2.pl reads_collapsed.fa genome.fa reads_collapsed_vs_genome.arf \
  miRBase_mmu_v14.fa miRBase_rno_v14.fa precursors_ref_this_species.fa \
  -t Mouse 2>report.log
```

This command will generate

* a directory with PDFs showing the structures, read signatures and score
  breakdowns of novel and known miRNAs in the data,
* an HTML webpage that links to all results generated (`result.html`),
* a copy of the novel and known miRNAs contained in the webpage but in text
  format which allows easy parsing (`result.csv`),
* a copy of the performance survey contained in the webpage but in text format
  (`survey.csv`), and
* a copy of the miRNA read signatures contained in the PDFs but in text format
  (`output.mrd`).

##### Example use 2

The user wishes to identify miRNAs in deep sequencing data from an animal with
no related species in miRBase:

```sh
miRDeep2.pl reads_collapsed.fa genome.fa reads_collapsed_vs_genome.arf \
  none none none 2>report.log
```

This command will generate the same type of files as example use 1 above.
Note that there it will in practice always improve miRDeep2 performance if
miRNAs from some related species is input, even if it is not closely related.

---


### `mapper.pl`

#### Description

Processes reads and/or maps them to the reference genome.

#### Input

Default input is

* a file in FASTA, `seq.txt` or `qseq.txt` format.

More input can be given depending on the options used.

#### Output

The output depends on the options used (see below).

Either

* a FASTA file with processed reads, or
* an ARF file with with mapped reads, or
* both

are output.

#### Options

##### Read input file

| option | description                     |
|--------|---------------------------------|
| `‑a`   | input file is `seq.txt` format  |
| `‑b`   | input file is `qseq.txt` format |
| `‑c`   | input file is FASTA format      |

##### Preprocessing/mapping

| option        | description                                                 |
|---------------|-------------------------------------------------------------|
| `‑h`          | parse to FASTA format                                       |
| `‑i`          | convert RNA to DNA alphabet (to map against genome)         |
| `‑j`          | remove all entries that have a sequence that contains letters other than `a`, `c`, `g`, `t`, u, `n`, `A`, `C`, `G`, `T`, `U`, or `N`. |
| `‑k <seq>`    | clip 3' adapter sequence                                    |
| `‑l <int>`    | discard reads shorter than `<int>` nts                      |
| `‑m`          | collapse reads                                              |
| `‑p <genome>` | map to genome (must be indexed by `bowtie-build`). The `genome` string must be the prefix of the bowtie index. For instance, if the first indexed file is called `h_sapiens_37_asm.1.ebwt` then the prefix is `h_sapiens_37_asm`. |
| `‑q`          | map with one mismatch in the seed (mapping takes longer)    |

##### Output files

| option    | description                        |
|-----------|------------------------------------|
| `‑s file` | print processed reads to this file |
| `‑t file` | print read mappings to this file   |

##### Other

| option | description                                  |
|--------|----------------------------------------------|
| `‑u`   | do not remove directory with temporary files |
| `‑v`   | outputs progress report                      |

#### Examples

The mapper module is designed as a tool to process deep sequencing reads and/or
map them to the reference genome. The module works in sequence space, and can
process or map data that is in sequence FASTA format.
A number of the functions of the mapper module are implemented specifically
with Solexa/Illumina data in mind. For example on how to post-process mappings
in color space, see example use 5:

##### Example use 1

The user wishes to parse a file in `qseq.txt` format to FASTA format, convert
from RNA to DNA alphabet, remove entries with non-canonical letters (letters
other than `a`, `c`, `g`, `t`, `u`, `n`, `A`, `C`, `G`, `T`, `U`, or `N`), clip
adapters, discard reads shorter than 18 nts and collapse the reads:

 ```sh
mapper.pl reads_qseq.txt -b -h -i -j -k TCGTATGCCGTCTTCTGCTTGT -l 18 -m \
  -s reads_collapsed.fa
```

##### Example use 2

The user wishes to map a FASTA file against the reference genome.
The genome has already been indexed by `bowtie-build`.
The first of the indexed files is named `genome.1.ebwt`:

```sh
mapper.pl reads_collapsed.fa -c -p genome -t reads_collapsed_vs_genome.arf
```

##### Example use 3

The user wishes to process the reads as in example use 1 and map the reads as
in example use 2 in a single step, while observing the progress:

```sh
mapper.pl reads_qseq.txt -b -h -i -j -k TCGTATGCCGTCTTCTGCTTGT -l 18 -m \
  -p genome -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v
```

##### Example use 4

The user wishes to parse a GEO file to FASTA format and process it as in
example use 1.
The GEO file is in tabular format, with the first column showing the sequence
and the second column showing the read counts:

```sh
geo2fasta.pl GSM.txt > reads.fa

mapper.pl reads.fa -c -h -i -j -k TCGTATGCCGTCTTCTGCTTGT -l 18 -m \
  -s reads_collapsed.fa
```

##### Example use 5

The user has already removed 3' adapters in color space and has mapped the
reads against the genome using the BWA tool. The BWA output file is named
`reads_vs_genome.sam`. Notice that the BWA output contains extra fields that
are not required for SAM format. Our converter requires these fields and thus
may not work with all types of SAM files. The user wishes to generate
`reads_collapsed.fa` and `reads_vs_genome.arf` to input to miRDeep2:

```sh
bwa_sam_converter.pl reads_vs_genome.sam reads.fa reads_vs_genome.arf

mapper.pl reads.fa -c -i -j -l 18 -m -s reads_collapsed.fa
```

---


### `quantifier.pl`

#### Description

The module maps the deep sequencing reads to predefined miRNA precursors and
determines by that the expression of the corresponding miRNAs.
First, the predefined mature miRNA sequences are mapped to the predefined
precursors. Optionally, predefined star sequences can be mapped to the
precursors too. By that the mature and star sequence in the precursors are
determined.
Second, the deep sequencing reads are mapped to the precursors. The number of
reads falling into an interval 2 nt upstream and 5 nt downstream of the
mature/star sequence is determined.

#### Input

* A FASTA file with precursor sequences,
* a FASTA file with mature miRNA sequences,
* a FASTA file with deep sequencing reads, and
* optionally a FASTA file with star sequences and the 3 letter code of the
  species of interest.

#### Output

* A 2 column table file called `miRNA_expressed.csv` with miRNA identifiers and
  its read count,
* a file called `miRNA_not_expressed.csv` with all miRNAs having 0 read counts,
* a signature file called `miRBase.mrd`,
* a file called `expression.html` that gives an overview of all miRNAs the
  input data, and
* a directory called `pdfs` that contains for each miRNA a PDF file showing its
  signature and structure.

#### Options

| option       | description                                                                                  |
|--------------|----------------------------------------------------------------------------------------------|
| -p [file.fa] | miRNA precursor sequences (around 70bp: One line per precursors sequence)
| -m [file.fa] | mature miRNA sequences (around 22nt)
| -P           | specify this option of your mature miRNA file contains 5p and 3p ids only
|	-c [file]    | config.txt file with different sample ids... or just the one sample id  -- deprecated
|	-s [star.fa] | optional star sequences from miRBase
|	-t [species] | e.g. Mouse or mmu
|	             | if not searching in a specific species all species in your files will be analyzed
|	             | else only the species in your dataset is considered
|	-y [time]    | optional otherwise its generating a new one
|	-d           | if parameter given pdfs will not be generated, otherwise pdfs will be generated
|	-o           | if parameter is given reads were not sorted by sample in pdf file, default is sorting
|	-k           | also considers precursor-mature mappings that have different ids, eg let7c
|	             | would be allowed to map to pre-let7a
|	-n           | do not do file conversion again
|	-x           | do not do mapping against precursor again
|	-g [int]     | number of allowed mismatches when mapping reads to precursors, default 1
|	-e [int]     | number of nucleotides upstream of the mature sequence to consider, default 2
|	-f [int]     | number of nucleotides downstream of the mature sequence to consider, default 5
|	-j           | do not create an output.mrd file and pdfs if specified
|	-W           | read counts are weighed by their number of mappings. e.g. A read maps twice so each position 
|              | gets 0.5 added to its read profile
|	-U           | use only unique read mappings; Caveat: Some miRNAs have multiple precursors. These will be 
|              | underestimated in their expression since the multimappers are excluded
| -u	         | list all values allowed for the species parameter that have an entry at UCSC

#### Example usage

```sh
quantifier.pl -p precursors.fa -m mature.fa -r reads.fa
```

---


### `make_html.pl`

#### Description

It creates a file called `result.html` that gives an overview of miRDeep2
detected miRNAs (known and novel ones). The HTML file lists up each detected
miRNA and provides among others information on its miRDeep2 score, reads mapped
to its mature, loop and star sequence, the mature, star and consensus precursor
sequences themselves and provides links to BLAST, BLAT, mirBase for miRBase
miRNAs and to a PDF file that shows the signature and structure.

#### Input

* A miRDeep2 output.mrd file and
* a miRDeep2 survey.csv file

#### Output

* A `result.html` file with an entry for each provisional miRNA that contains
  information about its assigned Id, miRDeep2 score, estimated probability that
  the miRNA candidate is a true positive, rfam alert, total read count, mature
  read count, loop read count, star read count, significant randfold p-value,
  miRBase miRNA, example miRBase miRNA with the same seed, BLAT, BLAST,
  consensus mature sequence, consensus star sequence and consensus precursor
  sequence. Furthermore, the miRBase miRNAs existent in the input data but not
  scored by miRDeep2 are listed.
* A directory called `pdfs` that contains for each provisional miRNA ID a PDF
  with its signature and structure.
* A file called `result.csv` (when option `-c` is used) that contains the same
  entries as the HTML file.

#### Options

| option                | description                                         |
|-----------------------|-----------------------------------------------------|
| `‑v <int>`            | only output hairpins with score above `<int>`       |
| `‑c`                  | also create overview in excel format                |
| `‑k <file>`           | supply file with known miRNAs                       |
| `‑s <file>`           | supply survey file if score cutoff is used to get information about how big is the confidence of resulting reads |
| `‑f <file>`           | miRDeep2 output MRD file                            |
| `‑e`                  | report complete survey file                         |
| `‑g`                  | report survey for current score cutoff              |
| `‑w <project_folder>` | automatically used when running webinterface, otherwise don't use it |
| `‑r <file>`           | Rfam file to check for already reported small RNA sequences |
| `‑q <file>`           | `miRBase.mrd` file produced by quantifier module    |
| `‑x <file>`           | `signature.arf` file with mapped reads to precursors |
| `‑t <org>`            | specify the organism from which your sequencing data was obtained |
| `‑u`                  | print all available UCSC input organisms            |
| `‑d`                  | do not generate PDFs                                |
| `‑y`                  | timestamp                                           |
| `‑z`                  | switch is automatically used when script is called by `quantifier.pl` |
| `‑o`                  | print reads in PDF signature sorted by their 3 letter code in front of their identifier |

#### Example usage

```sh
make_html.pl -f miRDeep_outfile -s survey.csv -c -e -y 123456789
```

---


### `clip_adapters.pl`

#### Description

Removes 3' end adaptors from deep sequenced small RNAs. The script searches for
occurrences of the six first nucleotides of the adapter in the read sequence,
starting after position 18 in the read sequence (so the shortest clipped read
will be 18 nts). If no matches to the first six nts of the adapter are
identified in a read, the 3' end of the read is searched for shorter matches to
the 5 to 1 first nts of the adapter.

#### Input

* A FASTA file with the deep sequencing reads and the adapter sequence (both in
  RNA or DNA alphabet).

#### Output

* A FASTA file with the clipped reads.

FASTA IDs are retained. If no matches to the adapter prefixes are identified in
a given read, the unclipped read is output.

#### Example usage

```sh
clip_adapters.pl reads.fa TCGTATGCCGTCTTCTGCTTGT > reads_clipped.fa
```

#### Notes

It is possible to clip adapters using more sophisticated methods.
Users are encouraged to test other methods with the miRDeep2 modules.

---


### `collapse_reads.pl`

#### Description

Collapses reads in the FASTA file to ensure that each sequence only occurs
once.
To indicate how many times reads the sequence represents, a suffix is added to
each FASTA identifier. *E.g.* a sequence that represents ten reads in the data
will have the `_x10` suffix added to the identifier.

#### Input

* A FASTA file, either in standard format or in the collapsed suffix format.

#### Output

* A FASTA file in the collapsed suffix format.

#### Options

| option | description      |
|--------|------------------|
| `‑a`   | outputs progress |

#### Example usage

```sh
collapse_reads.pl reads.fa > reads_collapsed
```

#### Notes

Since the script reads all FASTA entries into a hash using the sequence as key,
it can potentially use more than 3 GB memory when collapsing very big datasets,
\>50 million reads. In this case, the user can partition the reads
(for instance based on the 5' nucleotide), collapse separately and concatenate.

---


### `excise_precursors_iterative.pl`

#### Description

This script is a wrapper for `excise_precursors.pl`, which it calls one or more
times, incrementing the height of the read stack required for initiating
excision until the number of excised precursors falls below a given threshold.

#### Input

* The reference genome in FASTA format,
* the mapped reads in `.arf` format,
* a filename that the excised precursors will be written to, and
* the maximal number of precursors that should be reported.

#### Output

## The excised precursors in FASTA format.

#### Options

| option | description                |
|--------|----------------------------|
| `‑a`   | Output progress to screen. |

#### Example usage

```sh
excise_precursors_iterative.pl genome.fa reads_vs_genome.arf \
  potential_precursors.fa 50000 -a
```

---


### `excise_precursors.pl`

#### Description

Excises precursors from the genome using the mapped reads as guidelines.

#### Input

* The reference genome in FASTA format and
* the mapped reads in `.arf` format.

#### Output

* The excised precursors in FASTA format.

## Options

| option         | description                                                |
|----------------|------------------------------------------------------------|
| `‑a <integer>` | Only excise if the highest local read stack is `<integer>` reads high (default 2). |
| `‑b`           | Output progress to screen.                                 |

## Example usage

```sh
excise_precursors.pl genome.arf reads_vs_genome.arf -b
```

---


### `fastaparse.pl`

#### Description

Performs simple filtering of entries in a FASTA file.

#### Input

* A FASTA file.

#### Ouput

* A filtered FASTA file.

#### Options

| option     | description                                                    |
|------------|----------------------------------------------------------------|
| `‑a <int>` | only output entries where the sequence is minimum int nts long |
| `‑b`       | remove all entries that have a sequence that contains letters other than `a`, `c`, `g`, `t`, `u`, `n`, `A`, `C`, `G`, `T`, `U`, or `N`. |
| `‑s`       | output progress                                                |

#### Example usage

```sh
fastaparse.pl reads.fa -a 18 -s > reads_no_short.fa
```

---


### `fastaselect.pl`

#### Description

This script only prints out the FASTA entries that match an ID in the ID file.

#### Input

* A FASTA file and a file with IDs, one ID per line.

#### Output

* A FASTA file containing the FASTA entries that match an ID.

#### Options

| option | description                                                        |
|--------|--------------------------------------------------------------------|
| `‑a`   | only prints out entries that has an id that is not present in the ID file. |

#### Example usage

```sh
fastaselect.pl reads.fa reads_select.ids > reads_select.fa
```

---


### `find_read_count.pl`

#### Description

Scans a file searching for the suffixes that are generated by
`collapse_reads.pl` (e.g. `_x10`).
It sums up the integer values in the suffixes and outputs the sum. If a given
id occurs multiple times in the file, it will multi-count the integer value of
the ID. It will also only count the first integer occurrence in a given line.

#### Input

* Any file containing the suffixes that are generated by `collapse_reads.pl`.

This will typically be a FASTA file or a list of IDs.

#### Output

* The sum of integer values (the total read count).


#### Example usage

```sh
find_read_count.pl reads_collapsed.fa
```

---


### `geo2fasta.pl`

#### Description

Parses GSM format files into FASTA format.

#### Input

* GSM files in tabular format.

The first column should be sequences and the second column the number of times
the sequence occurs in the data.

#### Output

* A FASTA file, one sequence per line (the sequences are expanded).


#### Example usage

```sh
geo2fasta.pl GSM.txt > reads.fa
```

---


### `illumina_to_fasta.pl`

#### Description

Parses `seq.txt` or `qseq.txt` output from the Solexa/Illumina platform to
FASTA format.

#### Input

* A `seq.txt` or
* `qseq.txt` file.

By default `seq.txt`.

#### Output

* A FASTA file, one entry for each line of `seq.txt`.

The entries are named `seq` plus a running number that is incremented by one
for each entry. Any `.` characters in the `seq.txt` file is substituted with an
`N`.

#### Options

| option | description          |
|--------|----------------------|
| `‑a`   | format is `qseq.txt` |

#### Example usage

```sh
illumina_to_fasta.pl s_1.qseq.txt -a > reads.fa
```

---


### `miRDeep2_core_algorithm.pl`

#### Description

For each potential miRNA precursor input, the miRDeep2 core algorithm either
discards it or assigns it a log-odds score that reflects the probability that
the precursor is a genuine miRNA.

#### Input

Default input is

* an ARF file with the read signatures and
* an RNAfold output file with the structures of the potential miRNA precursors.

#### Output

* A .mrd file with all potential miRNA precursors that are scored.

#### Options

| option | description                                                        |
|--------|--------------------------------------------------------------------|
| `‑h`   | print this usage                                                   |
| `‑s`   | FASTA file with reference mature miRNAs from one or more related species |
| `‑t`   | print filtered                                                     |
| `‑u`   | limited output (only ids)                                          |
| `‑v`   | cut-off (default 1)                                                |
| `‑x`   | sensitive option for Sanger sequences                              |
| `‑y`   | file with randfold p-values                                        |
| `‑z`   | consider Drosha processing                                         |

#### Example usage

```sh
miRDeep2_core_algorithm.pl signature.arf potential_precursors.str \
  -s miRBase_related_species.fa -y potential_precursors.rand > output.mrd
```

#### Notes

The `-z` option has not been thoroughly tested.

---


### `parse_mappings.pl`

#### Description

Performs simple filtering of entries in an `.arf` file.

#### Input

Default input is

* an `.arf` file.

#### Output

* A filtered `.arf` file.

#### Options

| option      | description                                                   |
|-------------|---------------------------------------------------------------|
| `‑a <int>`  | Discard mappings of edit distance higher than this            |
| `‑b <int>`  | Discard mappings of read queries shorter than this            |
| `‑c <int>`  | Discard mappings of read queries longer than this             |
| `‑d <file>` | Discard read queries not in this file                         |
| `‑e <file>` | Discard read queries in this file                             |
| `‑f <file>` | Discard reference dbs not in this file                        |
| `‑g <file>` | Discard reference dbs in this file                            |
| `‑h`        | Discard remaining suboptimal mappings                         |
| `‑i <int>`  | Discard remaining suboptimal mappings and discard any reads that have more remaining mappings than this |
| `‑j`        | Remove any unmatched nts in the very 3' end                   |
| `‑k`        | Output progress to standard output                            |

#### Example usage

```sh
parse_mappings.pl reads_vs_genome.arf -a 0 -b 18 -c 25 -i 5 \
  > reads_vs_genome_parsed.arf
```

---


### `perform_controls.pl`

#### Description

Performs a designated number of rounds of permuted controls (for details, see
Friedländer et al., Nature Biotechnology, 2008).

#### Input

The permutation controls estimate the number of false positives produced by a
`miRDeep2_core_algorithm.pl` run.
The input to `perform_controls.pl` should be

* a file containing the exact command line used to initiate the
  `miRDeep2_core_algorithm.pl` run,
* the structure file input to `miRDeep2_core_algorithm.pl`, and
* the desired rounds of controls.

#### Output

* A file in `.mrd` format.

The output of each control run is separated by a line `permutation integer`.
The mean number of entries output by the control runs gives an estimate of the
false positives produced. The further contents (besides the number of entries)
of the `.mrd` output by `perform_controls.pl` is not biologically meaningful.

#### Options

| option | description               |
|--------|---------------------------|
| `‑a`   | Output progress to screen |

#### Example usage

```sh
perform_controls.pl line potential_precursors.str 100 \
  > output_controls.mrd
```

---


### `permute_structure.pl`

#### Description

In a file output by RNAfold, each entry can be partitioned into an 'id' part
and an 'other' part, consisting of the dot-bracket structure, sequence, mfe
etc. This scripts reads all 'id' parts into a hash and pairs them with 'other'
parts from random entries. This is used by the `perform_controls.pl` script.

#### Input

* An RNAfold output file.

#### Output

* An RNAfold output file with IDs moved to random entries.

#### Example usage

```sh
permute_structure.pl potential_precursors.str \
  > potential_precursors_permuted.str
```

---


### `prepare_signature.pl`

#### Description

Prepares the signature file to be input to the `miRDeep2_core_algorithm.pl`
script.

#### Input

* A FASTA file with deep sequencing reads and
* a FASTA file with precursors.

#### Output

* A signature file in `.arf` format.

#### Options

| option      | description                                                   |
|-------------|---------------------------------------------------------------|
| `‑a <file>` | FASTA file with the sequences of known mature miRNAs for the species. These sequences will not influence the miRDeep scoring, but will subsequently make it easy to estimate sensitivity of the run. |
| `‑b`        | Output progress to screen                                     |

#### Example usage

```sh
prepare_signature.pl reads_collapsed.fa potential_precursors.fa \
  -a miRBase_this_species.fa > signature.arf
```

---


### `rna2dna.pl`

#### Description

Substitutes `u`s and `U`s to `T`s.
This is useful since `bowtie` does not match `U`s to `T`s.

#### Input

* A FASTA file.

#### Output

* A substituted FASTA file.


#### Example usage

```sh
rna2dna.pl reads_RNA_alphabet.fa > reads_DNA_alphabet.fa
```

---


### `select_for_randfold.pl`

#### Description

This script identifies potential precursors whose structure is basically
consistent with Dicer recognition.
Since running randfold is time-consuming, it is practical only to estimate
p-values for those potential precursors that actually fold into hairpin
structures.

#### Input

* An ARF file with the read signatures and
* an RNAfold output file with the structures of the potential miRNA precursors.

#### Output

* A list of ids, separated by newlines.

#### Example usage

```sh
select_for_randfold.pl signature.arf potential_precursors.str \
  > potential_precursors_for_randfold.ids
```

---


### `survey.pl`

#### Description

Surveys miRDeep2 performance at score cut-offs from -10 to 10.

#### Input

Default input is

* a `.mrd` file output by the `miRDeep2_core_algorithm.pl` script.

#### Output

* A .csv file with performace statistics.

#### Options

| option      | description                                         |
|-------------|-----------------------------------------------------|
| `‑a <file>` | file outputted by controls                          |
| `‑b <file>` | mature miRNA FASTA reference file for the species   |
| `‑c <file>` | signature file                                      |
| `‑d <int>`  | read stack height necessary for triggering excision |

#### Example usage

```sh
survey.pl output.mrd -a output_controls.mrd -b miRBase_this_species.fa \
  -c signature.arf -d 2 > survey.csv
```

---


### `convert_bowtie_output.pl`

#### Description

It converts a `bowtie` `bwt` mapping file to a `mirdeep` `arf` file.

#### Input

* A file in `bwt` format.

#### Output

* A file in `mirdeep` `arf` format.

---


### `bwa_sam_converter.pl`

#### Description

It converts a `bwa` `sam` mapping file to a `mirdeep` `arf` file.

#### Input

* A `bwa` created file in `sam` format.

#### Output

* A file in `mirdeep` `arf` format.

---


### `samFLAGinfo.pl`

#### Description

It gives information about the `bwa` FLAG in a `bwa` created mapping file in
`sam` format.

#### Input

* A FLAG number created by `bwa`.

#### Output

* Information about the alignment created by `bwa`.

---


### `clip_adapters.pl`

#### Description

Removes 3' end adaptors from deep sequenced small RNAs.
The script searches for occurrences of the six first nucleotides of the adapter
in the read sequence, starting after position 18 in the read sequence (so the
shortest clipped read will be 18 nts). If no matches to the first six nts of
the adapter are identified in a read, the 3' end of the read is searched for
shorter matches to the 5 to 1 first nts of the adapter.

#### Input

* A FASTA file with the deep sequencing reads and
* the adapter sequence (both in RNA or DNA alphabet).

#### Output

* A FASTA file with the clipped reads.

FASTA IDs are retained. If no matches to the adapter prefixes are identified in
a given read, the unclipped read is output.

#### Example usage

```sh
clip_adapters.pl reads.fa TCGTATGCCGTCTTCTGCTTGT > reads_clipped.fa
```

#### Notes

It is possible to clip adapters using more sophisticated methods. Users are
encouraged to test other methods with the miRDeep2 modules.

---


### `sanity_check_genome.pl`

#### Description

It checks the supplied genome FASTA file for its correctness.
Identifier lines are not allowed to contain whitespaces and must be unique.
Sequence lines are not allowed to contain characters others than
`A`, `C`, `G`, `T`, `N`, `a`, `c`, `g`, `t`, or `n`.

#### Input

* A genome file in FASTA format

---


### `sanity_check_mapping_file.pl`

#### Description

It checks the supplied mapping file for its correctness.
Each line in the file must be in the ARF format.

#### Input

* A mapping file in ARF format.

---


### `sanity_check_mature_ref.pl`

#### Description

It checks the supplied `mature_miRNA` FASTA file for its correctness.
Identifier lines are not allowed to contain whitespaces and must be unique.
Sequence lines are not allowed to contain characters others than `A`, `C`, `G`,
`T`, `N`, `a`, `c`, `g`, `t`, or `n`.

#### Input

* A mature miRNA file in FASTA format.

---


### `sanity_check_reads_ready.pl`

#### Description

It checks the supplied reads file for its correctness.
Each identifier line must have the format of '>name_uniqueNumber_xnumber` e.g.
`>xyz_1_x20`. See also file `format_descriptions.txt` for more detailed
informations.

#### Input

* A mapping file in ARF format.

---


### `extract_miRNAs.pl`

#### Description

Extracts mature and precursor sequences from miRBase fasta files for 
species of interest.

#### Input

* A fasta file from miRBAase
* One or more species three letter code abbreviations

#### Output 
* A fasta file in a proper format usable by quantifier.pl and miRDeep2.pl.
* Multiline sequences from input files are put on a single line and MacOS and Windows linebreaks/carriage returns are removed

#### Example usage

```sh
extract_miRNAs.pl mature_miRBase.fa hsa > mature_hsa.fa
extract_miRNAs.pl hairpin_miRBase.fa hsa > hairpin_hsa.fa
extract_miRNAs.pl mature_miRBase.fa mmu,chi > mature_other.fa
```
