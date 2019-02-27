# miRDeep2 changelog

## Version 2.0.1.2

* fixed uninitialized value errors/warnings

* rephrased available species message

* fixed Travis CI (still used to use `~/.bash_profile` rather than `~/.bashrc`)

---

## Version 2.0.1.1

* improved detection of installation directory & Rfam indices

* instead of fractions, percantages are used in read statistics

* several minor improvemnts to the MRD output format

* negative repeat count warnings have been removed

* HTML output does not contain (now dead) links to (old) Rajewsky lab webpage
  anymore

* `install.pl` now modifies `~/.bashrc` rather than `~/.bash_profile`

* `miRDeep2.pl` now handles changed RNAfold output

---

## Version 2.0.1.0

* fixed bug resulting in negative repeat count due to too small `-l` value

* braces in regex were escaped [deprecated by upstream]

* miRDeep2 was licensed under GPL v3

* the documentation was improved and converted to Markdown

* minor improvements and code cleanup

* the tutorial is now run via Travis CI on Linux & Mac to test the build

* Sean Eddys squid libray was removed from the repository

---

## Version 2.0.0.8

* The `install.pl` script has been updated to download automatically the latest
  `PDF::API2` version and the latest `bowtie1` version.

* Source code has been cleaned up and pathes have been adapted. So when using
  this version please do a new installation of miRDeep2 and do not just copy
  scripts over old version. Otherwise your installation will likely be broken.

* RNAfold and randfold now compile again on MacOSX 10.10.5 due to adaptation of
  the `fold.c` source file in these packages.

* minor issues removed

* code is cleaned up and installation directories are also organzied more reasonable

* Sean Eddys squid libray is included in the package already since the download
  server is sometimes not to be working

---


## Version 2.0.0.6

* corrected bed file coordinates, off by one error fixed

* the zebrafish switch works now in the species specification

* novel detected miRNA sequences by the miRDeep2 module are now stored in
  separate files

* in the mapper module all tmp files are deleted now when the specific option
  is given

* install.pl script has been updated in regards of program download links

* additional file checks were included for the mapper module when using options
  `a`, `b`, `c`, or `e`

* novel detected miRNA sequences and bed coordinates are output automatically
  now

---


## Version 2.0.0.5

* The `quantifier.pl` has a new switch `-U` that causes the usage of only
  uniquely mapping reads to reference precursors.

* The `mapper.pl` reports now the number and percentage of mapped and unmapped
  reads from the collapsed reads file to a reference genome.

* miRDeep2 now aborts when no precursors were excised.

* A BED file of all miRDeep2 detected precursors is created in the end.

* An additional check of the ARF file ids and genome file ids was added to make
  sure that ids in both files really match.

* An additional check for the presence of the `Rfam` indices was added.

* A bug was fixed that occurred when the miRDNA IDs were too long.

* The HTML file in the end does not abort anymore when the `Rfam` index is
  missing.

---


## Version 2.0.0.4

* The `quantifier.pl` now calculates correctly the remaining read counts.

* Reads that map to more than one precursor are now added in a weighed way to
  each precursor. The normalization of read counts in the `quantifier.pl` has
  changed a bit. Each miRNA read count is now divided by the total number of
  sequenced miRNA read counts and then multiplied by 10E9.

* The `--strata` option is added to mapping of reads to the miRBase precursors
  in the quantifier module.

* The miRDeep2 html output has been updated so that for known miRNAs the field
`mature miRBase miRNA` now contains the miRNA that matches the predicted mature
sequence or in case that the predicted star sequence matches a known miRNA this
is then shown instead.

* The `install.pl` script can also use `curl` if `wget` is not available on the
  machine. This is useful for MACOSX users.

* The `mapper.pl` now also determines the number of cores on Mac machines
  correctly.

* A bug has been fixed in the `prepare_signature.pl` file. It crashed when no
 ` known_mature_mirna` file was supplied.

* RNAfold2 can be used as well now if already installed. Otherwise version
  1.8.4 of the Vienna package is downloaded when running the `install.pl`
  script.
