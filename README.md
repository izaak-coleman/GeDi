# GeDi: Generalised suffix array-based Direct SNV caller

Applying suffix arrays to increase the repertoire of detectable SNVs 
in tumour genomes.

**Ge**neralised suffix array-based **Di**rect SNV caller, or GeDi,
detects somatic Single Nucleotide Variants (SNVs) within paired tumour-control
NGS datasets. GeDi is capable of detecting SNVs reference- and mapping-free.
To achieve this, it compares the input sequencing data by
construction of multiple suffix arrays from which SNVs can be directly detected;
within the arrays, suffixes containing SNVs cluster into intervals enriched with
tumour-derived data. GeDi uses mapping only to determine SNV coordinates.
To do this, healthy-tissue derived sequences are mapped to the reference
genome as a proxy from which SNV coordinates are calculated. Accordingly,
not only does GeDi align SNVs reference- and mapping-free, it never
aligns tumour data to the reference genome. This enables GeDi
to detect SNVs at complex variant loci, such as sites of hypermutation,
with high sensitivity.

This work has been published in Coleman, I., Corleone, G., Arram, J. et al. GeDi: applying suffix arrays to increase the repertoire of detectable SNVs in tumour genomes. BMC Bioinformatics 21, 45 (2020). https://doi.org/10.1186/s12859-020-3367-3.


# Download
Either clone this repository, or download and unzip the zip file containing 
GeDi's source code. To clone the GeDi repo, issue the following on the command-line:

```
git clone https://github.com/izaak-coleman/GeDi
```
Then move the directory named `GeDi` (source code directory) 
to a location you are comfortable with it remaining; 
once GeDi is compiled, the source code directory cannot be moved 
without re-compilation (see compilation section below).


# Compilation
To reduce the number of dependencies requiring manual installation by the user, 
we included many of GeDi's dependencies within its source code directory. 
GeDi will be built directly within its source code directory and as a consequence of
the included dependencies, once GeDi is compiled the source code directory cannot
be moved without recompilation of GeDi. 

Nevertheless, a few dependencies and requirements (listed below)
remain to successfully compile GeDi.
  - GeDi can be run only on linux machine with x86\_64 architecture. 
  - GNU `make` must be installed 
  - `cmake` must be installed
  - `g++` version 4.8.1 or greater must be installed 
  - boost program options library must be installed. On CentOS, install with `yum install
  boost-devel.x86_64`, on Ubuntu with `sudo apt-get install libboost-all-dev`.
  - A bowtie2 index of a reference genome must be present. For example, an index
  of `ucsc.hg19.fasta` (Human Reference Genome 19) can be built by installing 
  bowtie2 and issuing command `bowtie2-build
  ucsc.hg19.fasta ucsc.hg19.fasta`, see 
  [bowtie2 reference manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
  for further details.

Once the dependencies have been installed and GeDi's source code
directory is placed in a preferred location, compile GeDi with the following:
```
cd GeDi
bash install.sh
```
`sudo` may be required. Compilation of GeDi will now initiate. Once compilation is complete,
test that GeDi compiled successfully by executing
```
./GeDi --help
```
which should display GeDi's help message.


# Running GeDi on example data

Within GeDi's source code directory is an `example_data` subdirectory.
It will be used in this section to provide an example of how to analyse a paired
tumour-control NGS dataset with GeDi in default mode. 

Within `example_data` are the following three files:
  - `example_dataset.file_list`
  - `example_control_data.fastq.gz`
  - `example_tumour_data.fastq.gz`

The tumour and control `fastq` files form the example paired tumour-control
NGS dataset. We assume that whatever fasta preprocessing (e.g quality analysis,
de-duplication and trimming) deemed necessary by the user has already been applied. 
`example_dataset.file_list` contains the list of `fastq` file names
that form the dataset. `example_dataset.file_list` will be read by GeDi to 
load the dataset. Each line of `example_dataset.file_list` has the structure `filename, type`, where `filename`
is the absolute path to a fastq file, and `type` is either `C` or `T`, denoting
whether the fastq file derives from control or tumour tissue respectively.
Hence, `example_dataset.file_list` contains the following:
```
./example_data/example_control_data.fastq.gz,C
./example_data/example_tumour_data.fastq.gz,T
```
Note that relative paths to fastq files are given, but the user
should *always use absolute paths in their* `.file_list`.

To run GeDi on the example data, issue the following command:
```
./GeDi -c chr22 -v 30 -t 4 -e 0 -i example_data/example_dataset.file_list -x /path/to/your/bowtie2_index/of/human_ref_genome/ucsc.hg19.fasta -o example_data
```
Once GeDi has finished executing, it will output three files:
 - `example_data.SNV_results`
 - `example_data.fastq.gz`
 - `example_data.sam`

The `.fastq.gz` and `.sam` files contain the healthy-tissue derived sequences
mapped to the reference genome as a proxy; they can be deleted. 
`example_data.SNV_results` contains the SNV calls made
by GeDi for the example dataset. This file has the following content:
```
Mut_ID	Type	Chr	Pos	Normal_NT	Tumor_NT
0	SNV	chr22	19613299	T	A

```
When analysing the toy dataset, GeDi detected and called a single SNV on
chromosome 22 at position 19613299, with control variant T and tumour variant
A. The column headers of an output `.SNV_results` file describe the following:
  - `MuT_ID` is an arbitrary unique number assigned to each called SNV.
  - `Type` specifies the type of variant GeDi called - current version of GeDi only calls SNVs.
  - `Chr` specifies the chromosome the called SNVs reside.
  - `Pos` specifies the 1-indexed chromosome position of the called SNVs.
  - `Normal_NT` specifies the control variant nucleotide of the called SNVs.
  - `Tumour_NT` specifies the tumour variant nucleotide of the called SNVs.

# Running GeDi on user data
By following the above example and making the appropriate substitutions
a user will be able to analyse their own data with GeDi.

Accordingly, to run GeDi in default mode, the main substitution is to ensure
the `.file_list` file contains the list of `fastq` files comprising
the users own dataset.

Although issuing `./GeDi --help` on the command-line will print information
regarding the arguments that must be passed to GeDi's parameters as well
as parameter function. Here, we explain what must be passed to the most important of GeDi's
parameters (each of which was specified in the above example):
 - `-i` The path and file name of the users `.file_list`.
 - `-c` If specified, GeDi will filter out all SNV calls that do not
        reside within the specified chromosome. Note that the passed argument
        should *exactly* match the fasta header of the desired chromosome 
        within the reference genome used for alignment. If not specified, GeDi will report
        all SNVs called regardless of the chromosome they reside in.
 - `-v` Expected or average coverage of the dataset.
 - `-t` Number of threads GeDi will execute with during parallel sections
 - `-x` The location and prefix of the bowtie2 index of Human Reference Genome.
        Note the argument passed to this parameter should be identical to the
        argument one would pass to `-x` if running bowtie2. See
        [bowtie2 reference manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) 
        for a further explanation of bowtie2's `-x` parameter.
 - `-o` Prefix of the output files produced by GeDi. This can include a file
        path, and GeDi will write files with the prefix at the location
        specified by the path.
 - `-e` Expected contamination: Proportion of tumour dataset sequencing reads
        expected to have derived from healthy tissue.

The function of remaining parameters can be found by issuing `GeDi --help` on
the command-line. For further details consult XXX.
