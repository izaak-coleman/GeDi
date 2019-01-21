# GeDi: Generalised suffix array-based Direct SNV caller

Applying suffix arrays to increase the repertoire of detectable SNVs 
in tumour genomes.

GeDi, or **Ge**neralised suffix array-based **Di**rect SNV caller,
detects somatic Single Nucleotide Variants (SNVs) within paired tumour-control
NGS datasets. GeDi is capable of detecting SNVs reference- and mapping-free.
To acheive this, it compares the input sequencing data by
construction of multiple suffix arrays from which SNVs can be directly detected;
within the arrays, suffixes containing SNVs cluster into intervals enriched with
tumour-derived data.

This work has been published in XXXX.


# Download
Either clone this repository, or download the zip file containing GeDi source
code and unzip the source code directory. To clone the GeDi repo, issue
the following on the command-line:

```
git clone https://github.com/izaak-coleman/GeDi
```
Move the directory named `GeDi` (source code directory) 
into a location where you are comforable with it remaining; 
once GeDi is compiled, the source code directory cannot be moved 
without re-compilation (see compilation section below).


# Compilation
To reduce the required dependencies for GeDi, we included many of the
libaries GeDi requires within its source code directory, which is where
the GeDi executable will be built. As a consequence, once GeDi is compiled
from within its source code directory, it cannot be moved without recompilation. 

Nevertheless, a few dependencies and requirements (listed below)
remain to successfully compile GeDi.
  - GeDi can be run only on linux machine with x86-64 architecture. 
  - GNU `make` must be installed 
  - `cmake` must be installed
  - `g++` version 4.8.1 or greater
  - boost program options library. On CentOS, install with `yum install
  boost-devel.x86_64`, on Ubuntu with `sudo apt-get install libboost-all-dev`.
  - A bowtie2 index of a reference genome. For example, `ucsc.hg19.fasta`
    can be build by installing bowtie2 and issuing `bowtie2-build
    ucsc.hg19.fasta ucsc.hg19.fasta`, see [bowtie2 reference manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
    for further details.

Once the depenendancies have been installed and GeDi's source code
directory is placed in a prefered location, compile GeDi with the following:
```
cd GeDi
bash install.sh
```
Compilation of GeDi will now initiate. Once compilation is complete,
test that GeDi compiled succesfully by executing
```
./GeDi --help
```
which should display GeDi's help message.


# Running GeDi on example data

Within GeDi's source code directory is an `example_data` subdirectory
that will provide an example of how to run GeDi in default mode.

Within `example_data` are the following three files:
  - `example_dataset.file_list`
  - `example_control_data.fastq.gz`
  - `example_tumour_data.fastq.gz`

The tumour and control `fastq` files form the example paired tumour-control
NGS dataset. `example_dataset.file_list` comprises a list of `fastq` files 
that form the dataset. This file is read by GeDi to load the dataset. Each line of 
`example_dataset.file_list` has the structure `filename, type`, where `filename`
is the absolute path to a fastq file, and `type` is either `C` or `T`, denoting
whether the fastq file derives from control or tumour tissue respectively.
Hence, `example_dataset.file_list` contains the following:
```
./example_data/example_control_data.fastq.gz,C
./example_data/example_tumour_data.fastq.gz,T
```
Note that relative paths to fastq files are given, but the user
should *always use absolute paths in their `.file_list`*

To run GeDi on the example data, issue the following command:
```
./GeDi -c chr22 -v 30 -t 4 -i example_data/example_datafiles.file_list -x /path/to/bowtie2_index/ucsc.hg19.fasta -o example_data
```
Once GeDi has finished executing, it will output three files:
 - `example_data.SNV_results`
 - `example_data.fastq.gz`
 - `example_data.sam`

The `.fastq.gz` and `.sam` file contain the healthy consensus sequences aligned
to the reference genome as a proxy used in the calculation of SNV coordinates;
they can be deleted. `example_data.SNV_results` contains the SNV calls made
by GeDi for the example dataset. This file has the following content:
```
Mut_ID	Type	Chr	Pos	Normal_NT	Tumor_NT
0	SNV	chr22	19613299	T	A

```
When analysing the toy dataset, GeDi detected and called a single SNV on
chromosome 22 at position 19613299, with control variant T and tumour variant
A. The column headers of an output `.SNV_results` file describe the following:
  - `MuT_ID` assigns to each all an arbitrary unique number. 
  - `Type` specifies the type of variant GeDi called - current version of GeDi only calls SNVs.
  - `Chr` specifies the chromosome the SNV resides.
  - `Pos` specifies the 1-indexed chromosome position of the SNV.
  - `Normal_NT` specifies the control variant of the SNV.
  - `Tumour_NT` specifies the tumour variant of the SNV.

# Running GeDi on user data
By following the above example and making the appropriate substitutions
a user will be able to analyse their own data with GeDi.

To run GeDi in default mode, the main substitution is to ensure
the `.file_list` file contains the list of `fastq` files comprising
the users own dataset.

Although issuing `./GeDi --help` on the command line will print information
regarding the arguments that must be passed to GeDi's parameters and parameter
function. Here, we explain what must be passed to the most important of GeDi's
parameters (each of which was specified in the above example):
 - `-i` must be given the path the users `.file_list`.
 - `-c` If specified, GeDi will filter out all SNV calls that do not
        reside within the specified chromosome. Note that, this argument
        should exactly match the fasta header within the reference genome
        used for alignment. If not specified, GeDI will report
        all SNVs called regardless of the chromosome they reside in.
 - `-v` Expected or average coverage of the dataset.
 - `-t` number of threads GeDi will execute with during parallel sections
 - `-x` the location and prefix of the bowtie2 human reference genome index.
        Note this argument passed to this parameter should be identical to the
        argument one would pass to `-x` if running bowtie2. See
        [bowtie2 reference manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) 
        for a further explaination of bowtie2's `-x` parameter.
 - `-o` prefix of the output files produced by GeDi. This can include a file
        path, and GeDi will write files with the prefix at the specified location.

