#GeDi: Generalised suffix array-based Direct SNV caller

Applying suffix array to increase the repertoire of detectable SNVs 
in tumour genomes.

GeDi, or *Ge*neralised suffix array-based *Di*rect SNV caller,
detects somatic Single Nucleotide Variants (SNVs) within paired tumour-control
NGS datasets. GeDi is capable of detecting SNVs reference- and mapping-free.
To acheive this, it compares the input sequencing data by
construction of multiple suffix arrays from which SNVs can be directly detected;
within the arrays, suffixes containing SNVs cluster into intervals enriched with
tumour-derived data.

This work has been published in XXXX.


# Download
Either clone this repository, or download the zip file containing GeDi source
code and unzip the source code directory. For example:

```
https://github.com/izaak-coleman/GeDi
```

Move the source code directory into a location where you are comforable with 
it remaining; once GeDi is compiled, the source code directory cannot be moved 
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


# Running GeDi

