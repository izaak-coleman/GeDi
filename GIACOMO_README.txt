Once you've compiled GeDi (see README.txt), you can run GeDi with the following

For the chromosome 22 dataset:
./GeDi -c chr22 -v 500 -t 32 -i chr22_input_list.txt -x \
/data/ic711/bowtie_indexes/ucsc.hg19.fasta -o chr22_output 

For the chromosome 22 dataset:
./GeDi -c chr17 -v 500 -t 32 -i chr17_input_list.txt -x \
/data/ic711/bowtie_indexes/ucsc.hg19.fasta -o chr17_output

You can see an explaination of all GeDi's parameters by running:
./GeDi --help

To briefly explain the required parameters here:

-x /path/to/bowtie/hg/index
This must be the path to the bowtie index of the reference genome that
GeDi will use for alignment of healthy data. It's the same path you would
give to the -x parameter of bowtie2

-v 500
expected coverage of the dataset

-t 32
launch with 32 threads

-i chr22_input_list.txt
This is a list of absolute paths to the tumour-control paired datasets. 
On the same line, each path must be followed by a comma and then a symbol, 
either C or T, denoting whether the file should be read into GeDi as
either control or tumour data respectively. For example, the file contents
may look like this:
/data/sorted_in_chr22_GOLDEN_blood.fastq.gz,C
/data/sorted_in_chr22_GOLDEN_primary.fastq.gz,T

-c chr22
The targeted chromosome in the analysis. 
Note that the string input into -c, in this case "chr22", must match
the fasta header for the corresponding chromosome in the
reference genome used to construct the bowtie2 index. In this case
it must be ">chr22", not ">22".

-o chr22_output
Prefix of GeDi's output filenames.
The sam and fastq file generated contain the healthy consensus
sequences used to align to bowtie2. These can be deleted.
The file chr22_output.SNV_results contains the list of SNV calls made
by GeDi.

