#ifndef READS_H
#define READS_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>


#include "util_funcs.h"

struct fastq_t {      // Struct only read needs to know about
  std::string id, seq, qual;
};

struct file_and_type {
  std::string first;    // filename
  bool second;          // data set
};


class ReadPhredContainer {
  // This class stores all the reads for the dataset. 
  // Functions are allowed direct access to reads. 

private:
  int minimum_suffix_size;
  int distal_trim_len;
  unsigned int maxLen;
  std::vector<std::string> HealthyReads;  // Container for healthy dataset
  std::vector<std::string> TumourReads;   // Container for cancer dataset 
  std::vector<std::string> HealthyPhreds;  // Read and phred containers correspond by index
  std::vector<std::string> TumourPhreds;


  void parseInputFile(std::string const& inputFile, std::vector<file_and_type> &datafiles);
  // Parses inputFile, extraction data file path and tissue subset

  void loadFastqRawDataFromFile(std::string filename, 
                              std::vector<std::string> & processed_reads,
                              std::vector<std::string> &processed_phreds);
  // Function loads DNA reads from filename.fasta.gz
  // either HealthyReads/TumourReads dep. on passes pointer

  void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                            std::vector<std::string> *processed_reads,
                            std::vector<std::string> *processed_phreds,
                            int from, int to, int tid);
  // **
  // Parameters are pointer based to allow compatibility with 
  // <thread> for potential multithreading. 
  // Quality processing each read involves:
  // 1) Discarding reads where percentage of positions having a 
  // phred score < PHRED_20 is > QUALITY_THRESH. 
  // 2) Removing 'N' chars, by slicing read into substrings that exclude any
  // 'N''s. Resultant substrings with length < min_suffix_size are discarded. 

  std::string performDistalTrim(std::string & s);
  // Function performs a distal trim of length distal_trim_len.
  // This is currently set to distance DISTAL_TRIM -> need
  // to switch to user input

public:
 // these arrays have a 1:1 mapping with the HealthyReads, TumourReads arrays
 // where their values express whether the read at the same index in 
 // the Healthy/TumourReads arrays if of LEFT or RIGHT type

 ReadPhredContainer(std::string const& inputFile);
 // Constructor for loading and processing reads

 char baseQuality(int index, int tissue, int pos);
 std::string & getPhredString(int index, int tissue);

 void printAllReads();
 // prints all reads to file reads_after_icsmufin.txt

  void writeContainer(std::vector<std::string> const& c, std::string const&
      fname);

 ~ReadPhredContainer();
 // Deletes the reads


 void printReads();
 // Prints the loaded reads REMOVE THIS FUNCTION

 unsigned int getSize(bool tissueType);
 // returns the size HealthyReads, or TumourReads dep. on tissueType

 unsigned int maxLengthRead();
 // returns maxLen


 std::string::iterator returnStartIterator(Suffix_t &suf);
// Function locates the read corresponding to suf.read_id and 
// Sets a pointer in that read starting at suf.offset
 
 std::string::iterator returnEndIterator(Suffix_t &suf);
// Function returns a pointer to the end of the read corresponding to 
// suf.read_id

 std::string returnSuffix(Suffix_t &suf);
// Function returns the suffix that the suffix_t represents

 std::string & getReadByIndex(int index, int tissue);
 // This function returns the read by adress from the given 
 // Reads vector

 int getMinSuffixSize();
 // returns the min suffix size of suffixes in the suffix array
 // as specified by the user


 void printRemainingReads(std::string const& filename);
 // prints the tuples of the read once loaded from files

 void printReadsAndId(int from, int to, int step) const;
 // iterates through the read containers over interval [from, to) 
 // with a step size of step. Prints the read, and for each printed
 // read, prints its id

};
#endif
