/*
Reads.h
Author: Izaak Coleman
*/


#ifndef READS_H
#define READS_H

#include <string>
#include <iostream>
#include <vector>
#include <utility>


#include "util_funcs.h"

struct fastq_t {
  std::string const id, seq, qual;
  fastq_t(std::string i, std::string s, std::string q): id(i), seq(s), qual(q){}
};

struct file_and_type {
  const std::string file;
  const bool tissue;
  file_and_type(std::string f, bool t): file(f), tissue(t){};
};

class ReadPhredContainer {
private:
  const int MIN_SUFFIX_SIZE;
  const int N_THREADS;

  // Length of the longest read
  unsigned int maxLen;          
 
  // Containers storing the reads and their associated phred strings.
  // Corresponding reads and phred strings share the same indicies in their
  // respective containers.
  std::vector<std::string> HealthyReads;
  std::vector<std::string> HealthyPhreds;
  std::vector<std::string> TumourReads;
  std::vector<std::string> TumourPhreds;

  void parseInputFile(std::string const& inputFile, 
           std::vector<file_and_type> &datafiles);
  // Parses the input file containing the file paths, names and
  // data types of the data files. 
  // Separates the fields and stores them in datafiles

  void loadFastqRawDataFromFile(std::string const& filename, 
                              std::vector<std::string> & processed_reads,
                              std::vector<std::string> &processed_phreds);
  // Loads the fastq data in filename into processed_read and processed_phreds,
  // processing them with qualityProcessRawData()

  void qualityProcessRawData(std::vector<fastq_t> const& r_data, 
                            std::vector<std::string> & processed_reads,
                            std::vector<std::string> & processed_phreds);
  // Performs the following in the order given:
  //  -- Discards reads with a greater than QUALITY_THRESH percentage
  //     of bases that have a phred score less than PHRED_20
  //  -- Removes 'N' characters by spliting a read into substrings that flank
  //     'N's. Resulting substrings with length < MIN_SUFFIX_SIZE are discarded.
  //  -- Appends TERM_CHAR to each read, for GSA construction.
public:

  ReadPhredContainer(std::string const& inputFile, int n);

  char baseQuality(int index, int tissue, int pos);
  // Returns Healthy/TumourPhreds[index][pos]

  std::string & getPhredString(int index, int tissue);

  std::string::iterator returnStartIterator(Suffix_t &suf);
  // Returns iterator to start of suffix represented by suf 
  std::string::iterator returnEndIterator(Suffix_t &suf);
  // Returns iterator to end of suffix represented by suf
  std::string returnSuffix(Suffix_t &suf);
  // Returns string of suffix represented by suf
  std::string & getReadByIndex(int index, int tissue);

  unsigned int getSize(bool tissueType);

  unsigned int maxLengthRead();

  int getMinSuffixSize();

  void free();
  // Function deallocates the memory allocated for containers:
  // - HealthyReads
  // - HealthyPhreds
  // - TumourReads
  // - TumourPhreds
};
#endif
/*
void printAllReads();
// prints all reads to file reads_after_icsmufin.txt
void printRemainingReads(std::string const& filename);
// prints the tuples of the read once loaded from files
void printReadsAndId(int from, int to, int step) const;
// iterates through the read containers over interval [from, to) 
// with a step size of step. Prints the read, and for each printed
// read, prints its id
void printReads();
// Prints the loaded reads REMOVE THIS FUNCTION

 void writeContainer(std::vector<std::string> const& c, std::string const& fname);
 std::string performDistalTrim(std::string & s);
 // Function performs a distal trim of length distal_trim_len.
 // This is currently set to distance DISTAL_TRIM -> need
 // to switch to user input
*/
