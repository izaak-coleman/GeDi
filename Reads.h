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
  unsigned int maxLen;
  // Index correspondence between Reads/Phreds
  std::vector<std::string> HealthyReads
  std::vector<std::string> HealthyPhreds;
  std::vector<std::string> TumourReads;
  std::vector<std::string> TumourPhreds;
  void parseInputFile(std::string const& inputFile, 
           std::vector<file_and_type> &datafiles);

  void loadFastqRawDataFromFile(std::string const& filename, 
                              std::vector<std::string> & processed_reads,
                              std::vector<std::string> &processed_phreds);

  void qualityProcessRawData(std::vector<fastq_t> *r_data, 
                            std::vector<std::string> *processed_reads,
                            std::vector<std::string> *processed_phreds,
                            int from, int to, int tid);
public:

  ReadPhredContainer(std::string const& inputFile);

  char baseQuality(int index, int tissue, int pos);

  std::string & getPhredString(int index, int tissue);

  std::string::iterator returnStartIterator(Suffix_t &suf);
  
  std::string::iterator returnEndIterator(Suffix_t &suf);

  std::string returnSuffix(Suffix_t &suf);

  std::string & getReadByIndex(int index, int tissue);

  unsigned int getSize(bool tissueType);

  unsigned int maxLengthRead();

  int getMinSuffixSize();
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
