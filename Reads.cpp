/*
Reads.cpp
Author: Izaak Coleman
*/


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstring>
#include <omp.h>
#include <zlib.h>

#include "kseq.h"
#include "util_funcs.h"
#include "Reads.h"
#include "benchmark.h"

KSEQ_INIT(gzFile, gzread);    // initialize .gz parser

using namespace std;
static const int      TERM_CHAR_CORRECTION = 1;
static const int      _MIN_SUFFIX_SIZE = 30;
static const double   QUALITY_THRESH = 0.1;
static const char     PHRED_20 = '5';
static const          string REMOVED_TOKENS = "N";
static const          string TERM_CHAR = "$"; // for GSA
static const          string HEALTHY_DATA = "H";
static const          string TUMOUR_DATA  = "T";

ReadPhredContainer::ReadPhredContainer(string const& inputFile, int n):
MIN_SUFFIX_SIZE(_MIN_SUFFIX_SIZE), N_THREADS(n) {
//START(ReadPhredContainer_ReadPhredContainer);
  vector<file_and_type> datafiles;
  maxLen = 0;
  parseInputFile(inputFile, datafiles);
  cout << "Loaded " << datafiles.size() << " data files." << endl;
#pragma omp parallel for num_threads(datafiles.size())
  for(int i=0; i < datafiles.size(); i++) {
#pragma omp critical 
    {
      cout << "Extracting data from " << datafiles[i].file << "..." << endl;
    }
    // thread-local containers
    vector<string> reads, phreds; 
    vector<unsigned int> lengths;
    loadFastqRawDataFromFile(datafiles[i].file, reads, phreds, lengths);
#pragma omp critical 
    {
      for (unsigned int l : lengths) {
        if (l > maxLen) maxLen = l;
      }

      if (datafiles[i].tissue == HEALTHY) {
        HealthyReads.insert(HealthyReads.end(), reads.begin(), reads.end());
        HealthyPhreds.insert(HealthyPhreds.end(), phreds.begin(), phreds.end());
      } else {
        TumourReads.insert(TumourReads.end(), reads.begin(), reads.end());
        TumourPhreds.insert(TumourPhreds.end(), phreds.begin(), phreds.end());
      }
    }
  }
  cout << "MAX READ LEN: " << maxLengthRead() << endl;
  cout << "Reads: " << TumourReads.size() + HealthyReads.size() << endl;
  cout << "Phreds: " << TumourPhreds.size() + HealthyPhreds.size() << endl;
//COMP(ReadPhredContainer_ReadPhredContainer);
}
void ReadPhredContainer::parseInputFile(string const& inputFile, 
                                        vector<file_and_type> &datafiles) {
//START(ReadPhredContainer_parseInputFile);
  ifstream sock;
  sock.open(inputFile.c_str());
  cout << "Gathering datafiles from " << inputFile << "." << endl;
  string file_string;
  while (getline(sock, file_string)) {
    vector<string> fields;
    split_string(file_string, ",\t ", fields);
    if (fields[1] == HEALTHY_DATA) {
      datafiles.push_back(file_and_type(fields[0], HEALTHY));
    } else if (fields[1] == TUMOUR_DATA) {
      datafiles.push_back(file_and_type(fields[0], TUMOUR));
    } else {
      cout << "Invalid input found in " << inputFile
           << ". Valid datatypes are either H or T." 
           << endl << "Program terminating." << endl;
      exit(1);
    }
    cout << "Input " << fields[0] << " as "
         << ((fields[1] == HEALTHY_DATA) ? "healthy" : "tumour") << " data "
         << endl;
  }
  sock.close();
//COMP(ReadPhredContainer_parseInputFile);
}

void ReadPhredContainer::loadFastqRawDataFromFile(string const& filename, 
                              vector<string> &processed_reads, 
                              vector<string> & processed_phreds,
                              vector<unsigned int> & lengths) {
//START(ReadPhredContainer_loadFastqFromRawDataFile);
  gzFile data_file;
  data_file = gzopen(filename.c_str(), "r");    // open stream to next fastq.gz 
  kseq_t *seq = kseq_init(data_file);           // init parser

  vector<fastq_t> fastq_elements;
  while (kseq_read(seq) >= 0) {
    if (seq->qual.l) {
      fastq_elements.push_back(fastq_t(seq->name.s, seq->seq.s, seq->qual.s));
    }
  }
  kseq_destroy(seq);
  gzclose(data_file);
  qualityProcessRawData(fastq_elements, processed_reads, processed_phreds, lengths);
//COMP(ReadPhredContainer_loadFastqFromRawDataFile);
}

//void ReadPhredContainer::qualityProcessRawData(vector<fastq_t> const& r_data, 
//                           vector<string> &processed_reads,
//                           vector<string> &processed_phreds){
//  processed_reads.reserve(r_data.size());
//  processed_phreds.reserve(r_data.size());
//  unsigned int elementsPerThread = r_data.size() / N_THREADS;
//  omp_set_num_threads(N_THREADS);
//#pragma omp parallel for 
//  for(unsigned int i = 0; i < r_data.size(); i += elementsPerThread) {
//    unsigned int to = i + elementsPerThread;
//    if (to > r_data.size()) {
//      to = r_data.size();
//    }
//    // local thread stores - used to avoid lock contention when copying
//    // to object data members
//    vector<string> reads, phreds;
//    vector<int> lengths; 
//
//    // begin quality check
//    for (unsigned int from = i; from < to; from++) {
//      double nLowQualBases = 0.0;
//      for(auto it = r_data[from].qual.begin(); it != r_data[from].qual.end(); it++) {
//        if (*it < PHRED_20) nLowQualBases++;
//      }
//      if ((nLowQualBases / r_data[from].qual.size()) > QUALITY_THRESH) continue;
//      // Remove 'N' chars
//      vector<string> read_substrs;
//      vector<string> phred_substrs;
//      int left_arrow{0}, right_arrow{0};
//      while((right_arrow = r_data[from].seq.find(REMOVED_TOKENS, left_arrow)) != string::npos) {
//        if (left_arrow == right_arrow) {
//          left_arrow++;
//          continue;
//        }
//        read_substrs.push_back(r_data[from].seq.substr(left_arrow, right_arrow - left_arrow));
//        phred_substrs.push_back(r_data[from].qual.substr(left_arrow, right_arrow - left_arrow));
//        left_arrow = right_arrow + 1;
//      }
//      read_substrs.push_back(r_data[from].seq.substr(left_arrow));
//      phred_substrs.push_back(r_data[from].qual.substr(left_arrow));
//      // Load data and store length of longest read
//      for (int i=0; i < read_substrs.size(); i++) {
//        if(read_substrs[i].length() >= MIN_SUFFIX_SIZE ) {
//          lengths.push_back(read_substrs[i].length());
//          reads.push_back(read_substrs[i] + TERM_CHAR);
//          phreds.push_back(phred_substrs[i]);
//        }
//      }
//    }
//#pragma omp critical 
//    {
//      for(int i=0; i < reads.size(); i++) {
//        processed_reads.push_back(reads[i]);
//        processed_phreds.push_back(phreds[i]);
//        if(lengths[i] > maxLen) maxLen = lengths[i];
//      }
//    }
//  }
//  processed_reads.shrink_to_fit();
//  processed_phreds.shrink_to_fit();
//}

void ReadPhredContainer::qualityProcessRawData(vector<fastq_t> const& r_data, 
                           vector<string> &processed_reads,
                           vector<string> &processed_phreds,
                           vector<unsigned int> &lengths) {
  processed_reads.reserve(r_data.size());
  processed_phreds.reserve(r_data.size());
  for(int i = 0; i < r_data.size(); i++) {
    // Discard low quality reads
    double nLowQualBases = 0.0;
    for(auto it = r_data[i].qual.begin(); it != r_data[i].qual.end(); it++) {
      if (*it < PHRED_20) nLowQualBases++;
    }
    if ((nLowQualBases / r_data[i].qual.size()) > QUALITY_THRESH) continue;
    // Remove 'N' chars
    vector<string> read_substrs;
    vector<string> phred_substrs;
    int left_arrow{0}, right_arrow{0};
    while((right_arrow = r_data[i].seq.find(REMOVED_TOKENS, left_arrow)) != string::npos) {
      if (left_arrow == right_arrow) {
        left_arrow++;
        continue;
      }
      read_substrs.push_back(r_data[i].seq.substr(left_arrow, right_arrow - left_arrow));
      phred_substrs.push_back(r_data[i].qual.substr(left_arrow, right_arrow - left_arrow));
      left_arrow = right_arrow + 1;
    }
    read_substrs.push_back(r_data[i].seq.substr(left_arrow));
    phred_substrs.push_back(r_data[i].qual.substr(left_arrow));
    // Load data and store length of longest read
    for (int i=0; i < read_substrs.size(); i++) {
      if(read_substrs[i].length() >= MIN_SUFFIX_SIZE ) {
        lengths.push_back(read_substrs[i].length());
        processed_reads.push_back(read_substrs[i] + TERM_CHAR);
        processed_phreds.push_back(phred_substrs[i]);
      }
    }
  }
}

unsigned int ReadPhredContainer::maxLengthRead() {
  return maxLen;
}
string::iterator ReadPhredContainer::returnStartIterator(Suffix_t &suf) {
//START(ReadPhredContainer_returnStartIterator);
  if(suf.type == HEALTHY) { 
    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnStartIterator() out of bounds " << endl;
      exit(1);
    }
    return HealthyReads[suf.read_id].begin() + suf.offset;
  }
  else {  // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnStartIterator() out of bounds " << endl;
      exit(1);
    }
    return TumourReads[suf.read_id].begin() + suf.offset;
  }
//COMP(ReadPhredContainer_returnStartIterator);
}

string::iterator ReadPhredContainer::returnEndIterator(Suffix_t &suf) {
//START(ReadPhredContainer_returnEndIterator);
  if (suf.type == HEALTHY) {
    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnEndIterator() out of bounds " << endl;
      exit(1);
    }
    return HealthyReads[suf.read_id].end();
  }
  else {  // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnEndIterator() out of bounds " << endl;
      exit(1);
    }
    return TumourReads[suf.read_id].end();
  }
//COMP(ReadPhredContainer_returnEndIterator);
}

string ReadPhredContainer::returnSuffix(Suffix_t &suf){
//START(ReadPhredContainer_returnSuffix);
  if (suf.type == HEALTHY) {
    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnSuffix() out of bounds " << endl;
      exit(1);
    }
    return HealthyReads[suf.read_id].substr(suf.offset);
  }
  else { // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnSuffix() out of bounds " << endl;
      exit(1);
    }
    return TumourReads[suf.read_id].substr(suf.offset);
  }
//COMP(ReadPhredContainer_returnSuffix);
}

unsigned int ReadPhredContainer::getSize(bool tissueType) {
  if (tissueType == HEALTHY) {
    return HealthyReads.size();
  }
  else {    // == TUMOUR
    return TumourReads.size();
  }
}

string & ReadPhredContainer::getReadByIndex(int index, int tissue) {
//START(ReadPhredContainer_getReadByIndex);
  if(tissue == HEALTHY) {
    if (index >= HealthyReads.size() || index < 0) {
      cout << "getReadByIndex() out of bounds" << endl;
      exit(1);
    }
    return HealthyReads[index];
  }
  else {  // tissue == TUMOUR || tissue == SWITCHED
    if (index >= TumourReads.size() || index < 0) {
      cout << "getReadsByIndex() out of bounds" << endl;
      exit(1);
    }
    return TumourReads[index];
  }
//COMP(ReadPhredContainer_getReadByIndex);
}
string & ReadPhredContainer::getPhredString(int index, int tissue) {
//START(ReadPhredContainer_getPhredString);
  if(tissue == HEALTHY) {
    if (index >= HealthyPhreds.size() || index < 0) {
      cout << "getReadByIndex() out of bounds" << endl;
      exit(1);
    }
    return HealthyPhreds[index];
  }
  else {  // tissue == TUMOUR || tissue == SWITCHED
    if (index >= TumourPhreds.size() || index < 0) {
      cout << "getReadsByIndex() out of bounds" << endl;
      exit(1);
    }
    return TumourPhreds[index];
  }
//COMP(ReadPhredContainer_getPhredString);
}

char ReadPhredContainer::baseQuality(int index, int tissue, int pos) {
  if (tissue == HEALTHY) {
    return HealthyPhreds[index][pos];
  }
  else { // tissue == TUMOUR || tissue == SWITCHED
    return HealthyPhreds[index][pos];
  }
}

int ReadPhredContainer::getMinSuffixSize() {
  return MIN_SUFFIX_SIZE;
}

void ReadPhredContainer::free() {
  HealthyReads.clear();
  HealthyPhreds.clear();
  TumourReads.clear();
  TumourPhreds.clear();
}

// end of file
/*

void ReadPhredContainer::printAllReads() {
  ofstream osock("/data/ic711/checking_read_trim_len.txt");
  for (string s : HealthyReads) {
    s.pop_back();
    osock << s << endl;
  }
  for (string s : TumourReads) {
    s.pop_back();
    osock << s << endl;
  }
  osock.close();
}

void ReadPhredContainer::printRemainingReads(std::string const& filename) {
  ofstream fileHandle(filename.c_str());

  for (unsigned int i=0; i < HealthyReads.size(); i++) {
    fileHandle << "(" << i << ",H)"  << std::endl;
  }
  for (unsigned int i=0; i < TumourReads.size(); i++) {
    fileHandle << "(" << i << ",T)"  << std::endl;
  }
  fileHandle.close();
} 

void ReadPhredContainer::printReadsAndId(int from, int to, int step) const {
  ofstream ofile("/data/ic711/readIdEquivICSmuFin.txt");
  if (from < 0 || to > HealthyReads.size() || to > TumourReads.size()) {
    cout << "Out of range" << endl;
    exit(1);   // should throw
  }
  ofile << "Healthy Reads" << endl;
  for (int it = from; it < to; it += step) {
    ofile << HealthyReads[it] << " : " << it << endl;
  }

  ofile << "Cancer Reads" << endl;
  for (int it = from; it < to; it += step) {
    ofile << TumourReads[it] << " : " << it << endl;
  }

  ofile.close();
}
void ReadPhredContainer::printReads(){
  std::cout << "Healthy file reads: " << endl;
  std::cout << "Size of HealthyReads: " << getSize(HEALTHY) << endl;
  for(string s : HealthyReads) {
    std::cout << s << endl;
  }
  std::cout << endl << endl;

  std::cout << "Tumour file reads: " << endl;
  std::cout << "Size of TumourReads: " << getSize(TUMOUR) << endl;
  for(string s : TumourReads) {
    std::cout << s << endl;
  }
}

void ReadPhredContainer::writeContainer(vector<string> const& c, string const& fname) {
  ofstream of(fname.c_str());
  for (string const& s : c) {
    of << s << endl;
  }
  of.close();
}

string ReadPhredContainer::performDistalTrim(string & s) {
  s.erase(0, distal_trim_len);
  s.erase(s.length() - distal_trim_len);
  return s;
}

*/
