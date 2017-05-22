// Reads.cpp
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstring>

// For parallel data processing 
#include <thread>

#include <zlib.h>   // gunzip parser
#include "kseq.h"   // fastq parser

#include "util_funcs.h"
#include "string.h" // split_string()
#include "Reads.h"

KSEQ_INIT(gzFile, gzread);    // initialize .gz parser


using namespace std;
static const int TERM_CHAR_CORRECTION = 1;
static const int MIN_SUFFIX_SIZE = 30;

// QUALITY_THRESH and PHRED_20 are used
// to discard reads where the percentage
// of positions that have a phred score < PHRED_20
// is > QUALITY_THRESH
static const double QUALITY_THRESH = 0.1;   // 10% 
static const char PHRED_20 = '5';           // lowest high quality phred score

static const string REMOVED_TOKENS = "N";   // Remove N input from algo
static const string TERM_CHAR = "$";        // GSA termination char.


ReadPhredContainer::ReadPhredContainer(string const& inputFile){
  minimum_suffix_size = MIN_SUFFIX_SIZE;
  vector<file_and_type> datafiles;
  maxLen = 0;

  parseInputFile(inputFile, datafiles);
  cout << "Loaded " << datafiles.size() << " data files." << endl;
  for(int i=0; i < datafiles.size(); i++) {
    cout << "Extracting data from " << datafiles[i].first << "..." << endl;
    if (datafiles[i].second) { // == HEALTHY
      loadFastqRawDataFromFile(datafiles[i].first, HealthyReads, HealthyPhreds);
    }
    else{ // filetype == TUMOUR
      loadFastqRawDataFromFile(datafiles[i].first, TumourReads, TumourPhreds);
    }
  }
  cout << "MAX READ LEN: " << maxLengthRead() << endl;
}

void ReadPhredContainer::writeContainer(vector<string> const& c, string const& fname) {
  ofstream of(fname.c_str());
  for (string const& s : c) {
    of << s << endl;
  }
  of.close();
}

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

void ReadPhredContainer::parseInputFile(string const& inputFile, 
                                        vector<file_and_type> &datafiles) {
  ifstream sock;
  sock.open(inputFile.c_str());
  cout << "Gathering datafiles from " << inputFile << "." << endl;
  string file_string;
  while (getline(sock, file_string)) {
    file_and_type file_info;
    vector<string> fields;
    split_string(file_string, ",\t ", fields);

    // load data
    file_info.first = fields[0];
    if (fields[1] == "H") file_info.second = HEALTHY;
    else if (fields[1] == "T") file_info.second = TUMOUR;
    else {
      cout << fields[1] 
           << " is not a valid datatype, either H or T." 
           << endl << "Program terminating." << endl;
      exit(1);
    }
    cout << "Input " << file_info.first << " as "
         << ((file_info.second) ? "healthy" : "tumour") << " data "
         << endl;
    datafiles.push_back(file_info);
  }
  sock.close();
}

string ReadPhredContainer::performDistalTrim(string & s) {
  s.erase(0, distal_trim_len);
  s.erase(s.length() - distal_trim_len);
  return s;
}

void ReadPhredContainer::loadFastqRawDataFromFile(string filename, 
                              vector<string> &processed_reads, 
                              vector<string> & processed_phreds) {

  gzFile data_file;
  data_file = gzopen(filename.c_str(), "r");    // open stream to next fastq.gz 
  kseq_t *seq = kseq_init(data_file);           // init parser

  vector<fastq_t> fastq_elements;

  // Load data from file into array
  int eof_check;
  fastq_t next_read;
  while ((eof_check = kseq_read(seq)) >=0) {

    // As we only want high quality reads, 
    // only reads with a quality score have potential to be added

    if (seq->qual.l) {
      next_read.id   = seq->name.s;
      // copy sequence
      next_read.seq  = seq->seq.s;
      // copy quality string
      next_read.qual = seq->qual.s;

      fastq_elements.push_back(next_read);     // load into global vector
    }
  }
  kseq_destroy(seq);
  gzclose(data_file);
  vector<fastq_t> *fastq_elements_p = &fastq_elements;
  qualityProcessRawData(fastq_elements_p, &processed_reads, &processed_phreds, 0, fastq_elements.size(), 0);
}

void ReadPhredContainer::qualityProcessRawData(vector<fastq_t> *r_data, 
                           vector<string> *processed_reads,
                           vector<string> *processed_phreds,
                           int from,
                           int to, int tid){

  processed_reads->reserve(to - from);
  processed_phreds->reserve(to - from);

  // Search through string, first determining quality, then 
  // removing substrings
  double n_low_qual_bases = 0.0;
  for(int i = from; i < to; i++) {
    // link iterators to quality read of the next fastq elem
    string::iterator iter  = (*r_data)[i].qual.begin();
    string::iterator end   = (*r_data)[i].qual.end();
    n_low_qual_bases = 0.0;                           // count low qual
    while (iter != end) {
      if (*iter < PHRED_20) {
        n_low_qual_bases++;
      }
      iter++;
    }
    if((n_low_qual_bases / (*r_data)[i].qual.size()) >  QUALITY_THRESH) {
      continue;   //, as read is low quality.
    }

    // REMEMBER TO REMOVE PERFORM DISTAL TRIM FOR PRODUCTION
    vector<string> read_substrs;  // remove N chars
    vector<string> phred_substrs;
    int left_arrow{0}, right_arrow{0};
    while((right_arrow = (*r_data)[i].seq.find(REMOVED_TOKENS, left_arrow)) != string::npos) {
      if (left_arrow == right_arrow) {
        left_arrow++;
        continue;
      }
      read_substrs.push_back((*r_data)[i].seq.substr(left_arrow, right_arrow - left_arrow));
      phred_substrs.push_back((*r_data)[i].qual.substr(left_arrow, right_arrow - left_arrow));
      left_arrow = right_arrow + 1;
    }
    read_substrs.push_back((*r_data)[i].seq.substr(left_arrow));
    phred_substrs.push_back((*r_data)[i].qual.substr(left_arrow));

    for (int i=0; i < read_substrs.size(); i++) {
      if(read_substrs[i].length() >= MIN_SUFFIX_SIZE /*- TERM_CHAR_CORRECTION */) {
        if (read_substrs[i].length() > maxLen) {   // determine longest read
          maxLen = read_substrs[i].length();
        }
        processed_reads->push_back(read_substrs[i] + TERM_CHAR);
        processed_phreds->push_back(phred_substrs[i]);
      }
    }
  }
}



ReadPhredContainer::~ReadPhredContainer() {
}

unsigned int ReadPhredContainer::maxLengthRead() {
  return maxLen;
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

string::iterator ReadPhredContainer::returnStartIterator(Suffix_t &suf) {
  // Use suf.type and suf.read_id to locate the read, and then set an iterator
  // pointing at suf.offset dist from begining

  string::iterator iter;
  if(suf.type == HEALTHY) { 

    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnStartIterator() out of bounds " << endl;
      exit(1);
    }
    iter = HealthyReads[suf.read_id].begin() + suf.offset;
  }
  else {  // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnStartIterator() out of bounds " << endl;
      exit(1);
    }
    iter = TumourReads[suf.read_id].begin() + suf.offset;
  }

  return iter;
}

string::iterator ReadPhredContainer::returnEndIterator(Suffix_t &suf) {
  // Use suf.type and suf.read_id to locate the read, then return an iterator to the 
  // end of that read

  string::iterator iter;
  if (suf.type == HEALTHY) {
    if (suf.read_id >= HealthyReads.size() || suf.read_id < 0) {
      cout << "returnEndIterator() out of bounds " << endl;
      exit(1);
    }
    iter = HealthyReads[suf.read_id].end();
  }
  else {  // suf.type == TUMOUR
    if (suf.read_id >= TumourReads.size() || suf.read_id < 0) {
      cout << "returnEndIterator() out of bounds " << endl;
      exit(1);
    }
    iter = TumourReads[suf.read_id].end();
  }

  return iter;
}

string ReadPhredContainer::returnSuffix(Suffix_t &suf){
  // return the string assoc. with suf
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
}
string & ReadPhredContainer::getPhredString(int index, int tissue) {
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
  return minimum_suffix_size;
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


// end of file


