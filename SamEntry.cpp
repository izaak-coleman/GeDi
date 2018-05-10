/*
SamEntry.cpp
Author: Izaak Coleman
*/


#include <string>
#include <vector>
#include "util_funcs.h"
#include "SamEntry.h"
#include "benchmark.h"

using namespace std;


// //COMPULSORY SAM FIELDS
const int SamEntry::TCNS  = 0;
const int SamEntry::FLAG  = 1;
const int SamEntry::RNAME = 2;
const int SamEntry::POS   = 3;
const int SamEntry::MAPQ  = 4;
const int SamEntry::CIGAR = 5;
const int SamEntry::RNEXT = 6;
const int SamEntry::PNEXT = 7;
const int SamEntry::TLEN  = 8;
const int SamEntry::SEQ   = 9;
const int SamEntry::QUAL  = 10;

// GeDi FIELDS
const int SamEntry::LEFT_OHANG = 11;
const int SamEntry::RIGHT_OHANG = 12;

SamEntry::SamEntry(string const& entry, vector<string> const & tumour_cns) {
  vector<string> fields;
  split_string(entry, "\t", fields);
  del = false;
  // ints
  flag = stoi(fields[SamEntry::FLAG]); 
  pos  = stoi(fields[SamEntry::POS]); 
  mapq = stoi(fields[SamEntry::MAPQ]); 
  pnext= stoi(fields[SamEntry::PNEXT]); 
  tlen = stoi(fields[SamEntry::TLEN]); 

  // Parse header
  vector<string> header_subs;
  split_string(fields[0], ";", header_subs);
  index = stoi(header_subs[0]);
  left_ohang = stoi(header_subs[1]);
  right_ohang = stoi(header_subs[2]);

  // strings
  t_cns = tumour_cns[index];
  rname = fields[SamEntry::RNAME]; 
  cigar = fields[SamEntry::CIGAR]; 
  rnext = fields[SamEntry::RNEXT]; 
  seq   = fields[SamEntry::SEQ]; 
  qual  = fields[SamEntry::QUAL]; 

}


bool SamEntry::containsIndel() {
  if (cigar.find('i') || cigar.find('d')) {
    return true;
  }
  return false;
}

void SamEntry::free() {
  del = true;
  SNVLocations.clear();
  t_cns.clear();
  rname.clear();
  cigar.clear();
  rnext.clear();
  seq.clear();
  qual.clear();
}

bool SamEntry::deleted() {
  return del;
}

bool SamEntry::indel() {
  return (cigar.find('D') != string::npos) || (cigar.find('I') != string::npos);
}

