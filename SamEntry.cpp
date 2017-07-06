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
const int SamEntry::HDR   = 0;
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
const int SamEntry::ID = 13;

SamEntry::SamEntry(string const& entry) {
  vector<string> fields;
  split_string(entry, "\t", fields);
  del = false;
  // ints
  flag = stoi(fields[SamEntry::FLAG]); 
  pos  = stoi(fields[SamEntry::POS]); 
  mapq = stoi(fields[SamEntry::MAPQ]); 
  pnext= stoi(fields[SamEntry::PNEXT]); 
  tlen = stoi(fields[SamEntry::TLEN]); 

  // strings
  rname = fields[SamEntry::RNAME]; 
  cigar = fields[SamEntry::CIGAR]; 
  rnext = fields[SamEntry::RNEXT]; 
  seq   = fields[SamEntry::SEQ]; 
  qual  = fields[SamEntry::QUAL]; 

  // Parse header
  vector<string> header_subs;
  // split off the cancer cns
  split_string(fields[SamEntry::HDR], "[", header_subs);
  hdr = header_subs[0];
  // cut of the end "]" char
  fields[SamEntry::HDR] = header_subs[1].substr(0, header_subs[1].length()-1);
  // split up the remaining fields
  header_subs.clear();
  split_string(fields[SamEntry::HDR], ";", header_subs);
  left_ohang = stoi(header_subs[0]);
  right_ohang = stoi(header_subs[1]);
  id = stoi(header_subs[2]);
}

void SamEntry::snv_push_back(int v) {
  SNVLocations.push_back(v);
}

int SamEntry::snvLocSize() {return SNVLocations.size();}
int SamEntry::snvLocation(int idx) {return SNVLocations[idx];}
bool SamEntry::containsIndel() {
  if (cigar.find('i') || cigar.find('d')) {
    return true;
  }
  return false;
}

void SamEntry::free() {
  del = true;
  SNVLocations.clear();
  hdr.clear();
  rname.clear();
  cigar.clear();
  rnext.clear();
  seq.clear();
  qual.clear();
}

bool SamEntry::deleted() {
  return del;
}

