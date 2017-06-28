/*
SamEntry.h
Author: Izaak Coleman
*/


#ifndef SAMENTRY_H
#define SAMENTRY_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <boost/any.hpp>

class SamEntry {
public:
  // SAM FIELDS
  // ints
  int flag;
  int mapq;
  int pos;
  int pnext;
  int tlen;
  int left_ohang;
  int right_ohang;

  // strings 
  std::string hdr;
  std::string rname;
  std::string cigar;
  std::string rnext;
  std::string seq;
  std::string qual;

  // SNV index
  std::vector<int> SNVLocations;

  // delete switch
  bool del;

  // keys
  static const int HDR;  
  static const int FLAG;
  static const int RNAME;
  static const int POS;
  static const int MAPQ;
  static const int CIGAR;
  static const int RNEXT;
  static const int PNEXT;
  static const int TLEN;
  static const int SEQ;
  static const int QUAL;
  static const int LEFT_OHANG;
  static const int RIGHT_OHANG;
  static const int BLOCK_ID;


  SamEntry(std::string const& line); 
  // Separates fields, discards any non compulsory fields. 



  int snvLocSize(); 
  // number of SNVs identified in the sam entry

  int snvLocation(int idx);
  // returns the location of an SNV 

  void snv_push_back(int v);

  bool containsIndel(); 
  // returns true if CIGAR string contains 'I' or 'D', false otherwise.

  void free();
  // frees memory allocated to string fields and SNVLocations

  bool deleted();
  // Returns true if free() was called on object, false otherwise

  void set(int k, int v) {
    switch (k) {
      case 1:  flag        = v; break;
      case 3:  pos         = v; break;
      case 4:  mapq        = v; break;
      case 7:  pnext       = v; break;
      case 8:  tlen        = v; break;
      case 11: left_ohang  = v; break;
      case 12: right_ohang = v; break;
    }

  }
  void set(int k, std::string const& v) {
    switch (k) {
      case 0:  hdr   = v; break;
      case 2:  rname = v; break;
      case 5:  cigar = v; break;
      case 6:  rnext = v; break;
      case 9:  seq   = v; break;
      case 10: qual  = v; break;
    }
  }
};
#endif
