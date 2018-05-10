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
  int index;

  // strings 
  std::string t_cns;
  std::string rname;
  std::string cigar;
  std::string rnext;
  std::string seq;  // which is the healthy cns
  std::string qual; 

  /* For the following three vectors SNVLocations[i] specifies an SNV at
   * genome location pos + SNVLocations[i], for which the healthy dataset
   * allele is nonMutatedBases[i], and the tumour dataset allele is
   * mutatedBases[i].
   */

  // SNV Positions relative to healthy consensus sequence
  std::vector<int> SNVLocations;
  // List of mutated SNV allele
  std::vector<char> mutatedBases;
  // LIst of non-mutated SNV allele
  std::vector<char> nonMutatedBases;


  // delete switch
  bool del;

  // keys
  static const int TCNS;  
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


  SamEntry(std::string const& line, std::vector<std::string> const & tumour_cns); 
  // Separates fields, discards any non compulsory fields. 

  bool containsIndel(); 
  // returns true if CIGAR string contains 'I' or 'D', false otherwise.

  void free();
  // frees memory allocated to string fields and SNVLocations

  bool deleted();
  // Returns true if free() was called on object, false otherwise

  bool indel();
  // Returns true of alignment contains indel relative to ref.

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
      case 0:  t_cns = v; break;
      case 2:  rname = v; break;
      case 5:  cigar = v; break;
      case 6:  rnext = v; break;
      case 9:  seq   = v; break;
      case 10: qual  = v; break;
    }
  }
};
#endif
