#ifndef SAMENTRY_H
#define SAMENTRY_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <boost/any.hpp>

class SamEntry {
public:
  // Sam compulsory fields
  static const int HDR;  // <int, string> K, V
  static const int FLAG;  // <int, int> K, V
  static const int RNAME;  // <int, string> K, V
  static const int POS;  // <int, int> K, V
  static const int MAPQ;  // <int, int> K, V
  static const int CIGAR;  // <int, string> K, V
  static const int RNEXT;  // <int, string> K, V
  static const int PNEXT;  // <int, string> K, V
  static const int TLEN;  // <int, string> K, V
  static const int SEQ;  // <int, string> K, V
  static const int QUAL; // <int, string> K, V

  // Sam optional fields
  static const int NM;    // <int, string> K, V 
  static const int MD;    // <int, string> K, V 
  static const int AS;    // <int, string> K, V 
  static const int BC;    // <int, string> K, V 
  static const int X0;    // <int, string> K, V 
  static const int X1;    // <int, string> K, V 
  static const int XN;    // <int, string> K, V 
  static const int XM;    // <int, string> K, V 
  static const int XO;    // <int, string> K, V 
  static const int XG;    // <int, string> K, V 
  static const int XT;    // <int, string> K, V 
  static const int XA;    // <int, string> K, V 
  static const int XS;    // <int, string> K, V 
  static const int XF;    // <int, string> K, V 
  static const int XE;    // <int, string> K, V 

  // GeDi-specific fields
  static const int LEFT_OHANG;   // <int, int> K,V
  static const int RIGHT_OHANG;  // <int, int> K, V
  static const int BLOCK_ID;     // <int, int> K, V


  SamEntry(std::string const& line); 
  template <typename RT>
    // RETURN BY CONST REF
  RT get(int key) {
    try {
      return boost::any_cast<RT> (fields[key]);
    }
    catch(...) {
      std::cout << "exception occured" << std::endl
                << "likely a cast to incorrect type, or\n"
                << "you tried to access an absent key." << std::endl;
    }
  }
  template <typename T>
  void set(int key, T value) {
    try {
      fields[key] = value;
    }
    catch (...) {
      std::cout << "exception occured\n"
                << "likely tried to access an abscent key." << std::endl;
    }
  }

  int snvLocSize();

  int snvLocation(int idx);

  void snv_push_back(int v);

  bool containsIndel(); 
private:
  std::vector<int> SNVLocations;
  unsigned int pair_id;
  std::map<int, boost::any> fields;

  std::string startsWith(std::string const& tok, std::vector<std::string> const&
      fields);
};
#endif
/*
  void setSNVLocation(int idx, int val);
*/
