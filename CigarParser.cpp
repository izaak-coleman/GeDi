// CigarParser.cpp
// Author: Izaak Coleman
#include <string>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <cstdlib>

#include "CigarParser.h"
using namespace std;

static const string CIGAR_OPS = "MIDHNPSX=";
static const string VALID_CHARS = "MIDHNPSX=0123456789";
static const string DIGITS = "0123456789";

CigarParser::CigarParser(string const & cig): CIGAR(cig) {
  assert(CIGAR.find_first_not_of(VALID_CHARS) == string::npos 
         && DIGITS.find(CIGAR[0]) != string::npos
         && CIGAR_OPS.find(CIGAR[CIGAR.size()-1]) != string::npos
         && CIGAR != ""); // Validity check

  sz = std::count_if(CIGAR.begin(), CIGAR.end(), 
          [] (char c) {
            return ((CIGAR_OPS.find(c) != string::npos) ? true : false);
          }
        );
}

CigarParser::CigarParser(string const & cig, bool const reverse): CIGAR(cig) {
  assert(CIGAR.find_first_not_of(VALID_CHARS) == string::npos 
         && DIGITS.find(CIGAR[0]) != string::npos
         && CIGAR_OPS.find(CIGAR[CIGAR.size()-1]) != string::npos
         && CIGAR != ""); // Validity check.

  sz = std::count_if(CIGAR.begin(), CIGAR.end(), 
          [] (char const & c){
            return ((CIGAR_OPS.find(c) != string::npos) ? true : false);
          }
        );
  if (reverse) {
    string rev = "";
    for (int i = sz - 1; i >= 0; --i) {
      rev += std::to_string(this->length_at(i)) + this->operation_at(i);
    }
    CIGAR = rev;
  }
}

int64_t CigarParser::length_at(int64_t const i) const {
  if (i >= sz || i < 0) return -1;
  int64_t last{0}, next{0}, cig_op{0};
  next = CIGAR.find_first_of(CIGAR_OPS, last);
  while (next != string::npos && cig_op < i) {
    last = next+1;
    next = CIGAR.find_first_of(CIGAR_OPS, last);
    ++cig_op;
  }
  try {
    return std::stoi(CIGAR.substr(last, next - last));
  }
  catch(std::invalid_argument) {
    cout << "Could not convert CIGAR operation's length field into integer. "
         << "It is therefore likely the input CIGAR (" << CIGAR
         << ") has invalid syntax. "
         << '\n' << "Program terminating with exit(-1).";
    exit(-1);
  }
}

char CigarParser::operation_at(int64_t const i) const {
  if (i >= sz || i < 0) return -1;
  int64_t last{0}, next{0}, cig_op{0};
  next = CIGAR.find_first_of(CIGAR_OPS, last);
  while (next != string::npos && cig_op < i) {
    last = next+1;
    next = CIGAR.find_first_of(CIGAR_OPS, last);
    ++cig_op;
  }
  return CIGAR[next];
}

int64_t CigarParser::size() const {
  return sz;
}
