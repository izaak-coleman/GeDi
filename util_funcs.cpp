// util_funcs.cpp
#include <string>
#include <cstring>
#include <vector>

#include "Suffix_t.h"
#include "util_funcs.h"
#include "Reads.h"
#include "benchmark.h"

using namespace std;

int computeLCP(Suffix_t &isuf, Suffix_t &jsuf, ReadPhredContainer &reads) {
START(Util_computeLCP);
  string::iterator isuf_iter  = reads.returnStartIterator(isuf);
  string::iterator isuf_end   = reads.returnEndIterator(isuf);
  string::iterator jsuf_iter  = reads.returnStartIterator(jsuf);
  string::iterator jsuf_end   = reads.returnEndIterator(jsuf);
  int lcp = 0;
  while (isuf_iter != isuf_end &&
         jsuf_iter != jsuf_end && 
         *isuf_iter == *jsuf_iter) {
    lcp++;      // matched next char, so extend lcp
    isuf_iter++; jsuf_iter++;
  }
  return lcp;   
COMP(Util_computeLCP);
}

void split_string(string s, string tokens, vector<string> &split_strings) {
START(Util_split_string);
  char *c_s = const_cast<char*>(s.c_str());
  char *c_tokens = const_cast<char*>(tokens.c_str());
  char *c_split = strtok(c_s, c_tokens); // split string into delimited cstrings
  while (c_split != NULL) {
    split_strings.push_back(c_split);
    c_split = strtok(NULL, c_tokens);
  }
COMP(Util_split_string);
}

string reverseComplementString(string const& s){
START(Util_reverseComplementString);
  string revcomp = "";
  for(int i = s.size()-1; i >= 0; i--) {
  // travel in reverse and switch for complementary
    switch(s[i]) {
      case 'A': revcomp += "T"; break;
      case 'T': revcomp += "A"; break;
      case 'C': revcomp += "G"; break;
      case 'G': revcomp += "C"; break;
    }
  }
  return revcomp;
COMP(Util_reverseComplementString);
}
