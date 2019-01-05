/*
util_funcs.cpp
Author: Izaak Coleman
*/


#include <string>
#include <cstring>
#include <vector>

#include "util_funcs.h"

using namespace std;

string get_dir_from_filename(string const & fname){
  int end_of_path = fname.find_last_of('/');
  if (end_of_path == string::npos) {
    return "./";
  }
  else {
    return fname.substr(0, end_of_path);
  }
}


void split_string(string const & s, string const & tokens, vector<string> & split_strings) {
  char *c_s = const_cast<char*>(s.c_str());
  char *c_tokens = const_cast<char*>(tokens.c_str());
  char *c_split = strtok(c_s, c_tokens); // split string into delimited cstrings
  while (c_split != NULL) {
    split_strings.push_back(c_split);
    c_split = strtok(NULL, c_tokens);
  }
}

string reverseComplementString(string const& s){
  string revcomp = "";
  for(int i = s.size()-1; i >= 0; i--) {
  // travel in reverse and switch for complementary
    switch(s[i]) {
      case 'A': revcomp += "T"; break;
      case 'T': revcomp += "A"; break;
      case 'C': revcomp += "G"; break;
      case 'G': revcomp += "C"; break;
      case 'N': revcomp += "N"; break;
    }
  }
  return revcomp;
}
