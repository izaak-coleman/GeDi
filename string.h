//string.h contains useful string functions beyond or extending
// c++ <string> and c <cstring> (<string.h>). 
// Functions often provide c++ interfaces for useful c functions

#ifndef STRING_H
#define STRING_H

#include <cstring>
#include <string>
#include <vector>

void split_string(std::string s, std::string tokens, 
                  std::vector<std::string> &split_strings);

std::string reverseComplementString(std::string s);

#endif
