#ifndef UTIL_FUNCS_H 
#define UTIL_FUNCS_H

#include <vector>
#include <string>
#include <utility>

#include "Suffix_t.h"

class ReadPhredContainer; // Forward decl. for use of returnStart/EndIterator()
enum {TUMOUR, HEALTHY, SWITCHED}; // Data type
enum {LEFT, RIGHT}; // Orientation in memory
static const int MIN_SUFFIX = 30;

int computeLCP(Suffix_t &isuf, Suffix_t &jsuf, ReadPhredContainer &reads);

struct mutation_classes{
  std::vector<int> SNV_pos;
};

void split_string(std::string s, std::string tokens, std::vector<std::string>
    &fields);

std::string reverseComplementString(std::string const& s);

struct consensus_pair {
  std::string mutated;
  std::string non_mutated;
  std::string mqual;
  std::string nqual;
  unsigned int pair_id;
  int mut_offset;
  int nmut_offset;
  int left_ohang;
  int right_ohang;
  mutation_classes mutations;
};

#endif
