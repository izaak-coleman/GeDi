#ifndef UTIL_FUNCS_H 
#define UTIL_FUNCS_H

#include <vector>
#include <string>
#include <utility>



#include "Suffix_t.h"

class ReadPhredContainer;       // forward decl.

enum {TUMOUR, HEALTHY, SWITCHED};
enum {LEFT, RIGHT};         // distinguish, paired end type

static const int MIN_SUFFIX = 30;

int computeLCP(Suffix_t &isuf, Suffix_t &jsuf, ReadPhredContainer &reads);
// Returns the longest common prefix between isuf and jsuf suffixes
struct mutation_classes{
  std::vector<int> SNV_pos;
};

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
