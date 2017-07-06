/*
util_funcs.h
Author: Izaak Coleman
*/


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
// Compute longest common prefic between suffixes
void split_string(std::string s, std::string tokens, std::vector<std::string>
    &fields);
// Splits string at token into fields

std::string reverseComplementString(std::string const& s);
std::string* reverseComplementStringHeap(std::string const& s);
// Returns pointer to reverse complement of s stored on the heap.

struct consensus_pair {
  std::string mutated;      // Tumour derived consensus sequence 
  std::string non_mutated;  // Healthy derived consensus sequence
  std::string mqual;        // Tumour derived quality string
  std::string nqual;        // Healthy derived quality string
  int mut_offset;           // Offsets for alignments
  int nmut_offset;
  int left_ohang;           // 0 index calibration after trimming
  int right_ohang;
  int id;
};

#endif
