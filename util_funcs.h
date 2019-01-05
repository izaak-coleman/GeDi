/*
util_funcs.h
Author: Izaak Coleman
*/

#ifndef UTIL_FUNCS_H 
#define UTIL_FUNCS_H

#include <vector>
#include <string>

enum {TUMOUR, HEALTHY, SWITCHED}; // Data type
enum {LEFT, RIGHT}; // Orientation in memory

std::string get_dir_from_filename(std::string const & fname);

void split_string(std::string const & s, std::string const & tokens, std::vector<std::string>
    & fields);
// Splits string at token into fields

std::string reverseComplementString(std::string const& s);

struct consensus_pair {
  std::string mutated;      // Tumour derived consensus sequence 
  std::string non_mutated;  // Healthy derived consensus sequence
  std::string mqual;        // Tumour derived quality string
  std::string nqual;        // Healthy derived quality string
  int mut_offset;           // Offsets for alignments
  int nmut_offset;
  int left_ohang;           // 0 index calibration after trimming
  int right_ohang;
};

#endif
