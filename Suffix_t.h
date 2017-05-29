/*
Suffix_t.h
Author: Izaak Coleman
*/


#ifndef SUFFIX_T_H
#define SUFFIX_T_H

struct Suffix_t {
  // Coloured GSA element.
  unsigned int read_id; // Index of read in container
  uint16_t offset;      // Index of start of suffix in read
  bool type;  // Associates element with either HealthyReads or TumourReads
};
#endif

