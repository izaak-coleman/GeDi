// util_funcs.cpp
#include <string>

#include "Suffix_t.h"
#include "util_funcs.h"
#include "Reads.h"

using namespace std;

int computeLCP(Suffix_t &isuf, Suffix_t &jsuf, ReadPhredContainer &reads) {

  // Get suffix pointers in reads
  string::iterator isuf_iter  = reads.returnStartIterator(isuf);
  string::iterator isuf_end   = reads.returnEndIterator(isuf);
  string::iterator jsuf_iter  = reads.returnStartIterator(jsuf);
  string::iterator jsuf_end   = reads.returnEndIterator(jsuf);

  // computes lcp
  int lcp = 0;

  while (isuf_iter != isuf_end &&
         jsuf_iter != jsuf_end && 
         *isuf_iter == *jsuf_iter) {

    lcp++;      // matched next char, so extend lcp

    isuf_iter++; jsuf_iter++;
  }

  return (lcp);   
}
