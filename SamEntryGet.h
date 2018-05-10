/* 
SamEntryGet.h
Author: Izaak Coleman
*/

#include "SamEntry.h"

template <typename T>
T get(int, SamEntry const *){
  // Dummy function. This function will (and should correctly) cause
  // compile time error if instantiated. Functions purpose is to
  // enable template specialization by return type of get()
  const int k = 0; k = 1;
}
// int spec.
template <>
int get<int>(int k, SamEntry const * s) {
  switch (k) {
    case 1:  return s->flag;
    case 3:  return s->pos;
    case 4:  return s->mapq;
    case 7:  return s->pnext;
    case 8:  return s->tlen;
    case 11: return s->left_ohang;
    case 12: return s->right_ohang;
  }
}
// string spec.
template <>
std::string get<std::string>(int k, SamEntry const * s) {
  switch (k) {
    case 0:  return s->t_cns;
    case 2:  return s->rname;
    case 5:  return s->cigar;
    case 6:  return s->rnext;
    case 9:  return s->seq;
    case 10: return s->qual;
  }   
}

