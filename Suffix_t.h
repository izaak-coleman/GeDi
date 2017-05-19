// Suffix_t.h
#ifndef SUFFIX_T_H
#define SUFFIX_T_H

struct Suffix_t {

  /* The three variables, read_id, offset and type are required to
   * co-ordinate the Suffix_t object, to the exact location of the suffix
   * it points to in the Reads object. To achive this, the three
   * variables are used as follows:
   * - type: Associates the Suffix_t object with one of the two vectors, 
             either Reads.TumourReads or Reads.HealthyReads
   * - read_id: Associates Suffix_t object with the read at index read_id, 
                of the vector specified by type. 
   * - offset: Associates the Suffix_t object with the suffix, of the read
               pointed to by read_id, at the vector pointed to by type, such 
               that, the suffix is R[offset..m]$ */

  unsigned int read_id;
  uint16_t offset;
  bool type;
};
#endif

