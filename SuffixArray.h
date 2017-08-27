/*

SuffixArray.h
Author: Izaak Coleman
*/


#ifndef SUFFIXARRAY_H
#define SUFFIXARRAY_H

#include <iostream>
#include <vector>
#include <string>

#include "Suffix_t.h"
#include "Reads.h"

class SuffixArray {
private:
  const int N_THREADS;
  const int MIN_SUFFIX_SIZE;
  ReadPhredContainer *reads;
  std::vector<Suffix_t> SA;

  void constructColouredGSA(int min_suffix);
  // Primary function call to generate GeDi's coloured Generalized Suffix Array.

  void constructSuffixArray(int64_t * &dSA, 
                             int64_t *startOfTumour, 
                             int64_t *sizeOfRadixSA);
  // Constructs a suffix array via RadixSA (Rajasekaran, Nicolae 2014).
  // The input for RadixSA is a concatenation of all the reads in the
  // data set. The first half of the concatenation consists of reads
  // from the healthy data set, the second, reads from the cancer data set.
  // The position in the concatenation where cancer reads begin is stored
  // in startOfTumour.

  void buildBinarySearchArrays(
      std::vector<std::pair<unsigned int, int64_t> > *healthyBSA, 
      std::vector<std::pair<unsigned int, int64_t> > *tumourBSA);
  // Wrapper calling generateBSA() for tumour and healthy data sets
  void generateBSA(std::vector<std::pair<unsigned int, int64_t> > &BSA,
                   bool type);
  // Produces a sorted vector that is used for binary searches (binary
  // search array). Elements are (r,c) tuples, of which there is one
  // element per read in the binary search array, where:
  // -- r = index of read in its container (either HealthyReads or TumourReads).
  // -- c = index of read in the concatenation
  // Rank of the tuples implies same rank for tuple fields: For all tuples in 
  // the binary search array, if (r1,c1) < (r2,c2) then r1 < r2 and c1 < c2

  void transformSuffixArrayBlock(std::vector<Suffix_t> *block, 
      std::vector<std::pair<unsigned int, int64_t> > *healthyBSA,
      std::vector<std::pair<unsigned int, int64_t> > *tumourBSA,
      int64_t *SA, int64_t from, 
      int64_t to, int64_t startOfTumour, int min_suf);
  // Function transforms a suffix array element, an integer,
  // into a coloured GSA element, a tuple of type Suffix_t.
  // For a given suffix array element, the fields of the Suffix_t it is 
  // transformed into are given by:
  // -- Suffix_t::type: If the suffix array element is less or greater than
  //    StartOfTumour
  // -- Suffix_t::read_id = r from the (r,c) returned by a
  //    SuffixArray:binarySearch() for the suffix array element
  // -- Suffix_t::offset = (suffix array element - c) for the same (r,c) used
  //    by Suffix_t::read_id. 

  std::string concatenateReads(bool type);
  // Returns a concatenation of the reads in either HealthyReads or TumourReads.

public:
  SuffixArray(ReadPhredContainer &reads, int min_suffix, int n_threads);

  Suffix_t & getElem(unsigned int index);

  unsigned int getSize();

  std::pair<unsigned int, int64_t> 
  binarySearch(std::vector<std::pair<unsigned int,  int64_t> > &BSA, 
               unsigned int suffix_index);
  // Searches over a binary search array returning the (r,c) where
  // c is less than suffix_index and there is no other (r,c) with a
  // smaller for (suffix_index - c)

  void free();
  // Deallocated memory allocated for SA
};
/*
  void buildGSAFile(std::vector<Suffix_t> &GSA, std::string filename);

  void constructGSAFromFile(std::vector<Suffix_t> &GSA, std::string filename);

  void constructTotalRadixSA(uint8_t min_suffix);
  // this function uses radixSA from (Sanguthevar Rajasekaran, Marius Nicolae, 
  // 2014). The function contruct in parallel two suffix arrays, 
  // one for each tissue type using radixSA, it then computes a genealizes
  // suffix array using Suffix_t types, then finally
  // merges the two suffix arraus

  void generalizedRadixSA(std::vector<Suffix_t> *TissueSA, 
      bool type, uint8_t min_suffix);
  // Function generates tissue specific suffix array using 
  // radixSA and wrapper algorithms to transform to generalized suffix array


  // MERGE SORT FUNCTIONS
  
  bool lexCompare(Suffix_t &lhs, Suffix_t &rhs);
  // Function lexiographically compares two suffixes, returning true if 
  // lhs is before rhs, false otherwise. 


  void loadUnsortedSuffixes(uint8_t min_suffix);
  // Function reads through Text (DNA reads) generating a suffix_t for 
  // each suffix of each read, and loads into SA

  void lexMergeSort();
  // Functin lexicographically sorts the suffixes in SA, by the end of this
  // function SA has been transformed into an actual suffix Array

  void sort(unsigned int from, unsigned int to, unsigned int level);
  // recursive divide function or lexMergeSort
  
  void merge(unsigned int from, unsigned int mid, unsigned int to);
  // conquer function of lexMergeSort
  
  std::vector<Suffix_t> *copyOf(unsigned int from, unsigned int to);
  // SA copy function used by merge. Returns a pointer to a subsection
  // of SA spanning from indicies [from, to). On the heap

  void printSuffixArray(std::string const& filename);
  // print the suffixes of the suffix array to file
  void printSuffixes();
  // Function prints out the suffix strings assoc. with each element
  // in the suffix array

  void printReadsInGSA(std::string const& filename);

  // Outputs (idx,type) tuples

  void printSuffixData();
*/


#endif 

