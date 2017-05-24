// SuffixArray.h
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
  const int MIN_SUFFIX;
  ReadPhredContainer *reads;
  std::vector<Suffix_t> SA;      // pointer to suffix array
 
  void parallelGenRadixSA(int min_suffix);

  void transformSuffixArrayBlock(std::vector<Suffix_t> *block, 
      std::vector<std::pair<unsigned int, unsigned int> > *healthyBSA,
      std::vector<std::pair<unsigned int, unsigned int> > *tumourBSA,
      unsigned long long *radixSA, unsigned int from, 
      unsigned int to, unsigned int startOfTumour, int min_suf);

  void buildBinarySearchArrays(
      std::vector<std::pair<unsigned int, unsigned int> >
    *healthyBSA, std::vector<std::pair<unsigned int, unsigned int> > *tumourBSA);

  void generateBSA(std::vector<std::pair<unsigned int, unsigned int> > &BSA,
      bool type);

  void generateParallelRadix(unsigned long long **radixSA, 
      unsigned int *startOfTumour, unsigned int *sizeOfRadixSA);

  std::string concatenateReads(bool type);

public:
  SuffixArray(ReadPhredContainer &reads, int min_suffix, int n_threads);

  Suffix_t & getElem(int index);

  unsigned int getSize();

  std::pair<unsigned int, unsigned int> 
    binarySearch(std::vector<std::pair<unsigned int, 
        unsigned int> > &BSA, unsigned int suffix_index);
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

