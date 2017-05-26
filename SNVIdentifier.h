#ifndef SNVIDENTIFIER_H 
#define SNVIDENTIFIER_H 

#include <vector>
#include <string>
#include <iostream>
#include <set>      // contain reads
#include <utility>  // need coordinates to define read index, bool
#include <mutex>

#include "util_funcs.h"
#include "Suffix_t.h"
#include "SuffixArray.h"

// a read_tag is extracted for each read in a break ppoint
// block determining now the read aligns to the other
// reads in its block, and its aligned orientation
struct read_tag{
  unsigned int read_id;
  mutable int offset;
  mutable bool orientation;
  int tissue_type;
};

struct read_tag_compare{
  // this functor only compares the read id from hashtag
  // the other data is considered metadata that I can use later

  bool operator() (const read_tag &a, const read_tag &b) const {

    // the comparison needs to avoid adding duplicate reads.
    // duplicate reads will have the same read_id val, and
    // derive from the same data set, either H,H or T,S.
    // As such, when comparisons between H,T and H,S are made,
    // it is perfectly valid for both elems to have same read_id val
    // this comparison allows for this
    if( (a.tissue_type == TUMOUR || a.tissue_type == SWITCHED) &&
        (b.tissue_type == TUMOUR || b.tissue_type == SWITCHED) ) {
        return a.read_id < b.read_id;
    }
    else if (a.tissue_type == HEALTHY && b.tissue_type == HEALTHY) {
        return a.read_id < b.read_id;
    }
    else {
        return (a.tissue_type % 2) < (b.tissue_type % 2); // mod to keep SWITCHED and TUMOUR together
    }
  }
};

using bpBlock = std::set<read_tag, read_tag_compare>;
//struct bp_block {
//  std::set<read_tag, read_tag_compare> block;
//  unsigned int id;
//
//  // wrap insert func to avoid refactoring
//  void insert(read_tag r) {
//    block.insert(r);
//  }
//
//  int size() {
//    return block.size();
//  }
//
//  void clear() {
//    block.clear();
//    id = 0;
//  }
//};


class SNVIdentifier {
private:
  const char MIN_PHRED_QUAL;
  const int GSA1_MCT;
  int GSA2_MCT;
  const int COVERAGE_UPPER_THRESHOLD;
  const int N_THREADS;
  const int MAX_LOW_CONFIDENCE_POS;
  const double ECONT;
  const double ALLELIC_FREQ_OF_ERROR;

  std::mutex cancer_extraction_lock;
  std::mutex cout_lock;
  std::mutex extractGroupsWorkerLock;
  std::mutex buildConsensusPairLock;

  ReadPhredContainer *reads;
  SuffixArray *SA;    // store a pointer to SA for access
  std::set<unsigned int> CancerExtraction;
  std::vector<bpBlock> SeedBlocks;  // only contain tumour read subblock
  std::vector<bpBlock> BreakPointBlocks;
  std::vector<consensus_pair> consensus_pairs;
  
  unsigned int backUpSearchStart(unsigned int seed_index);

  void transformBlock(unsigned long long *from, unsigned long long *to,
                      std::vector< std::pair<unsigned int, unsigned int> > *bsa,
                      std::vector<read_tag> *block);
  
  void extractCancerSpecificReads();

  bool excessLowQuality(consensus_pair & pair);

  bool groupSuffixes(int direction, int seedIndex, bpBlock &block, 
      bool orienation, int calibration);

  bool getSuffixesFromLeft(int seed_index, 
                           bpBlock &block,
                           bool orientation, int calibration);
  
  bool getSuffixesFromRight(int seed_index, 
                           bpBlock &block, 
                           bool orientation, int calibration);

  void extractNonMutatedAlleles(bpBlock &block, consensus_pair &pair);

  bool extendBlock(int seed_index, bpBlock 
      &block, bool orientation, int calibration);

  bool lexCompare(std::string const& l, std::string const& r, unsigned int min_lr);

  int lcp(std::string const& l, std::string const& r, unsigned int mlr);

  int minVal(int a, int b);

  long long int binarySearch(std::string query); //O(n + nm)

  void extractionWorker(unsigned int to, unsigned int from);

  void buildConsensusPairsWorker(bpBlock* block, bpBlock* end);

  void seedBreakPointBlocks();

  void extractGroups(std::vector<read_tag> const& gsa);

  void extractGroupsWorker(unsigned int seed_index, unsigned int to,
                           std::vector<read_tag> const* gsa_ptr);

  void buildQualityString(std::string & qual, std::vector<int> const&
      freq_matrix, int width, std::string const& cns, bool tissue);

  int computeLCP(read_tag const& a, read_tag const& b);

  void mergeBlocks(bpBlock & to, bpBlock & from);

  void unifyBlocks(std::vector<bpBlock> & seedBlocks);

  int convertOffset(read_tag const& tag);

  void trimHealthyConsensus(consensus_pair & pair);

  void trimCancerConsensus(consensus_pair & pair);

  void maskLowQualityPositions(consensus_pair & pair);

  void generateConsensusSequence(bool tissue, bpBlock const& block, int &
      cns_offset, std::string & cns, std::string & qual);
  
  std::string addGaps(int ngaps);

  std::string readTagToString(read_tag const& tag);

public:
  SNVIdentifier(SuffixArray &SA, 
                    ReadPhredContainer &reads, 
                    char min_phred, int gsa1_mct, int gsa2_mct,
                    int coverage_upper_threshold,
                    int n_threads, int max_low_confidence_pos,
                    double econt, double allelic_freq_of_error);

  unsigned int getSize();

  consensus_pair & getPair(int i);

  int cnsPairSize();

};




#endif
/*
  bool generateConsensusSequence(unsigned int block_id, int &cns_offset, 
      bool tissue_type, unsigned int &pair_id,
      std::string &cns, std::string & qual);
  long long int backUpToFirstMatch(long long int bs_hit, std::string query);
  // binarySearch() will it a match, but it may not be the first match.
  // Once binarySearch() finds the match, it calls backUpToFirstMatch()
  // which finds the smallest indexed suffix that matches the query

  void extractNonMutatedAlleles();
  // Re-interrogate the suffix array, extracting healthy reads
  // looping through each block and searching for each 30bp substring
  // of each read
  ~BreakPointBlocks();
  // Destructor dealocates BPG

  void outputExtractedCancerReads(std::string const& filename);
  // Outputs (idx,type) tuples from CancerExtraction

  void outputFromBPB(std::string const& filename);
  void maskLowQualityPositions(consensus_pair & pair, bool & low_quality);
  char revCompCharacter(char ch, bool rc);

  std::string reverseComplementString(std::string s);
  // returns the reverse complement of a string

  void printAlignedBlocks();
  void printAlignedBlock(bpBlock block);
  void printBreakPointBlocks();
  // Prints out each of the breakpoint groups, type and read index
//  void printCancerExtractionGroups();
// prints out the cancer extraction groups from the initial extractions

//  Redundant functions 
//
//
//  void extractNonMutatedAllelesDeep();
//
//  void unifyComplementaryGroups();
//
//
//  void discardSingleOrientationBlocks();
//  // After breakpoint blocks have been generated remove blocks that just
//  // contain one orientation
//
//  bool thirtyBasePairOverlap(unsigned int lhs, unsigned int rhs);
//  // check if two reads share a 30bp overlap
//
//bool sequenceMatch(std::string right, std::string left);
//
//  void unifyReverseComplementaryReads();
//  // Function performs the final stage to initialze the breakpoint groups.
//  // This function merges groups that cover the same mutation, 
//  // but in different orientations
  void makeBreakPointBlocks();
  void invalidatePosition(std::vector< std::vector<int> > &alignment_counter, 
      int pos);
  // Where number of reads is < TRIM_VALUE this sets bases rows at the
  // insufficient position to -1, invalidating the position, as
  // the number of reads must be at least 0.
 */
