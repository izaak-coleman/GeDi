#ifndef BRANCHPOINTGROUPS_H

#define BRANCHPOINTGROUPS_H

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



// Moved the break point blocks from vector<set<read_tag, read_tag_compare>>
// to a struct, allowing a block id to carried with each break point block
// this is used for development purposes and will be redundant
// once the algorithm has been developed

struct bp_block {
  std::set<read_tag, read_tag_compare> block;
  unsigned int id;

  // wrap insert func to avoid refactoring
  void insert(read_tag r) {
    block.insert(r);
  }

  int size() {
    return block.size();
  }

  void clear() {
    block.clear();
    id = 0;
  }
};


class BranchPointGroups {
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
  // CancerExtraction contains sets of ints, where each set
  // covers some mutation that was gathered, and the ints represent
  // the read_id that covered the mutation. 
  // CancerExtraction is used as a storage container, that is
  // used to unify reads that cover the same mutation and are
  // in the same orientation. This is done by unifyComplementaryGroups()
  // The unified groups are then stored in ComplementaryUnified

  std::vector<bp_block> SeedBlocks;  // only contain tumour read subblock
  std::vector<bp_block> BreakPointBlocks;
  std::vector<consensus_pair> consensus_pairs;
  

  void makeBreakPointBlocks();

  unsigned int backUpSearchStart(unsigned int seed_index);

  void transformBlock(unsigned long long *from, unsigned long long *to,
                      std::vector< std::pair<unsigned int, unsigned int> > *bsa,
                      std::vector<read_tag> *block);
  

  void extractCancerSpecificReads();
  // Function makes a pass through SA, grouping suffixes that have an LCP

  // of >= 30. With such an LCP, we are confident they are covering the same
  // genomic location. We can then check the econt ratio of this group, 
  // identifying if the group consists of cancer specific reads.



  bool getSuffixesFromLeft(int seed_index, 
                           std::set<read_tag, read_tag_compare> &block,
                           bool orientation, int calibration);
  // Function gathers suffixes from left (towards 0) in the array
  // that share an lcp of >= 30 with the suffix at SA[seed_index]
  
  
  bool getSuffixesFromRight(int seed_index, 
                           std::set<read_tag, read_tag_compare> &block, 
                           bool orientation, int calibration);
  // Function gathers suffixes from right (towards end) in the array
  // that share an lcp of >= 30 with the suffix at SA[seed_index]

  // Function called directly by makeReadGroup if max LCP is between seed
  // and seed + 1 in LCP. Adds all the suffixes indcies 
  // that have lcp == lcp(seed+1, seed)



  void extractNonMutatedAlleles();
  // Re-interrogate the suffix array, extracting healthy reads
  // looping through each block and searching for each 30bp substring
  // of each read

  void extractNonMutatedAlleles(bp_block &block, consensus_pair &pair);

  bool extendBlock(int seed_index, std::set<read_tag, read_tag_compare> 
      &block, bool orientation, int calibration);
  // Once a read covering a mutated allele
  // has been found, extract reads with >= 30bp lcp in common with seed_index

  bool lexCompare(std::string l, std::string r, unsigned int min_lr);
  // perform a lexographical comparison of strings l, r. 
  // However, to avoid redundant comps, comparison starts from 
  // position min_lr

  int lcp(std::string l, std::string r, unsigned int mlr);
  // avoid redund comps with mlr

  int minVal(int a, int b);
  // return smallest of a, b

  long long int binarySearch(std::string query);
  // Function performs a string search for query in SA.
  // The function performs this search with, in practice
  // O(n + log m) comparisons, rather than O(n log m)
  // by avoiding redundant lex comparisons

  long long int backUpToFirstMatch(long long int bs_hit, std::string query);
  // binarySearch() will it a match, but it may not be the first match.
  // Once binarySearch() finds the match, it calls backUpToFirstMatch()
  // which finds the smallest indexed suffix that matches the query

  void extractionWorker(unsigned int to, unsigned int from);

  void invalidatePosition(std::vector< std::vector<int> > &alignment_counter, 
      int pos);
  // Where number of reads is < TRIM_VALUE this sets bases rows at the
  // insufficient position to -1, invalidating the position, as
  // the number of reads must be at least 0.
  void buildConsensusPairsWorker(bp_block* block, bp_block* end);
  void seedBreakPointBlocks();
  // Input: CancerExtraction
  // Output: Loads the groups of cancer specific reads into BreakPointBlocks
  // Details: Uses a GSA in order to group the reads. Groups are formed
  // from reads contiguous in the array that have LCP >= 30
  void extractGroups(std::vector<read_tag> const& gsa);
  // Input: GSA of cancer specific reads
  // Output: BreakPointBlocks with loaded blocks
  // Details: Forms the groups of contiguous reads with LCP >= 30 
  // using seed and extension 

  void extractGroupsWorker(unsigned int seed_index, unsigned int to,
                           std::vector<read_tag> const* gsa_ptr);

  void buildQualityString(std::string & qual, std::vector<std::vector<int> > const&
      freq_matrix, std::string const& cns, bool tissue);
  // Function steps through each position of the string and
  // determines whether a position should be masked if:
  //  -- Cancer: if the number of bases contributing to the 
  //             consensus base is < CTR then mask. Or, if the 
  //             number of bases with frequency above the error threshold
  //             (ALLELIC_ERROR_THRESH) is > 1 then mask
  //  -- Health: if the number of bases with frequency above the error threshold
  //             (ALLELIC_ERROR_THRESH) is > 1 then mask

  int computeLCP(read_tag const& a, read_tag const& b);

  char revCompCharacter(char ch, bool rc);

  void mergeBlocks(bp_block & to, bp_block & from);

  void unifyBlocks(std::vector<bp_block> & seedBlocks);

  int convertOffset(read_tag const& tag);
  // converts offsets from LEFT and RIGHT orientation

  void printAlignedBlock(bp_block block);

  void trimHealthyConsensus(consensus_pair & pair);
  void trimCancerConsensus(consensus_pair & pair);
  void maskLowQualityPositions(consensus_pair & pair, bool & low_quality);

  
public:


  BranchPointGroups(SuffixArray &SA, 
                    ReadPhredContainer &reads, 
                    char min_phred, int gsa1_mct, int gsa2_mct,
                    int coverage_upper_threshold,
                    int n_threads, int max_low_confidence_pos,
                    double econt, double allelic_freq_of_error);
  // Construtor: Uses SA to load data. Constructor called 
  // generateReadGroups, and so, the object is functional after this call

  ~BranchPointGroups();
  // Destructor dealocates BPG

  unsigned int getSize();
  // returns size of BreakPointBlocks

  std::string reverseComplementString(std::string s);
  // returns the reverse complement of a string

  void printBreakPointBlocks();
  // Prints out each of the breakpoint groups, type and read index

  void printAlignedBlocks();



  bool generateConsensusSequence(unsigned int block_id, int &cns_offset, 
      bool tissue_type, unsigned int &pair_id,
      std::string &cns, std::string & qual);
  // This function performs pileup and returns the consensus sequence
  // for either the TUMOUR or HEALTHY sequence from the block indexed at 
  // block_id

  void generateConsensusSequence(bool tissue, bp_block const& block, int &
      cns_offset, unsigned int & pair_id, std::string & cns, std::string & qual);

  std::string addGaps(int ngaps);
  // Function returns a string of lenth ngaps, where gaps are '-'


  void outputExtractedCancerReads(std::string const& filename);
  // Outputs (idx,type) tuples from CancerExtraction

  void outputFromBPB(std::string const& filename);
  std::string readTagToString(read_tag const& tag);

  consensus_pair & getPair(int i);
  int cnsPairSize();
};



//  void printCancerExtractionGroups();
// prints out the cancer extraction groups from the initial extractions

//  Redundant functions 
//
//  void extractMutationSites();
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


#endif
