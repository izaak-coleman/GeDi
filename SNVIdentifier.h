/*
SNVIdentifier.h
Author: Izaak Coleman
*/


#ifndef SNVIDENTIFIER_H 
#define SNVIDENTIFIER_H 

#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <utility>
#include <mutex>

#include "util_funcs.h"
#include "Suffix_t.h"
#include "SuffixArray.h"

// read_tag is used for elements of second GSA
// to determine which cancer specific reads
// derive from the same genomic locations.
struct read_tag{
  unsigned int read_id;
  mutable int offset;
  mutable bool orientation;
  int tissue_type;
};

struct read_tag_compare{
  // Functor for read_tag comparison logic.
  bool operator() (const read_tag &a, const read_tag &b) const {
    // Within the set read_tags are:
    // -- First sorted by tissue type TUMOUR = SWITCHED < HEALTHY
    // -- Then sorted by read_id
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

class SNVIdentifier {
private:
  const char MIN_PHRED_QUAL;
  // Read characters with a phred score < MIN_PHRED_QUAL will not contribute to 
  // the consensus sequence.

  const int GSA1_MCT;
  // Blocks of suffixes covering the same genomic location in the coloured GSA,
  // that are mostly cancer reads with a cancer read coverage of < GSA1_MCT
  // are not extracted.

  int GSA2_MCT;
  // Blocks of suffixes covering the same genomic location in the second GSA
  // that have a read coverage of < GSA2_MCT are not extracted.
  // Consensus sequence positions generated from < GSA2_MCT reads are
  // masked. Tumour consensus sequence positions with < GSA2_MCT reads
  // supporting the chosen consensus base are masked.

  const int COVERAGE_UPPER_THRESHOLD;
  // Blocks with a coverage > COVERAGE_UPPER_THRESHOLD are discarded.
  const int N_THREADS;
  const int MAX_LOW_CONFIDENCE_POS;
  // Blocks with > MAX_LOW_CONFIDENCE_POS number of low confidence postions
  // are discarded.
  const double ECONT;
  const double ALLELIC_FREQ_OF_ERROR;
  // Aligned positions with > 1 base with a frequency above
  // ALLELIC_FREQ_OF_ERROR are considered low confidence positions and masked.

  std::mutex cancer_extraction_lock;
  std::mutex cout_lock;
  std::mutex extractGroupsWorkerLock;
  std::mutex buildConsensusPairLock;

  ReadPhredContainer *reads;
  SuffixArray *SA; 
  std::set<unsigned int> CancerExtraction;
  // Elements are indicies of cancer specific reads extracted from the
  // coloured GSA
  std::vector<bpBlock*> SeedBlocks;
  // Elements are initial break point blocks, generated from 
  // seedBreakPointBlocks(), each block consists of the the cancer specific
  // reads covering single location.

  std::vector<bpBlock> BreakPointBlocks;
  // Fully formed break point blocks. Consisting of a set of cancer
  // specific reads covering a potential mutation, and the corresponding
  // healthy reads that cover the same location, along with the relative alignment
  // information, used to align each read against its own block.

  std::vector<consensus_pair> consensus_pairs;
  // Elements are consensus sequences generated from break point blocks.

  unsigned int backUpSearchStart(unsigned int seed_index);
  // The worker threads that seed break point blocks from the second gsa
  // may begin block generation from the middle of block. 
  // backUpSearchStart() moves the start of block generation to the
  // start of the block a thread may have landed in the middle of.

  void extractCancerSpecificReads();
  // Divides coloured GSA into even chunks and deploys threads
  // with each thread computing extractGroupsWorker() over an allocated chunk.

  void extractionWorker(unsigned int to, unsigned int from);
  // Passes through allocated GSA chunk and extracts reads as cancer specific 
  // reads if:
  // -- Suffixes of cancer specific reads form a block that covers the
  //    same genomic location, with the number of cancer reads >= GSA_MCT1
  // -- The ratio of healthy / cancer reads in the block is <= ECONT

  void seedBreakPointBlocks();
  // Generates the seed break point blocks from the cancer specific reads
  // in CancerExraction, storing result in SeedBlocks.
  // The second GSA is firsly constructed. Again, RadixSA (Rajasekran, Nicolae
  // 2014) is employed to generate a suffix array, which is then split
  // into even chunks, each deployed on a thread and 
  // transformed with transformBlock() into a second GSA chunk.
  // extractGroups() then generates seed blocks from second GSA.

  void transformBlock(unsigned long long *from, unsigned long long *to,
                      std::vector< std::pair<unsigned int, unsigned int> > *bsa,
                      std::vector<read_tag> * block);
  // Transforms a block of the suffix array elements into a block of the
  // second GSA elements.
  
  void extractGroups(std::vector<read_tag> const& gsa);
  // Divides second GSA gsa into even chunks and deploys threads
  // with each thread computing extractGroupsWorker() over an allocated chunk.

  void extractGroupsWorker(unsigned int seed_index, unsigned int to,
                           std::vector<read_tag> const* gsa_ptr);
  // Passes through the allocated second GSA chunk extracting sets
  // of reads (seed break point blocks), represented by read_tag if:
  // -- Suffixes of reads cover same genomic location
  // -- Number of reads in seed break point block is > GSA_MCT2

  void buildConsensusPairsWorker(bpBlock** block, bpBlock** end);
  // Builds consensus pairs out of an allocated group of seed break point
  // blocks. 

  void generateConsensusSequence(bool tissue, bpBlock const& block, int &
      cns_offset, std::string & cns, std::string & qual);
  // For a given break point block, function builds a consensus sequence,
  // and quality string for tumour or healthy derived reads.
  // And returns the alignment position of the consensus sequence (cns_offset)

  void buildQualityString(std::string & qual, std::vector<int> const&
      freq_matrix, std::string const& cns, bool tissue);
  // Generates the consensus sequence quality string with various
  // masking codes

  void extractNonMutatedAlleles(bpBlock &block, consensus_pair &pair);
  // Using the cancer consensus sequence extracts healthy reads 
  // covering the same location and inserts them into the block.
  // This is done by searching for a 30 character sequence of the consensus
  // sequence and extracting healthy reads sharing the sequence
  // If search fails, the two regions flanking of the failed search
  // sequence are searched.

  long long int binarySearch(std::string query);
  // Binary search for query in coloured GSA, returns:
  //  -- failed search: -1 
  //  -- success: index of matching suffix

  int lcp(std::string const& l, std::string const& r, unsigned int mlr);
  // optimized lcp computer. Optimization: Starts search lcp count and comparison
  // mlr and l[mlr],r[mlr] respectively. As characters at index < mlr in
  // suffixes of the GSA will be identical

  bool lexCompare(std::string const& l, std::string const& r, unsigned int min_lr);
  // Performs an optimized comparison, with same optimization as lcp()

  int minVal(int a, int b);

  bool extendBlock(int seed_index, bpBlock 
      &block, bool orientation, int calibration);
  // seed_index, the location where binarSearch() found a match to 
  // the query is used as the start point to gather healthy reads
  // that also match the query towards the left and the right in the
  // coloured GSA.
  // Returns false if search failed to extract any healthy reads,
  // true otherwise
  bool getSuffixesFromLeft(int seed_index, 
                           bpBlock &block,
                           bool orientation, int calibration);
  // called by extendBlock(). Returns false if failed to 
  // extract any healthy reads, true otherwise
  bool getSuffixesFromRight(int seed_index, 
                           bpBlock &block, 
                           bool orientation, int calibration);
  // called by extendBlock(). Returns false if failed to 
  // extract any healthy reads, true otherwise

  bool excessLowQuality(consensus_pair & pair);
  // Returns true if number of low quality positions is >
  // MAX_LOW_CONFIDENCE_POS, false otherwise.

  int computeLCP(read_tag const& a, read_tag const& b);
  // Computes the lcp between two read_tag suffixes, returning the lcp

  int convertOffset(read_tag const& tag);
  // converts alignment offset between forward and reverse orientation,
  // returning the converted offset

  void trimHealthyConsensus(consensus_pair & pair);
  // Trims distal regions of the healthy consensus sequence
  // up to the last consensus character that has a corresponding
  // quality string value of 'T' or 'B', and succeeds a stretch
  // of consensus characters that have a contigous stretch of 'T' or 'B'
  // quality string values, starting from left and right ends.

  void trimCancerConsensus(consensus_pair & pair);
  // Trims the distal regions of the cancer consensus sequence,
  // such that the ends of the cancer consensus sequence never extend
  // beyond the healthy consensus sequence in the aligned pair.

  void maskLowQualityPositions(consensus_pair & pair);
  // In the aligned consensus pair, at a given aligned position,
  // switches the cancer consensus sequence character 
  // with the healthy consensus sequence character,
  // if the character in either one of the quality strings at the aligned
  // position is not '-'. 
  
  std::string addGaps(int ngaps);
  // adds gaps preceeding a read to give aligned output.

  std::string readTagToString(read_tag const& tag);

  void mergeBlocks(bpBlock & to, bpBlock & from);

  void unifyBlocks(std::vector<bpBlock> & seedBlocks);

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
  void free();
  // Deallocates memosy allocated for consensus_pairs

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
