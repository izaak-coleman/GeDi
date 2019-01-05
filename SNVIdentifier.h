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
#include <list>
#include <utility>
#include <memory>

#include "util_funcs.h"
#include "gsa.h"

// read_tag is used for elements of second GSA
// to determine which cancer specific reads
// derive from the same genomic locations.
struct read_tag{
  int64_t read_id;
  mutable int16_t offset;
  mutable bool orientation;
  unsigned char tissue_type;
  read_tag (){};
  read_tag(int64_t r, int16_t o, bool ori, char t):read_id(r), offset(o),
  orientation(ori), tissue_type(t){};
};

struct read_tag_compare{
  // Functor for read_tag comparison logic.
  bool operator() (const read_tag &a, const read_tag &b) const {
    // Within the set read_tags are:
    // -- First sorted by tissue type TUMOUR = SWITCHED < HEALTHY
    // -- Then sorted by read_id
    return a.read_id < b.read_id;
  }
};

using bpBlock = std::set<read_tag, read_tag_compare>;

struct vBlockComp{
  bool operator() (std::shared_ptr<bpBlock const> const a, std::shared_ptr<bpBlock const> const b) {
    if (a->size() < b->size()) return true;
    if (a->size() > b->size()) return false;
    else {
      for(bpBlock::const_iterator a_it = a->cbegin(), b_it = b->cbegin();
          a_it != a->cend(); ++a_it, ++b_it) {
        if (a_it->read_id == b_it->read_id) continue;
        return a_it->read_id < b_it->read_id;
      }
    }
    return false; // same reads in blocks
  }
};

class SNVIdentifier {
private:
  GSA * gsa;
  const char MIN_PHRED_QUAL;
  // Read characters with a phred score < MIN_PHRED_QUAL will not contribute to 
  // the consensus sequence.

  const int PRI_MSS;
  // Blocks of suffixes covering the same genomic location in the coloured GSA,
  // that are mostly cancer reads with a cancer read coverage of < PRI_MSS
  // are not extracted.

  int AUX_MSS;
  // Blocks of suffixes covering the same genomic location in the second GSA
  // that have a read coverage of < AUX_MSS are not extracted.
  // Consensus sequence positions generated from < AUX_MSS reads are
  // masked. Tumour consensus sequence positions with < AUX_MSS reads
  // supporting the chosen consensus base are masked.

  const int COVERAGE_UPPER_THRESHOLD;
  // Blocks with a coverage > COVERAGE_UPPER_THRESHOLD are discarded.

  int N_THREADS;

  const int MAX_SNPS;
  // Blocks with > MAX_SNPS number of low confidence postions
  // are discarded.

  const double ECONT;

  const double ALLELIC_FREQ_OF_ERROR;
  // Aligned positions with > 1 base with a frequency above
  // ALLELIC_FREQ_OF_ERROR are considered low confidence positions and masked.

  std::set<int64_t> CancerExtraction;
  // Elements are indicies of cancer specific reads extracted from the
  // coloured GSA

  std::vector<std::shared_ptr<bpBlock> > SeedBlocks;
  // Elements are initial break point blocks, generated from 
  // seedBreakPointBlocks(), each block consists of the the cancer specific
  // reads covering single location.


  int64_t backUpSearchStart(int64_t seed_index);
  // The worker threads that seed break point blocks from the second gsa
  // may begin block generation from the middle of block. 
  // backUpSearchStart() moves the start of block generation to the
  // start of the block a thread may have landed in the middle of.

  void extractCancerSpecificReads();
  // Divides coloured GSA into even chunks and deploys threads
  // with each thread computing extractGroupsWorker() over an allocated chunk.

  void extractionWorker(int64_t to, int64_t from, std::set<int64_t> & threadExtr);
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

  void buildConsensusPairs(std::shared_ptr<bpBlock> * block,
      std::shared_ptr<bpBlock> * end,
      std::vector<consensus_pair> & threadWork);
  // Builds consensus pairs out of an allocated group of seed break point
  // blocks. 

  void writeConsensusPairs(std::string const & consensus_pairs, std::string const & fname);

  void generateConsensusSequence(bool tissue, bpBlock const& block,
      int & cns_offset, std::string & cns, std::string & qual);
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

  int64_t binarySearch(std::string const & query);
  // Binary search for query in coloured GSA, returns:
  //  -- failed search: -1 
  //  -- success: index of matching suffix

  int lcp(std::string const & l, std::string const & r, unsigned int mlr);
  // optimized lcp computer. Optimization: Starts search lcp count and comparison
  // mlr and l[mlr],r[mlr] respectively. As characters at index < mlr in
  // suffixes of the GSA will be identical

  bool lexCompare(std::string const & l, std::string const & r, unsigned int min_lr);
  // Performs an optimized comparison, with same optimization as lcp()

  int minVal(int a, int b);

  bool extendBlock(int64_t seed_index, bpBlock 
      &block, bool orientation, int calibration);
  // seed_index, the location where binarSearch() found a match to 
  // the query is used as the start point to gather healthy reads
  // that also match the query towards the left and the right in the
  // coloured GSA.
  // Returns false if search failed to extract any healthy reads,
  // true otherwise
  bool getSuffixesFromLeft(int64_t seed_index, 
                           bpBlock &block,
                           bool orientation, int calibration);
  // called by extendBlock(). Returns false if failed to 
  // extract any healthy reads, true otherwise
  bool getSuffixesFromRight(int64_t seed_index, 
                           bpBlock &block, 
                           bool orientation, int calibration);
  // called by extendBlock(). Returns false if failed to 
  // extract any healthy reads, true otherwise

  bool excessLowQuality(consensus_pair & pair);
  // Returns true if number of low quality positions is >
  // MAX_SNPS, false otherwise.

  int convertOffset(read_tag const& tag);
  // converts alignment offset between forward and reverse orientation,
  // returning the converted offset

  char rc(char c, int d);

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

  bool noSNV(consensus_pair const & pair);
  // Returns true if pair.mutated and pair.non_mutated do not contain
  // a mismatch representative of an SNV. 

  int64_t computeLCP(std::string::const_iterator a, std::string::const_iterator b);

  void buildVariantBlocks(int64_t const * dSA, int64_t const dSA_sz, int64_t * seed_idx, int64_t const
      * const to, std::vector< std::pair<int64_t, int64_t> > const & bsa, 
      std::string const & concat, std::set<std::shared_ptr<bpBlock>, vBlockComp > &
      twork);

  read_tag constructReadTag(int64_t const * dSA, int64_t const dSA_sz, 
      std::vector<std::pair<int64_t, int64_t> > const & bsa, 
      std::string const & concat, int64_t const * ptr);

  void xorSwap(int64_t *x, int64_t *y);

  int64_t bubbleRemove(int64_t * const a, int64_t const sz, int64_t const
      invalid);

  void remove_short_suffixes(int64_t* &sa, int64_t &sa_sz, int64_t
    min_suffix_length, std::string const & concat);

  int64_t suffixLen(int64_t const i, std::string const & concat);

  int assignBaseDisregardingPhred(int const pos, std::vector<read_tag> const & block,
      int const max_offset);


  std::pair<int64_t, int64_t> binarySearch(std::vector< std::pair<int64_t,
      int64_t> > const & bsa, int64_t sa_pos);

  void buildFastqAndTumourCNSData(std::list<consensus_pair> & consensus_pairs, std::string & fastq);

  bool singleIndel(std::string const & cigar);

  bool containsIndel(std::string const & tumour, std::string const & control);

public:
  SNVIdentifier(GSA & gsa,
                    std::string const& basename, 
                    char min_phred, int gsa1_mct, int gsa2_mct,
                    int coverage_upper_threshold,
                    int n_threads, int max_low_confidence_pos,
                    double econt, double allelic_freq_of_error);

  std::vector<std::string> tumour_cns;

  // DEBUG
  void printAlignedBlock(bpBlock block);
  void printExtractedCancerReads();
  std::string addGaps(int ngaps);
};
#endif
