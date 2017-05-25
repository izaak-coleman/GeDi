#ifndef GENOMEMAPPER_H
#define GENOMEMAPPER_H

#include <string>
#include <vector>
#include <fstream>

#include "BreakPointBlocks.h"
#include "Reads.h"
#include "SamEntry.h"

struct snv_aln_info {
 std::vector<int> SNV_pos;
 int flag;
 int chr;
 int position;
 int left_ohang;
 int right_ohang;
 std::string non_mutated_cns;
 std::string mutated_cns;
 unsigned int pair_id;
};

struct single_snv {
 int flag;
 std::string chr;
 int position;
 char mutation_base;
 char healthy_base;
 unsigned int pair_id;
};

class GenomeMapper {
private:
  const int MIN_MAPQ;
  const std::string CHR;
  BreakPointBlocks *BPB; 
  ReadPhredContainer *reads;

  void identifySNVs(std::vector<SamEntry> &alignments);

  void countSNVs(SamEntry &alignment, int left);

  void constructSNVFastqData(std::string const& samName);

  void parseSamFile(std::vector<SamEntry> &alignments, std::string filename);

  static bool compareSNVLocations(const single_snv &a, const single_snv &b);

  void outputSNVToUser(std::vector<SamEntry> &alignments, std::string reportFilename);


public:
    GenomeMapper(BreakPointBlocks &bpb, ReadPhredContainer &reads,
                 std::string outpath, std::string const& basename,
                 std::string const& chr, std::string const& bwt_idx,
                 int min_mapq);
};

/*
  void callBWA();
  // Function first calls bwa aln:
  // cmd ./bwa/bwa aln -t 16 hg19.fa cns_pairs.fastq > cns_pairs.sai
  // Then cals to bwa samse
  // ./bwa/bwa samse hg19.fa cna_pairs.sai cns_pair.fastq > cns_pairs.sam
  void callBowtie2();
  // Function calls bowtie2 with following command
  // "./bowtie2 -x /data/ic711/insilico_data/bowtie_index/hg19 -U \
  // /data/ic711/result/cns_pairs.fastq -S /data/ic711/result/cns_pairs.sam"
  void buildConsensusPairs();
  // Function fills the consensus_pairs vector with 
  // consensus pairs generated from the breakpoint blocks
  // of BPG
  void maskLowConfidencePositions(consensus_pair &pair, 
      std::vector< std::vector<int> > &healthy_base_freq,
      std::vector< std::vector<int> > &tumour_base_freq,
      bool &discard);
  // Function scans through frequency vectors. If 
  // more than two bases have a frequency over ALLELIC_FREQ_OF_ERROR,
  // then the position is "masked" by writing the 
  // cancer cns position to the healthy cns position

  void maskLowQualityPositions(consensus_pair & pair, bool &low_quality);
  void trimCancerConsensus(consensus_pair &pair);
  // Trims the cancer consensus sequence, so its length is at most
  // the length of the healthy consensus sequence.
    void printConsensusPairs();
    // print out each mutated and non mutated string
    void printMutation(char healthy, char cancer, std::ofstream &mut_file);
    void printAlignmentStructs(std::vector<SamEntry> &alignments);
  void printGaps(int gaps);

  void printAllAlignments(std::vector<SamEntry> &alignments);

  void printSingleAlignment(SamEntry &snv);
  void correctReverseCompSNV(std::vector<SamEntry> &alignments);
  std::string generateParseString(mutation_classes &m);

  // Function generates a string of the following string (mutation string)
  // [SNV:a;b;c][SSV:e;f;g][LSV:x;y;z]
  // this string is stored as the header for each
  // aligned healthy read. This allows identification
  // of the mutation indexes directly from the SAM file
 */
#endif
