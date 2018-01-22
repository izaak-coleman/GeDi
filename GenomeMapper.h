/*
GenomeMapper.h
Author: Izaak Coleman
*/

#ifndef GENOMEMAPPER_H
#define GENOMEMAPPER_H

#include <string>
#include <vector>
#include <fstream>

#include "SNVIdentifier.h"
#include "SamEntry.h"


// In order to sort the mulitple potential SNVs present in a 
// SamEntry, each SNV is stored as a single_snv along with
// relevant information from the SamEntry.
struct single_snv {
 std::string chr;
 int position;
 char mutation_base;
 char healthy_base;
 int pair_id;//  REMOVE!!
};

class GenomeMapper {
private:
  const int MIN_MAPQ;       // Discard SamEntries with MAPQ < MIN_MAPQ
  const std::string CHR;    // Discard SamEntries not aligned to CHR
  SNVIdentifier *snvId;     

  void identifySNVs(std::vector<SamEntry*> &alignments);
  // Calls countSNVs() for all valid SamEntries, translating the
  // consensus sequence complementarity to that of the reference genome.

  void countSNVs(SamEntry* &alignment, int left);
  // Compares consensus sequences identifying SNVs.

  void constructSNVFastqData(std::string const& fastqName);
  // Builds the fastq file from the consensus pairs. Fastq file will
  // be used to align against the reference genome.
  // Fastq element structure:
  //
  // @TumourConsensusSequence[left_ohang;right_ohang]
  // HealthyConsensusSequence
  // +
  // DummyQualityString (!)'s of length HealthyConsensusSequence

  void parseSamFile(std::vector<SamEntry*> &alignments, std::string filename);
  // Builds list of SamEntries from sam file

  static bool compareSNVLocations(const single_snv &a, const single_snv &b);
  // Comparison logic for single_snv data type.

  static bool equalSNVs(single_snv const& a, single_snv const& b);
  // Returns true of two single_snv's have equivalent fields

  void outputSNVToUser(std::vector<SamEntry*> &alignments, std::string reportFilename);
  // Reports the list of SNVs

public:
    GenomeMapper(SNVIdentifier &snv,
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
