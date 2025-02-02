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
#include "CigarParser.h"


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

  void countSNVs(SamEntry* &alignment, int left, CigarParser const & cp);
  // Compares consensus sequences identifying SNVs.

  void parseSamFile(std::vector<SamEntry*> &alignments, std::string filename);
  // Builds list of SamEntries from sam file

  static bool compareSNVLocations(const single_snv &a, const single_snv &b);
  // Comparison logic for single_snv data type.

  static bool equalSNVs(single_snv const& a, single_snv const& b);
  // Returns true of two single_snv's have equivalent fields

  void outputSNVToUser(std::vector<SamEntry*> &alignments, std::string reportFilename);
  // Reports the list of SNVs

  char operationOfSNV(int const SNVpos, CigarParser const & cp);
  // Returns the CIGAR operation in which the SNV at SNVpos is located.

  int calibrateWithCIGAR(int SNVpos, CigarParser const & cp);
  // Calibrates SNVPos according to the insertion/deletion information
  // contained within CIGAR, such that the SNVPos is referencing the 
  // correct reference genome index.

  bool contiguousMismatch(std::string const & tumour, std::string const & control, int ohang);

  

public:
    GenomeMapper(SNVIdentifier &snv,
                 std::string const& basename,
                 std::string const& chr, std::string const& bwt_idx,
                 int min_mapq);
};
#endif
