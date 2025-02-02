/*
GenomeMapper.cpp
Author: Izaak Coleman
*/


#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

#include "GenomeMapper.h"
#include "SNVIdentifier.h"
#include "util_funcs.h"
#include "SamEntry.h"
#include "SamEntryGet.h"
#include "CigarParser.h"

#define STR(str) #str
#define STRING(str) STR(str)

#define STR(str) #str
#define STRING(str) STR(str)


using namespace std;

static const int REVERSE_FLAG = 16;
static const int FORWARD_FLAG = 0;

GenomeMapper::GenomeMapper(SNVIdentifier &snv, 
                           string const& basename,
                           string const& chr,
                           string const& bwt_idx,
                           int min_mapq):
                           MIN_MAPQ(min_mapq),
                           CHR(chr) {
  cout << endl << endl << "SNV calling." << endl;
  this->snvId = &snv;
  string fastqName(basename + ".fastq.gz"),
         samName(basename + ".sam"),
         outName(basename + ".SNV_results");

  cout << "Aligning control consensus sequences as proxy..." << endl;
  cout << "Bowtie2 output:" << endl;
  string pwd(STRING(PWD));
  string command_aln(pwd + "/bowtie2-2.3.4/bowtie2 -p 16 -x " + bwt_idx + " -U "
      + fastqName + " -S " + samName);
  system(command_aln.c_str());
  vector<SamEntry*> alignments;

  parseSamFile(alignments, samName);
  cout << "Calculating SNV locations..." << endl;
  identifySNVs(alignments);
  snvId->tumour_cns.clear();
  
  cout << "Calling..." << endl;
  outputSNVToUser(alignments, outName);
}

void GenomeMapper::parseSamFile(vector<SamEntry*> &alignments, string filename) {
  ifstream snvSam(filename);
  string line;
  while(getline(snvSam, line)) {
    if (line[0] == '@') continue;
    SamEntry * entry = new SamEntry(line, snvId->tumour_cns); 
    if (get<string>(SamEntry::RNAME, entry) != CHR && CHR != "") {
      delete entry;
      entry = nullptr;
      continue;
    } 
    if(get<int>(SamEntry::MAPQ, entry) < MIN_MAPQ) { // MAPQ filter
      delete entry;
      entry = nullptr;
      continue;
    }
    alignments.push_back(entry);
  }
  snvSam.close();
}


void GenomeMapper::identifySNVs(vector<SamEntry*> &alignments) {
  for (SamEntry* &entry : alignments) {
    if(get<int>(SamEntry::FLAG, entry) == FORWARD_FLAG) {
      CigarParser cp(get<string>(SamEntry::CIGAR, entry));
      countSNVs(entry, get<int>(SamEntry::LEFT_OHANG, entry), cp);
    }
    else if (get<int>(SamEntry::FLAG, entry) == REVERSE_FLAG) {
      CigarParser cp(get<string>(SamEntry::CIGAR, entry));
      entry->set(SamEntry::TCNS, reverseComplementString(get<string>(SamEntry::TCNS, entry))); 
      countSNVs(entry, get<int>(SamEntry::RIGHT_OHANG, entry), cp); // invert overhangs due to rev comp
    }
  }
}

char GenomeMapper::operationOfSNV(int const SNVPos, CigarParser const & cp) {
  int64_t idx = 0;
  for(int64_t count = 0; count < SNVPos && idx < cp.size(); ++idx) {
    char op = cp.operation_at(idx);
    if (op == 'M' || op == 'I') {
      count += cp.length_at(idx);
    }
  }
  return cp.operation_at(idx-1);
}

int GenomeMapper::calibrateWithCIGAR(int SNVPos, CigarParser const & cp) {
  SNVPos++; // change to base-1 index in line with CIGAR
  if (cp.operation_at(0) == 'M' && cp.length_at(0) >= SNVPos) {
    // Then CIGAR has no effect.
    return SNVPos-1;  // -1 returns base-0 index
  }
  if (operationOfSNV(SNVPos, cp) == 'I') { 
    // Then the SNV is located within an insertion that is not found within
    // the reference genome. Accordingly, no logical position of the SNV
    // exists.
    return -1;
  }

  // Otherwise, calibrate SNVPos according to CIGAR
  int64_t origSNVPos = SNVPos;
  for (int64_t i = 0, count = 0; count < origSNVPos && i < cp.size(); ++i) {
    char op = cp.operation_at(i);
    if (op == 'I') {
      SNVPos -= cp.length_at(i);
    }
    else if (op == 'D') {
      SNVPos += cp.length_at(i);
    }

    // count incrementation ensures operations beyond the original
    // SNV position are not incorrectly calibrated for.
    if (op == 'M' || op == 'I') {
      count += cp.length_at(i);
    }
  }
  return SNVPos-1; // -1 returns base-0 index
}


bool GenomeMapper::contiguousMismatch(string const & mutated, string const & nonMutated, int ohang) {
  bool noSNVs = true;
  for(int i=0; i < mutated.size() - 1; i++) {
    if (mutated[i] != nonMutated[i + ohang] &&
        mutated[i+1] != nonMutated[i+1 + ohang]) {
      return true;
    }
  }
  return false;
}


void GenomeMapper::countSNVs(SamEntry * &alignment, int ohang, CigarParser const
    & cp) {
  string mutated = get<string>(SamEntry::TCNS, alignment);
  string nonMutated = get<string>(SamEntry::SEQ, alignment);
  bool noSNVs = true;
  if (contiguousMismatch(mutated, nonMutated, ohang)) {
    delete alignment;
    alignment = nullptr;
    return;
  }

  // SNV at start
  if (mutated[0] != nonMutated[0 + ohang] &&
      mutated[1] == nonMutated[1 + ohang]) {
      alignment->mutatedBases.push_back(mutated[0]);
      alignment->nonMutatedBases.push_back(nonMutated[0+ohang]);
      alignment->SNVLocations.push_back(calibrateWithCIGAR(0 + ohang, cp)); 
      noSNVs = false;
  }
  // SNV at end 
  int cnsLen = mutated.size();
  if (mutated[cnsLen-1] != nonMutated[cnsLen-1 + ohang] &&
      mutated[cnsLen-2] == nonMutated[cnsLen-2 + ohang]) {
      alignment->mutatedBases.push_back(mutated[cnsLen-1]);
      alignment->nonMutatedBases.push_back(nonMutated[cnsLen-1+ohang]);
      alignment->SNVLocations.push_back(calibrateWithCIGAR(cnsLen-1 + ohang, cp));
      noSNVs = false;
  }
  // SNV in body
  for (int i=1; i < mutated.size() - 1; i++) {
    if (mutated[i-1] == nonMutated[i-1 + ohang] &&
        mutated[i] != nonMutated[i + ohang] &&
        mutated[i+1] == nonMutated[i+1 + ohang] ) {
      alignment->mutatedBases.push_back(mutated[i]);
      alignment->nonMutatedBases.push_back(nonMutated[i+ohang]);
      alignment->SNVLocations.push_back(calibrateWithCIGAR(i + ohang, cp));

      noSNVs = false;
    }
  }
  if (noSNVs) {
    delete alignment;
    alignment = nullptr;
  }
}
bool GenomeMapper::compareSNVLocations(const single_snv &a, const single_snv &b) {
  return a.position < b.position;
}

bool GenomeMapper::equalSNVs(single_snv const& a, single_snv const& b) {
  return (a.position == b.position) &&
         (a.healthy_base == b.healthy_base) &&
         (a.mutation_base == b.mutation_base) &&
         (a.chr == b.chr);
}

void GenomeMapper::outputSNVToUser(vector<SamEntry*> &alignments, string outName) {
  vector<single_snv> separate_snvs;
  for(SamEntry* & entry : alignments) {
    if (entry == nullptr) {
      continue;
    }
    for(int i=0; i < entry->SNVLocations.size(); i++) {
      int alignment_index = get<int>(SamEntry::POS, entry);
      single_snv snv;
      snv.chr = get<string>(SamEntry::RNAME, entry);
      snv.position = alignment_index + entry->SNVLocations[i];;
      snv.healthy_base = entry->nonMutatedBases[i];
      snv.mutation_base = entry->mutatedBases[i];
      separate_snvs.push_back(snv);
    }
    delete entry;
    entry = nullptr;
  }
  // sort the snvs 
  std::sort(separate_snvs.begin(), separate_snvs.end(), compareSNVLocations);
  auto last = std::unique(separate_snvs.begin(), separate_snvs.end(), &equalSNVs);
  separate_snvs.erase(last, separate_snvs.end());
  
  ofstream report(outName);
  report << "Mut_ID\tType\tChr\tPos\tNormal_NT\tTumor_NT\n";
  int64_t i=0;
  for(single_snv &snv : separate_snvs) {
    report << i << "\t" << "SNV\t" << snv.chr << "\t"
           << snv.position << "\t"
           << snv.healthy_base << "\t" 
           << snv.mutation_base
           << "\n";
    i++;
  }
  report.close();
}



/*
void GenomeMapper::callBWA() {
  cout << "Calling bwa..." << endl;

  string command_aln = 
    "./bwa/bwa aln -t 16 /data/ic711/insilico_data/smufin_provided_data/ref_genome/hg19.fa /data/ic711/result/cns_pairs.fastq > /data/ic711/result/cns_pairs.sai";

  string command_sampe = 
    "./bwa/bwa samse /data/ic711/insilico_data/smufin_provided_data/ref_genome/hg19.fa /data/ic711/result/cns_pairs.sai /data/ic711/result/cns_pairs.fastq > /data/ic711/result/cns_pairs.sam";

  system(command_aln.c_str());  
  system(command_sampe.c_str());

  cout << "Finished bwa call." << endl;
}
void GenomeMapper::buildConsensusPairs() {
  // generate consensus pair for each breakpoint block
  // and then add starting gaps to align sequence pair

  consensus_pairs.reserve(snvId->getSize()); // make room
  int continued{0};
  for (int i=0; i < snvId->getSize(); ++i) {
    consensus_pair pair;
    pair.left_ohang = pair.right_ohang = 0;   // default to no overhang

    bool skip_mutated{false}, skip_non_mutated{false};
    skip_mutated = snvId->generateConsensusSequence(i, pair.mut_offset, TUMOUR, pair.pair_id, pair.mutated, pair.mqual);
    skip_non_mutated = snvId->generateConsensusSequence(i, pair.nmut_offset, HEALTHY, pair.pair_id, pair.non_mutated, pair.nqual);

    // discard sequences that do not contain both a non-mutated
    // and mutated cns pair
    if(skip_mutated == true  || skip_non_mutated == true) {
      continued++;
      continue;
    }

    //cout << "Cancer: " << endl;
    //cout << pair.mutated << endl
    //     << pair.mqual << endl;
    //cout << "Healthy: " << endl;
    //cout << pair.non_mutated << endl
    //     << pair.nqual << endl;
    trimCancerConsensus(pair);                // trim extra cancer sequence

    // TURN LQF OFF MASK
    bool low_quality_block = false;
    maskLowQualityPositions(pair, low_quality_block);
    if (low_quality_block) {
      continue;
    }
    consensus_pairs.push_back(pair);
  }
  cout << "Skipped " << continued << endl;
  consensus_pairs.shrink_to_fit();
}

void GenomeMapper::maskLowQualityPositions(consensus_pair & pair, bool &low_quality) {
  int num_low_quality_positions{0};
  for (int pos=0; pos < pair.mutated.size(); pos++) {
    if (pair.mqual[pos] != '-' || pair.nqual[pos + pair.left_ohang] != '-') {
        pair.mutated[pos] = pair.non_mutated[pos + pair.left_ohang];
        num_low_quality_positions++;
    }
  }

  //if (num_low_quality_positions > MAX_LOW_CONF_POSITIONS) low_quality = true;
}
//void GenomeMapper::buildConsensusPairs() {
//  // generate consensus pair for each breakpoint block
//  // and then add starting gaps to align sequence pair
//
//  consensus_pairs.reserve(BPB->getSize()); // make room
//  int continued{0};
//  for (int i=0; i < BPB->getSize(); ++i) {
//
//    vector< vector<int> > tumour_base_frequency, healthy_base_frequency;
//    consensus_pair pair;
//    pair.left_ohang = pair.right_ohang = 0;   // default to no overhang
//
//    pair.mutated = BPB->generateConsensusSequence(i,
//      pair.mut_offset, TUMOUR, pair.pair_id, tumour_base_frequency);
//
//    pair.non_mutated = BPB->generateConsensusSequence(i,
//        pair.nmut_offset, HEALTHY, pair.pair_id, healthy_base_frequency);
//
//    // discard sequences that do not contain both a non-mutated
//    // and mutated cns pair
//    if(pair.mutated == "\0" || pair.non_mutated == "\0") {
//      continued++;
//      continue;
//    }
//
//    trimCancerConsensus(pair);                // trim extra cancer sequence
//
//    //// TURN OFF MASK

//    bool low_quality_block = false;

//    maskLowConfidencePositions(pair, healthy_base_frequency, 
//        tumour_base_frequency, low_quality_block);
//
//    if (low_quality_block) {
//      continue;
//    }
//    consensus_pairs.push_back(pair);
//  }
//  cout << "Skipped " << continued << endl;
//  consensus_pairs.shrink_to_fit();
//}

//void GenomeMapper::maskLowConfidencePositions(consensus_pair &pair,
//                                vector< vector<int> > &healthy_base_freq,
//                                vector< vector<int> > &tumour_base_freq,
//                                bool &discard) {
//  discard = false;
//  unsigned int start_h= 0, start_t= 0;
//
//  // start points do not take into account the updated
//  // matrix - which should have been cleaved
//  for(int i=0; i < healthy_base_freq[0].size(); i++) {
//    if(healthy_base_freq[0][i] != -1) {
//      start_h = i; // + left_ohang
//      break;    
//    }
//  }
//  for(int i=0; i < tumour_base_freq[0].size(); i++) {
//    if(tumour_base_freq[0][i] != -1) {
//      start_t = i; // should be corrected based on -1 updated matrix
//      break;
//    }
//  }
//
//  // mask based on tumour cns
//  for(int pos = start_t; pos < pair.mutated.size() + start_t; pos++) {
//    int n_tumour_bases_above_err_freq = 0;
//
//    if(tumour_base_freq[0][pos] == -1) {    // is this squared to the final base
//      break;
//    }
//
//    // get total reads
//    double total_reads = 0;
//    for(int base=0; base < 4; base++) {
//      total_reads += tumour_base_freq[base][pos];
//    }
//
//    // calc number of bases over the error frequency
//    for(int base=0; base < 4; base++) {
//      if(tumour_base_freq[base][pos] / total_reads > ALLELIC_FREQ_OF_ERROR) {
//        n_tumour_bases_above_err_freq++;
//      }
//    }
//
//    // if the number of bases with a high allelic frequency is above
//    // one, then the position is of low condifence, so mask
//    if (n_tumour_bases_above_err_freq > 1) {
//      pair.mutated[pos - start_t] = pair.non_mutated[pos - start_t + pair.left_ohang];
//    }
//  }
//
//  int number_of_low_conf_positions = 0; 
//  // mask based on healthy cns
//  for(int pos = start_h + pair.left_ohang; pos < pair.mutated.size() +
//      start_h + pair.left_ohang; pos++) {
//
//    int n_healthy_bases_above_err_freq = 0;
//    if(healthy_base_freq[0][pos] == -1) {
//      break;
//    }
//
//    // get total reads
//    double total_reads = 0;
//    for(int base = 0; base < 4; base++) {
//      // am I not counting over the non-updated start position (should be 
//      // pos + left_ohang, so i'm counting the cancer positions)
//      total_reads += healthy_base_freq[base][pos];
//    }                                                 
//
//    // cals number of bases over the error frequency
//    for(int base=0; base < 4; base++) {
//      if(healthy_base_freq[base][pos] / total_reads > ALLELIC_FREQ_OF_ERROR) {
//        n_healthy_bases_above_err_freq++;
//      }
//    }
//
//    // if n bases with high allelic freq. is above one, then 
//    // position is low confidence so mask
//    if(n_healthy_bases_above_err_freq > 1) {
//      pair.mutated[pos - pair.left_ohang - start_h] = pair.non_mutated[pos -
//        start_h];
//      number_of_low_conf_positions++;
//    }
//  }
//
//
//  // if there were too many low condidence positions, 
//  // the block is low quality, so discard
//  if(number_of_low_conf_positions > MAX_LOW_CONF_POSITIONS) {
//    discard = true;
//  }
//  cout << "number of low conf: " << number_of_low_conf_positions << endl;
//}
void GenomeMapper::trimCancerConsensus(consensus_pair & pair) {
  // Trim portions of the cancer consensus sequence that are
  // longer than heathy. Leave healthy if longer.
  // left cleave
  if (pair.mut_offset > pair.nmut_offset) {
    pair.mutated.erase(0, pair.mut_offset - pair.nmut_offset);
    pair.mqual.erase(0, pair.mut_offset - pair.nmut_offset);
  }
  else if (pair.mut_offset < pair.nmut_offset) {
    pair.left_ohang = pair.nmut_offset - pair.mut_offset;
  }

  // right cleave
  if (pair.mutated.size() > (pair.non_mutated.size() - pair.left_ohang)) {
    int dist = pair.mutated.size() - (pair.non_mutated.size() - pair.left_ohang);
    pair.mutated.erase(pair.mutated.size() - dist, dist);
    pair.mqual.erase(pair.mqual.size() - dist, dist);
  }
  else if (pair.mutated.size() < (pair.non_mutated.size() - pair.left_ohang)) {
    pair.right_ohang = (pair.non_mutated.size() - pair.left_ohang) - pair.mutated.size();
  }
  // else, equal length. Do nothing
}
void GenomeMapper::printMutation(char healthy, char cancer, ofstream &mut_file) {
  mut_file << healthy <<  ", " << cancer << endl;
}

void GenomeMapper::printAlignmentStructs(vector<SamEntry> &alignments) {
  for (SamEntry &entry: alignments) {
    cout << "Cancer seq:" << endl;
    if(entry.get<int>(SamEntry::FLAG) == FORWARD_FLAG) {
      printGaps(entry.get<int>(SamEntry::LEFT_OHANG));
    }
    else {
      printGaps(entry.get<int>(SamEntry::RIGHT_OHANG));
    }

    cout << entry.get<string>(SamEntry::HDR) << endl;
    cout << "Healthy seq:" << endl;
    cout << entry.get<string>(SamEntry::SEQ) << endl;
    cout << "Left ohang:  " << entry.get<int>(SamEntry::LEFT_OHANG) << endl;
    cout << "Right ohang: " << entry.get<int>(SamEntry::RIGHT_OHANG) << endl;
    cout << "Pair id: " << entry.get<int>(SamEntry::BLOCK_ID) << endl;
    for(int i = 0; i < entry.snvLocSize(); i++) {
      cout << entry.snvLocation(i) << ", ";
    }
    cout << endl << endl;
  }
}

void GenomeMapper::printConsensusPairs() {
  for(consensus_pair cns_pair : consensus_pairs) {
    cout << "Mutation Sequence:" << endl;
    cout << cns_pair.mutated << endl;
    cout << "Healthy Sequence:" << endl;
    cout << cns_pair.non_mutated << endl;

    cout << "SNV locations" << endl;
    for(int pos : cns_pair.mutations.SNV_pos) {
      cout << pos << ", ";
    }
    cout << endl << endl;
  }
}
void GenomeMapper::printGaps(int gaps) {
  for(; gaps > 0; gaps--) {
    cout << "-";
  }
}

void GenomeMapper::printAllAlignments(vector<SamEntry> &alignments){
  for(SamEntry &entry: alignments) {
    cout << "FLAG  :" << entry.get<int>(SamEntry::FLAG) << endl;
    cout << "CHR  :" << entry.get<string>(SamEntry::RNAME) << endl;
    cout << "POS :" << entry.get<int>(SamEntry::POS) << endl;
    cout << "HEALTHY: " << entry.get<string>(SamEntry::SEQ) << endl;
    cout << "TUMOUR : " << entry.get<string>(SamEntry::HDR) << endl;
    for(int i = 0; i < entry.snvLocSize(); i++) {
      cout << entry.snvLocation(i)  << ", ";
    }
    cout << endl << endl;
  }
}

void GenomeMapper::printSingleAlignment(SamEntry &entry) {
  cout << "FLAG  :" << entry.get<int>(SamEntry::FLAG) << endl;
  cout << "CHR  :" << entry.get<string>(SamEntry::RNAME) << endl;
  cout << "POS :" << entry.get<int>(SamEntry::POS) << endl;
  cout << "HELATHY :" << entry.get<string>(SamEntry::SEQ) << endl;
  cout << "TUMOUR :" << entry.get<string>(SamEntry::HDR) << endl;
  for(int i = 0; i < entry.snvLocSize(); i++) {
    cout << entry.snvLocation(i)  << ", ";
  }
  cout << endl << endl;
}


void GenomeMapper::correctReverseCompSNV(vector<SamEntry> &alignments) {

  for(SamEntry &entry : alignments) {

    if(entry.get<int>(SamEntry::FLAG) == FORWARD_FLAG) {
      for(int i = 0; i < entry.snvLocSize(); i++) {
        entry.setSNVLocation(i, entry.snvLocation(i) - 1);
      }
    }
    else if(entry.get<int>(SamEntry::FLAG) == REVERSE_FLAG) { // convert indecies to rev comp and rev comp cns
      entry.get<string>(SamEntry::HDR) = reverseComplementString(entry.get<string>(SamEntry::HDR));
      for(int i = 0; i < entry.snvLocSize(); i++) {
        entry.setSNVLocation(i, entry.get<string>(SamEntry::HDR).size() - entry.snvLocation(i));
      }
    }
  }
}

string GenomeMapper::generateParseString(mutation_classes &m) {
  // this codes the mutations found in a consensus sequence
  // the format is [SNV:a;b;c][SSV:e;f;g][LSV:x;y;z]
  // lower case characters represent the index values of
  // various mutations. 

  string mutation_string;
  mutation_string = "[SNV:";
  // code all the SNV positions into the string
  // at the same time, converting the positions 
  // from 0-based to 1-based

  for (int snv_index : m.SNV_pos) {
    mutation_string += to_string(snv_index+1) + ";";
  }
  mutation_string.pop_back();	// remove trailing ";"
  mutation_string += "][SSV:][LSV:]";
  return mutation_string;
}


string GenomeMapper::reverseComplementString(string const& s){
  string revcomp = "";
  for(int i = s.size()-1; i >= 0; i--) {
  // travel in reverse and switch for complementary
    switch(s[i]) {
      case 'A': revcomp += "T"; break;
      case 'T': revcomp += "A"; break;
      case 'C': revcomp += "G"; break;
      case 'G': revcomp += "C"; break;
    }
  }
  return revcomp;
}
//void GenomeMapper::countSNVs() {
//
//
//
//  for(consensus_pair &p : consensus_pairs) {
//    bool discard = false; 
//    // remove any consensus pairs that
//    // do not contain SNV style sequence differences
//    // these could result from indels
//
//    for(int i=0; i < p.mutated.size() - 1; i++){
//      if (p.mutated[i] != p.non_mutated[i] && 
//	  p.mutated[i+1] != p.non_mutated[i+1]) {
//          discard = true;
//      }
//    }
//    
//    if(discard) {
//      continue;
//    }
//
//
//    // identify SNV occuring at the start of the string
//    if (p.mutated[0] != p.non_mutated[0] && 
//        p.mutated[1] == p.non_mutated[1]) {
//
//      p.mutations.SNV_pos.push_back(0);
//    }
//
//    // identify SNVs occuring within the string
//    for(int i=1; i < p.mutated.size() - 1; i++) {
//      if(p.mutated[i] != p.non_mutated[i] &&        // char diffrent but
//          p.mutated[i+1] == p.non_mutated[i+1] &&   // ...next char same
//          p.mutated[i-1] == p.non_mutated[i-1]) {   // ...prev char same
//
//        // then count as SNV
//        p.mutations.SNV_pos.push_back(i);   // store index of variants
//      }
//    }
//
//    // identify SNV occuring at the very end of the string
//    if (p.mutated[p.mutated.size()-1] !=
//        p.non_mutated[p.non_mutated.size()-1] &&
//        p.mutated[p.mutated.size()-2] == 
//        p.non_mutated[p.non_mutated.size()-2]){
//      
//      p.mutations.SNV_pos.push_back(p.mutated.size()-1);
//    }
//  }
//}

 */
