/*
SNVIdentifier.cpp
Author: Izaak Coleman
*/

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <set>
#include <utility>
#include <iterator>
#include <map>
#include <list>
#include <limits>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <memory>
#include <iterator>

#include <zlib.h>
#include <omp.h>

#include "divsufsort64.h"
#include "util_funcs.h"
#include "SNVIdentifier.h"
#include "GenomeMapper.h"

#include "benchmark.h"
using namespace std;

const char TERM = '$';



SNVIdentifier::SNVIdentifier(GSA &_gsa,
                                     string outpath,
                                     string const & basename,
                                     char mpq, int g1, int g2, int cut,
                                     int t, int mlcp, double e, double a):
                                     MIN_PHRED_QUAL(mpq), 
                                     GSA1_MCT(g1), 
                                     GSA2_MCT(g2), 
                                     COVERAGE_UPPER_THRESHOLD(cut),
                                     N_THREADS(t),
                                     MAX_LOW_CONFIDENCE_POS(mlcp),
                                     ECONT(e), 
                                     ALLELIC_FREQ_OF_ERROR(a) {
  if (GSA1_MCT > GSA2_MCT) GSA2_MCT = GSA1_MCT;
  if (outpath[outpath.size()-1] != '/') outpath += "/";
  gsa = &_gsa;    
  cout << "n_threads" << N_THREADS << endl;
  cout << "READ_LENGTH: " << gsa->get_max_read_len() << endl;
  cout << "GSA1_MCT : " << GSA1_MCT  << endl;
  cout << "GSA2_MCT: " << GSA2_MCT << endl;
  cout << "UBT: " << COVERAGE_UPPER_THRESHOLD << endl;
  cout << "ALLELIC_FREQ_OF_ERROR: " << ALLELIC_FREQ_OF_ERROR << endl;
  cout << MIN_PHRED_QUAL << endl;
  cout << "Extracting cancer-specific reads..." << endl;
  START(extractCancerSp);
  extractCancerSpecificReads(); 
  COMP(extractCancerSp);
  cout << "No of extracted reads: " << CancerExtraction.size() << endl;
  // Group blocks covering same mutations in both orientations
  cout << "Seeding breakpoint blocks by constructing GSA2..." << endl;
  START(seedBreakPointBlocks);
  seedBreakPointBlocks();
  COMP(seedBreakPointBlocks);
  CancerExtraction.clear();
  cout << "Number of seed blocks: " << SeedBlocks.size() << endl;
  cout << "<<<<<<<<<<<<<<building consensus pairs " << endl;

  START(CNS_construction);
  int64_t elementsPerThread;
  int n_threads;
  if (SeedBlocks.size() < N_THREADS) {
    omp_set_num_threads(1);
    n_threads = 1;
    elementsPerThread = SeedBlocks.size();
  } else {
    omp_set_num_threads(N_THREADS);
    n_threads = N_THREADS;
    elementsPerThread = SeedBlocks.size() / n_threads;
  }

  std::list<consensus_pair> consensus_pairs;
#pragma omp parallel for ordered
  for (int i = 0; i < n_threads; i++) {
    shared_ptr<bpBlock> * from = (&SeedBlocks[0] + i*elementsPerThread);
    shared_ptr<bpBlock> * to;
    if (i == n_threads - 1) {
      to = (&SeedBlocks[0] + SeedBlocks.size());
    } else {
      to = (&SeedBlocks[0] + (i+1)*elementsPerThread);
    }
    vector<consensus_pair> threadWork;
    buildConsensusPairsWorker(from, to, threadWork);
#pragma omp ordered
    {
      consensus_pairs.insert(consensus_pairs.end(), threadWork.begin(), threadWork.end());
    }
  }
  COMP(CNS_construction);
  // pop the consensus_pair elements, and for each pop, add an element to
  // tumour_cns and fastq_string. 
  string fastq_data;
  buildFastqAndTumourCNSData(consensus_pairs, fastq_data);

  cout << "<<<<<<<<<<<<<<Finished break point block construction" << endl;
  cout << "Writing consensus pairs as fastq" << endl;
  writeConsensusPairs(fastq_data, outpath + basename + ".fastq.gz");

}

void SNVIdentifier::buildFastqAndTumourCNSData(std::list<consensus_pair> & consensus_pairs, string & fastq) {
  int64_t index = 0;
  tumour_cns.reserve(consensus_pairs.size());
  while (!consensus_pairs.empty()) {
    consensus_pair pair = consensus_pairs.front();
    tumour_cns.push_back(pair.mutated);
    string qual(pair.non_mutated.size(), '!');
    fastq += "@" + std::to_string(index)  + ";" + to_string(pair.left_ohang) + 
             ";" + to_string(pair.right_ohang) + "\n" + pair.non_mutated 
             + "\n+\n" + qual + "\n";
    ++index;
    consensus_pairs.pop_front();
  }
}

void SNVIdentifier::writeConsensusPairs(string const & fastq_data, string const & fname) {
  gzFile outfile = gzopen(fname.c_str(), "wb");
  int64_t src_sz = fastq_data.size();
  int chunk = 0x4000;
  char * src = const_cast<char*>(fastq_data.c_str());
  char * it = src;
  for (; (it + chunk)  < (src + src_sz); it += chunk) {
    if(gzwrite(outfile, it, chunk) != chunk) {
      cout << "WRITE ERROR" << endl;
    }
  }
  gzwrite(outfile, it, ((src + src_sz) - it));
  gzclose(outfile);
}
void SNVIdentifier::buildConsensusPairsWorker(shared_ptr<bpBlock> * block,
    shared_ptr<bpBlock> * end,
    vector<consensus_pair> & threadWork) {
//START(SNVIdentifier_buildConsensusPairsWorker);
  for(; block < end; block++) {
    if ((*block)->size() > COVERAGE_UPPER_THRESHOLD) {
      *block = nullptr;
      continue;
    }
    consensus_pair pair;
    pair.left_ohang = pair.right_ohang = 0;
    generateConsensusSequence(TUMOUR, **block, pair.mut_offset, pair.mutated, pair.mqual);


    extractNonMutatedAlleles(**block, pair);

    if ((*block)->size() > COVERAGE_UPPER_THRESHOLD) {
      *block = nullptr;
      continue;
    }
    generateConsensusSequence(HEALTHY, **block, pair.nmut_offset, pair.non_mutated, pair.nqual);
    *block = nullptr;
    if (pair.mutated.empty() || pair.non_mutated.empty()) {
      continue;
    }
    if (excessLowQuality(pair)) {
      continue;
    }
    trimHealthyConsensus(pair); // MUST trim healthy first
    trimCancerConsensus(pair);
    if (pair.mutated.empty() || pair.non_mutated.empty()) {
      continue;
    }
    maskLowQualityPositions(pair);

    if (noSNV(pair)) {
      continue;
    }
    threadWork.push_back(pair);
    //string qual(pair.non_mutated.size(), '!');
    //threadWork += "@" + pair.mutated + "[" + to_string(pair.left_ohang) + 
    //           ";" + to_string(pair.right_ohang) + "]\n" + pair.non_mutated 
    //           + "\n+\n" + qual + "\n";
  }
}

void SNVIdentifier::printAlignedBlock(bpBlock block) {
  if (block.size() > COVERAGE_UPPER_THRESHOLD || block.empty()) return;
  int max = std::numeric_limits<int>::min();

  for (read_tag tag : block) {
    if (tag.orientation == LEFT) {
      tag.offset = convertOffset(tag);
    }
    if (tag.offset > max) {
      max = tag.offset;
    }
  }

  cout << "TUMOUR subblock " << endl;
  for (read_tag tag : block) {
    if (tag.tissue_type == HEALTHY || tag.tissue_type == SWITCHED) continue;
    string read;
    int offset = tag.offset;
    if (tag.orientation == LEFT) {
      read = reverseComplementString(gsa->get_suffix_string(tag.read_id));
      offset = convertOffset(tag);
    } 
    else {
      read = gsa->get_suffix_string(tag.read_id);
      read.pop_back();
    }

    cout << addGaps(max - offset) << read << " - " 
         << tag.read_id << " - " 
         << tag.offset << " - " 
         << ((tag.orientation) ? "R" : "L") << " - "
         << ((tag.tissue_type == HEALTHY) ? "H" : ((tag.tissue_type == SWITCHED) ? "S" : "T")) << endl;
  }

  cout << "HEALTHY subblock " << endl;
  for (read_tag tag : block) {
    if (tag.tissue_type == TUMOUR) continue;
    string read;
    int offset = tag.offset;
    if (tag.orientation == LEFT) {
      read = reverseComplementString(gsa->get_suffix_string(tag.read_id));
      offset = convertOffset(tag);
    } 
    else {
      read = gsa->get_suffix_string(tag.read_id);
      read.pop_back();
    }

    cout << addGaps(max - offset) << read << " - " 
         << tag.read_id << " - " 
         << tag.offset << " - " 
         << ((tag.orientation) ? "R" : "L") << " - "
         << ((tag.tissue_type == HEALTHY) ? "H" : ((tag.tissue_type == SWITCHED) ? "S" : "T")) << endl;
  }
}
string SNVIdentifier::addGaps(int n_gaps) {
  string gaps(n_gaps, '-');
  return gaps;
}

bool SNVIdentifier::noSNV(consensus_pair const & pair) {
  if (pair.mutated[0] != pair.non_mutated [0 + pair.left_ohang] &&
      pair.mutated[1] == pair.non_mutated[1 + pair.left_ohang]) {
    return false;
  }
  // SNV at end 
  int cnsLen = pair.mutated.size();
  if (pair.mutated[cnsLen-1] != pair.non_mutated[cnsLen-1 + pair.left_ohang] &&
      pair.mutated[cnsLen-2] == pair.non_mutated[cnsLen-2 + pair.left_ohang]) {
      return false;
  }
  // SNV in body
  for (int i=1; i < pair.mutated.size() - 1; i++) {
    if (pair.mutated[i-1] == pair.non_mutated[i-1 + pair.left_ohang] &&
        pair.mutated[i] != pair.non_mutated[i + pair.left_ohang] &&
        pair.mutated[i+1] == pair.non_mutated[i+1 + pair.left_ohang] ) {
      return false;
    }
  }
  return true;
}


bool SNVIdentifier::excessLowQuality(consensus_pair & pair) {
  int numLowQuality{0};
  for (int i = 0; i < pair.mutated.size(); i++) {
    numLowQuality += (pair.mqual[i] == 'L');
  }
  for (int i = 0; i < pair.non_mutated.size(); i++) {
    numLowQuality += (pair.nqual[i] == 'L' || pair.nqual[i] == 'B');
  }
  if (numLowQuality > MAX_LOW_CONFIDENCE_POS) return true;
  return false;
}

void SNVIdentifier::maskLowQualityPositions(consensus_pair & pair) {
  for (int pos=0; pos < pair.mutated.size(); pos++) {
    if (pair.mqual[pos] != '-' || pair.nqual[pos + pair.left_ohang] != '-') {
      pair.mutated[pos] = pair.non_mutated[pos + pair.left_ohang];
    }
  }
}

void SNVIdentifier::trimCancerConsensus(consensus_pair & pair) {
  if (pair.mut_offset > pair.nmut_offset) {
    pair.mutated.erase(0, pair.mut_offset - pair.nmut_offset);
    pair.mqual.erase(0, pair.mut_offset - pair.nmut_offset);
  }
  else if (pair.mut_offset < pair.nmut_offset) {
    pair.left_ohang = pair.nmut_offset - pair.mut_offset;
  }
  if (pair.mutated.size() > (pair.non_mutated.size() - pair.left_ohang)) {
    int dist = pair.mutated.size() - (pair.non_mutated.size() - pair.left_ohang);
    pair.mutated.erase(pair.mutated.size() - dist, dist);
    pair.mqual.erase(pair.mqual.size() - dist, dist);
  }
  else if (pair.mutated.size() < (pair.non_mutated.size() - pair.left_ohang)) {
    pair.right_ohang = (pair.non_mutated.size() - pair.left_ohang) - pair.mutated.size();
  }
}

void SNVIdentifier::trimHealthyConsensus(consensus_pair & pair) {
  int left_arrow = 0, right_arrow = pair.nqual.size() - 1;
  for(; right_arrow >= 0 && 
        (pair.nqual[right_arrow] == 'T' ||
         pair.nqual[right_arrow] == 'B'); right_arrow--);
  right_arrow++;
  pair.nqual.erase(right_arrow);
  pair.non_mutated.erase(right_arrow);

  for(; left_arrow < pair.nqual.size() && 
        (pair.nqual[left_arrow] == 'T' ||
         pair.nqual[left_arrow] == 'B'); left_arrow++);
  pair.nqual.erase(0, left_arrow);
  pair.non_mutated.erase(0, left_arrow);
  pair.nmut_offset -= left_arrow;
}

void SNVIdentifier::extractNonMutatedAlleles(bpBlock &block,
    consensus_pair &pair) {
  bool LEFT{false}, RIGHT{true};
  string query = pair.mutated.substr(pair.mut_offset, gsa->get_min_suf_size());
  string rcquery = reverseComplementString(query);
  int64_t fwd_idx = binarySearch(query);
  int64_t rev_idx = binarySearch(rcquery);

  bool success_left{false}, success_right{false};
  if (fwd_idx != -1) {
    success_right = extendBlock(fwd_idx, block, RIGHT, 0);
  }
  if (rev_idx != -1) {
    success_left = extendBlock(rev_idx, block, LEFT, 0);
  }
  if (!(success_left || success_right)) {
    for (int i = pair.mut_offset - gsa->get_min_suf_size();
        i <= pair.mut_offset + gsa->get_min_suf_size();
        i += gsa->get_min_suf_size() * 2) { // * 2 skips center search
      if (i < 0 || i > pair.mutated.size() - gsa->get_min_suf_size()) {
        continue;
      }
      query = pair.mutated.substr(i, gsa->get_min_suf_size());
      rcquery = reverseComplementString(query);
      fwd_idx = binarySearch(query);
      rev_idx = binarySearch(rcquery);
      if (fwd_idx != -1) {
        extendBlock(fwd_idx, block, RIGHT, pair.mut_offset - i);
      }
      if (rev_idx != -1) {
        extendBlock(rev_idx, block, LEFT, pair.mut_offset - i);
      }
    }
  }
}
char SNVIdentifier::rc(char c, int dir) {
  if (dir == 1) {
    return c;
  }
  switch (c) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';
  }
  return 'x';
}

void SNVIdentifier::generateConsensusSequence(bool tissue,
    bpBlock const& block, int & cns_offset, string & cns,
    string & qual) {
  vector<read_tag> subBlock;            // work with subset
  for (read_tag const& tag : block) {
    if (tissue == HEALTHY  && 
       (tag.tissue_type == HEALTHY || tag.tissue_type == SWITCHED)) {
      subBlock.push_back(tag);
    }
    else if(tissue == TUMOUR && tag.tissue_type == TUMOUR) {
      subBlock.push_back(tag);
    }
  }
  int max_offset = 0;
  int min_offset = std::numeric_limits<int>::max();
  for (read_tag &tag : subBlock) {
    if (tag.orientation == LEFT) {
      tag.offset = convertOffset(tag);
    }
    if (tag.offset > max_offset) {
      max_offset = tag.offset;
    }
    if (tag.offset < min_offset) {
      min_offset = tag.offset;
    }
  }
  if(min_offset == std::numeric_limits<int>::max()) {
    min_offset = 0;
  }


  unsigned int width = max_offset + gsa->get_max_read_len() - min_offset;
  vector<int> cnsCount(4*width, 0);
  for (read_tag const & tag : subBlock) {

    string::const_iterator readPtr = gsa->suffix_at(tag.read_id);
    string::const_iterator readPtr_end = gsa->suffix_at(tag.read_id) + gsa->len(tag.read_id) - 2;
    string::const_iterator phredPtr = gsa->phred_at(tag.read_id);
    string::const_iterator phredPtr_end = gsa->phred_at(tag.read_id) + gsa->len(tag.read_id) - 2;
    int d= 1;
    if (tag.orientation == LEFT) {
      // swap pointers
      string::const_iterator rtemp, ptemp;
      rtemp = readPtr; readPtr = readPtr_end; readPtr_end = rtemp;
      ptemp = phredPtr; phredPtr = phredPtr_end; phredPtr_end = ptemp;
      d= -1;
    }

    int i = 0;
    for (; (readPtr + i*d) != readPtr_end; i++) {
        cnsCount[
          ((max_offset - tag.offset + i)*4    ) * (rc(*(readPtr + i*d),d) == 'A') +  
          ((max_offset - tag.offset + i)*4 + 1) * (rc(*(readPtr + i*d),d) == 'T') +
          ((max_offset - tag.offset + i)*4 + 2) * (rc(*(readPtr + i*d),d) == 'C') +
          ((max_offset - tag.offset + i)*4 + 3) * (rc(*(readPtr + i*d),d) == 'G') 
        ] += 1 * (*(phredPtr + i*d) >= MIN_PHRED_QUAL);
    }
    cnsCount[
      ((max_offset - tag.offset + i)*4    ) * (rc(*(readPtr + i*d),d) == 'A') +  
      ((max_offset - tag.offset + i)*4 + 1) * (rc(*(readPtr + i*d),d) == 'T') +
      ((max_offset - tag.offset + i)*4 + 2) * (rc(*(readPtr + i*d),d) == 'C') +
      ((max_offset - tag.offset + i)*4 + 3) * (rc(*(readPtr + i*d),d) == 'G') 
    ] += 1 * (*(phredPtr + i*d) >= MIN_PHRED_QUAL);
  }

  for (int pos=0; pos < width; pos++) {
    int maxVal = 0, maxInd = 0, m = -1;
    for (int base=0; base < 4; ++base) {
      bool cond = (cnsCount[pos*4 + base] > maxVal);
      maxVal = (cnsCount[pos*4 + base] & m*cond) | (maxVal & m*!cond);
      maxInd = (base & m*cond) | (maxInd & m*!cond);
    }
    if (maxVal == 0) {
      maxInd = assignBaseDisregardingPhred(pos, subBlock, max_offset);
    }
    switch(maxInd) {
      case 0: cns += "A"; break;
      case 1: cns += "T"; break;
      case 2: cns += "C"; break;
      case 3: cns += "G"; break;
    }
  }
  cns_offset = max_offset;
  string q_str(width, '-');
  for (int pos=0; pos < width; pos++) {
    // Mask position if there is < GSA_MCT2 reads contributing to consensus
    int nreads=0;
    for (int base=0; base < 4; base++) nreads += cnsCount[pos*4 + base];
    if (nreads < GSA2_MCT) {
      if (tissue == TUMOUR) {
        q_str[pos] = 'X';
      }
      else { // tissue == SWITCHED || tissue == HEALTHY
        q_str[pos] = 'T';
      }
    }
  }
  buildQualityString(q_str, cnsCount, cns, tissue);
  qual = std::move(q_str);
}

int SNVIdentifier::assignBaseDisregardingPhred(int const pos, vector<read_tag> const &
    block, int const max_offset) {
  int c[4] = {};
  string a = "ATCG";
  for (read_tag const & r : block) {
    int idx = r.offset - max_offset + pos;
    if (r.orientation == LEFT) {
      idx = gsa->len(r.read_id) - idx - 2;
    }
    if (idx < 0 || idx > gsa->len(r.read_id)-1) continue;
    string::const_iterator p = gsa->suffix_at(r.read_id) +  idx;
    ++c[a.find(rc(*p, (r.orientation == RIGHT)))];
  }
  int maxVal = 0, maxInd = 0;
  for (int i = 0; i < 4; ++i) {
    if (c[i] > maxVal) {
      maxVal = c[i];
      maxInd = i;
    }
  }
  return maxInd;
}

void SNVIdentifier::buildQualityString(string & qual, vector<int> const& freq_matrix,
    string const& cns, bool tissue) {
  for (int pos=0; pos < cns.size(); pos++) {
    if (tissue == TUMOUR && qual[pos] == '-') { // do not override 'X'
      int base{0};
      switch (cns[pos]) {
        case 'A': base = 0; break;
        case 'T': base = 1; break;
        case 'C': base = 2; break;
        case 'G': base = 3; break;
      }
      if (freq_matrix[pos*4 + base] < GSA2_MCT) {
        qual[pos] = 'C';
        continue;
      }
    }
    double total_bases = 0;
    for (int base = 0; base < 4; base++) total_bases += freq_matrix[pos*4 + base];
    int n_bases_above_err_freq{0};
    for(int base = 0; base < 4; base++) {
      if ((freq_matrix[pos*4 + base] / total_bases) > ALLELIC_FREQ_OF_ERROR) {
        n_bases_above_err_freq++;
      }
    }
    if (n_bases_above_err_freq > 1) {
      if (qual[pos] == 'X' || qual[pos] == '-') {
        qual[pos] = 'L';
      }
      else if (qual[pos] == 'T') {
        qual[pos] == 'B';
      }
    }
  }
}

void SNVIdentifier::extractCancerSpecificReads() {
  omp_set_num_threads(N_THREADS);
  int64_t elementsPerThread = (gsa->size() / N_THREADS);
#pragma omp parallel for
  for (unsigned int i=0; i < N_THREADS; i++) {
    int64_t from = i*elementsPerThread;
    int64_t to;
    if (i == N_THREADS - 1) {
      to = gsa->size();
    } else {
      to = (i+1)*elementsPerThread;
    }
    set<int64_t> threadWork;
    extractionWorker(from, to, threadWork);
#pragma omp critical
    {
      CancerExtraction.insert(threadWork.begin(), threadWork.end());
      threadWork.clear();
    }
  }
}

void SNVIdentifier::printExtractedCancerReads() {
  ofstream of("/data/ic711/extracted_cancer_read_ids.txt");
  for (unsigned int rid : CancerExtraction) of << rid << '\n';
  of << endl;
  of.close();
}

int64_t SNVIdentifier::backUpSearchStart(int64_t seedIndex) {
  if (seedIndex == 0) {
    return seedIndex;
  }
  int64_t startPoint = seedIndex;
  while (computeLCP(gsa->suffix_at(gsa->sa_element(startPoint)),
        gsa->suffix_at(gsa->sa_element(seedIndex - 1)))
      >= gsa->get_min_suf_size()) {
    seedIndex--;
    if (seedIndex == 0) break;
  }
  return seedIndex;
}

int64_t SNVIdentifier::computeLCP(string::const_iterator a, string::const_iterator b) {
  int64_t lcp{0};
  while (*a != TERM && *b != TERM && *a == *b) {lcp++; a++; b++;}
  return lcp;
}

void SNVIdentifier::extractionWorker(int64_t seed_index, int64_t to, set<int64_t> & threadExtr) {
  seed_index = backUpSearchStart(seed_index);
  int64_t extension {seed_index + 1};
  while (seed_index < to && seed_index != gsa->size() - 1) {
    double c_reads{0}, h_reads{0};    // reset counts
    // Assuming that a > 2 group will form, start counting from seed_index
    if (gsa->tissuetype(gsa->sa_element(seed_index)) == HEALTHY) h_reads++;
    else c_reads++;
    while (computeLCP(gsa->suffix_at(gsa->sa_element(seed_index)),
          gsa->suffix_at(gsa->sa_element(extension)))
        >= gsa->get_min_suf_size()) {
      // tally tissue types of group
      if(gsa->tissuetype(gsa->sa_element(extension)) == HEALTHY) h_reads++;
      else c_reads++;
      extension++;
      if (extension == gsa->size()) break;    // bound check GSA
    }
    // Group size == 1 and group sizes of 1 permitted and groups is cancer read
    //cout << h_reads/c_reads << endl;

    if (extension - seed_index == 1 && GSA1_MCT   == 1 && gsa->tissuetype(seed_index) == TUMOUR) {
      threadExtr.insert(gsa->read_id_of_suffix(gsa->sa_element(seed_index)));           // extract read
    }
    else if (c_reads >= GSA1_MCT  && (h_reads / (c_reads + h_reads)) <= ECONT)  {
     // cout << "Inserted" << endl;
      for (int64_t i = seed_index; i < extension; i++) {
        if (gsa->tissuetype(gsa->sa_element(i)) == TUMOUR) {
          threadExtr.insert(gsa->read_id_of_suffix(gsa->sa_element(i)));
        }
      }
    }
    seed_index = extension++; 
  }

  // On very rare cases, the very last value in the GSA is unique. 
  // In such a case, this value will not have been analysed 
  // by the loop. In these cases, seed_index == SA.size()-1, rather
  // than SA.size() - which is the case if the extension proceeds to the
  // end.  Therefore, if the case condition is true, we need to check
  // it separately.
  if (seed_index == gsa->size() -1) {
    if (GSA1_MCT  == 1 && gsa->tissuetype(gsa->sa_element(seed_index)) == TUMOUR) {
      threadExtr.insert(gsa->read_id_of_suffix(gsa->sa_element(seed_index)));
    }
  }
}


void SNVIdentifier::seedBreakPointBlocks() {
  string concat("");
  vector<pair<int64_t, int64_t>> bsa;
  set<int64_t>::iterator it = CancerExtraction.begin();
  bsa.push_back(pair<int64_t, int64_t>(*it,0));
  concat += gsa->get_suffix_string(*it);
  concat += reverseComplementString(gsa->get_suffix_string(*it));
  concat += '#';

  cout << "Cancer Extraction size: " <<  CancerExtraction.size() << endl;
  int64_t concat_idx = 0;
  it++;
  for (; it != CancerExtraction.end(); it++) {
    concat += gsa->get_suffix_string(*it);
    concat += reverseComplementString(gsa->get_suffix_string(*it));
    concat += '#';
    concat_idx += gsa->len(*std::prev(it)) * 2;
    bsa.push_back(pair<int64_t, int64_t>(*it, concat_idx));
  }

  cout << "Building cancer specific sa" << endl;
  int64_t dSA_sz = concat.size();
  int64_t * dSA = (int64_t*) std::malloc(dSA_sz*sizeof(int64_t));
  divsufsort64((uint8_t*)const_cast<char*>(concat.c_str()), dSA,(int64_t)dSA_sz);
  cout << "Size of dSA before removes_short_suffixes: " << dSA_sz << endl;

  // void function changes both dSA and dSA_sz
  remove_short_suffixes(dSA, dSA_sz, gsa->get_min_suf_size(), concat);
  cout << "Size of dSA after removes_short_suffixes: " << dSA_sz << endl;

  omp_set_num_threads(N_THREADS);
  int64_t elementsPerThread = dSA_sz / N_THREADS;
  set<shared_ptr<bpBlock>, vBlockComp> result;
#pragma omp parallel for ordered
  for (int i = 0; i < N_THREADS; i++) {
    int64_t * from = (dSA + i*elementsPerThread);
    int64_t * to = nullptr;
    if (i == N_THREADS  - 1) {
      to = (dSA + dSA_sz);
    }
    else {
      to = (dSA + (i+1)*elementsPerThread);
    }
    set<shared_ptr<bpBlock> , vBlockComp> twork;
    buildVariantBlocks(dSA, dSA_sz, from, to, bsa, concat, twork);
#pragma omp ordered
  //  SeedBlocks.insert(SeedBlocks.end(), twork.begin(), twork.end());
    result.insert(twork.begin(), twork.end());
  }
  START(seed_b_ins);
  SeedBlocks.insert(SeedBlocks.end(), result.begin(), result.end());
  COMP(seed_b_ins);
  std::free(dSA);
}

void SNVIdentifier::buildVariantBlocks(int64_t const * dSA, int64_t const
    dSA_sz,int64_t * seed_idx, int64_t const *
    const to, vector< pair<int64_t, int64_t> > const & bsa, string const &
    concat, set< shared_ptr<bpBlock> , vBlockComp> & twork) {

  // push seed_idx back to start of group
  if (seed_idx > dSA) {
    while (computeLCP(concat.cbegin() + *seed_idx, concat.cbegin() +
          *(seed_idx-1)) >= gsa->get_min_suf_size()) {
      --seed_idx;
      if (seed_idx == dSA) break;
    }
  }
  int64_t * ext = seed_idx + 1;
  while (seed_idx < to && seed_idx != (dSA + dSA_sz - 1)) {
    while(computeLCP(concat.cbegin() + *seed_idx, concat.cbegin() + *ext) >=
        gsa->get_min_suf_size()) {
      ++ext;

      if (ext == (dSA + dSA_sz)) break;
    }
    if (ext - seed_idx >= GSA2_MCT) {
      // Then load into a variant block
      shared_ptr<bpBlock> block(new bpBlock);
      for (int64_t * it = seed_idx; it < ext; ++it) {
        read_tag rt = constructReadTag(dSA, dSA_sz, bsa, concat, it);
        block->insert(rt);
      }
      twork.insert(block);
    }
    seed_idx = ext++;
  }
  if (seed_idx == (dSA + dSA_sz - 1) && GSA2_MCT == 1) {
    shared_ptr<bpBlock> block;
    block.reset(new bpBlock);
    block->insert(constructReadTag(dSA, dSA_sz, bsa, concat, (dSA + dSA_sz -1)));
    twork.insert(block);
  }
}

read_tag SNVIdentifier::constructReadTag(int64_t const * dSA, int64_t const dSA_sz,
    vector<pair<int64_t, int64_t> > const & bsa, string const &
    concat, int64_t const *ptr) {
  auto orientation = [dSA, dSA_sz, &concat] (int64_t const *sa_pos) {
    if (sa_pos < dSA || sa_pos > (dSA+dSA_sz)) exit(1);
    string::const_iterator it = concat.cbegin() + *sa_pos;
    for (; it < concat.cend(); ++it) {
      if (*it == TERM) return RIGHT;
      else if (*it == '#') return LEFT;
    }
  };
  pair<int64_t, int64_t> readConcat = binarySearch(bsa, *ptr);
  int offset = *ptr - readConcat.second;
  int readSize = gsa->len(readConcat.first);
  bool ori = orientation(ptr);
  if (ori == LEFT) {
    offset -= readSize;
    offset = readSize - offset - gsa->get_min_suf_size() - 1;
  }
  read_tag r(readConcat.first, offset, ori, TUMOUR);
  return r;
}


void SNVIdentifier::xorSwap(int64_t *x, int64_t *y) {
  if (x != y) {
    *x ^= *y;
    *y ^= *x;
    *x ^= *y;
  }
}

int64_t SNVIdentifier::bubbleRemove(int64_t * const a, int64_t const sz, int64_t const invalid) {
  int64_t *it = a; ++it;
  int64_t *base = a;
  while (it < (a + sz)) {
    if (*base == invalid) {
      while (it < (a + sz) && *it == invalid) ++it;
      if (it >= (a + sz)) break;
      xorSwap(base, it);
    }
    ++base; ++it;
  }
  return base - a;
}

int64_t SNVIdentifier::suffixLen(int64_t const i, string const & concat) {
  if (i < 0 || i >= concat.size()) return -1;
  int64_t len = 1;
  for (string::const_iterator it = concat.cbegin() + i; it < concat.cend() &&
      *it != TERM  && *it != '#'; ++it) ++len;
  return len;
}

void SNVIdentifier::remove_short_suffixes(int64_t* &sa, int64_t &sa_sz, int64_t
    min_suffix_length, string const & concat) {
  omp_set_num_threads(N_THREADS);
  int64_t elementsPerThread = sa_sz / N_THREADS;
#pragma omp parallel for
  for (int i = 0; i < N_THREADS; i++) {
    int64_t* from = (sa + i*elementsPerThread);
    int64_t* to = nullptr;
    if (i == N_THREADS - 1) {
      to = (sa + sa_sz);
    }
    else {
      to = (sa + (i+1)*elementsPerThread);
    }
    while (from < to) {
      if (suffixLen(*from, concat) <= min_suffix_length) {
        *from = -1;
      }
      from++;
    }
  }
  sa_sz = bubbleRemove(sa, sa_sz, -1);
  sa = (int64_t*) std::realloc(sa, sa_sz*sizeof(int64_t));
}

pair<int64_t, int64_t> SNVIdentifier::binarySearch(vector< pair<int64_t,
    int64_t> > const & bsa, int64_t sa_pos) {
  int64_t r = bsa.size();
  int64_t l = 0;
  while (l < r) {
    int64_t m = (l+r)/2;
    if (sa_pos > bsa[m].second) l = m+1;
    else if (sa_pos < bsa[m].second) r = m;
    else return bsa[m];
  }
  l--;
  return bsa[l];
}



string SNVIdentifier::readTagToString(read_tag const& tag) {
  string read = gsa->get_suffix_string(tag.read_id);
  int offset = tag.offset;
  string dollar = "";
  if (tag.orientation == LEFT) {
    offset = read.size() - (tag.offset + gsa->get_min_suf_size() + 1);
    read = reverseComplementString(read);
    dollar = TERM;
  }
  return read.substr(offset) + dollar;
}

int SNVIdentifier::computeLCP(read_tag const& a, read_tag const& b) {
  string a_str = readTagToString(a);
  string b_str = readTagToString(b);
  int lcp = 0;
  while (lcp < a_str.length() &&
         lcp < b_str.length() &&
         a_str[lcp] == b_str[lcp]) lcp++;
  return lcp;
}

void SNVIdentifier::extractBlocks(vector<read_tag> const& gsa) {
  omp_set_num_threads(N_THREADS);
  int64_t elementsPerThread = gsa.size() / N_THREADS;
#pragma omp parallel for
  for (int i=0; i < N_THREADS; i++) {
    int64_t from = i*elementsPerThread;
    int64_t to;
    if (i == N_THREADS - 1) {
      to = gsa.size();
    } else {
      to  = (i+1)*elementsPerThread;
    }
    vector<shared_ptr<bpBlock> > threadWork;
    extractBlocksWorker(from, to, &gsa, threadWork);
#pragma omp critical
    SeedBlocks.insert(SeedBlocks.end(), threadWork.begin(), threadWork.end());
  }
}

void SNVIdentifier::extractBlocksWorker(int64_t seedIndex,
                                            int64_t to,
                                            vector<read_tag> const* gsa_ptr,
                                            vector<shared_ptr<bpBlock> > & localThreadStore) {
  vector<read_tag> const& aux_gsa = *gsa_ptr;
  // backup to block start
  if (seedIndex !=  0) {
    int64_t startPoint = seedIndex;
    while (computeLCP(aux_gsa[startPoint], aux_gsa[seedIndex-1]) >=
        gsa->get_min_suf_size()) {
      seedIndex--;
      if (seedIndex == 0) break;
    }
  }
  int64_t extension{seedIndex + 1};
  while (seedIndex < to && seedIndex != aux_gsa.size() - 1) {
    shared_ptr<bpBlock> block;
    while (computeLCP(aux_gsa[seedIndex], aux_gsa[extension]) >= gsa->get_min_suf_size()) {
      extension++;
      if (extension == aux_gsa.size()) break;
    }
    if (extension - seedIndex >= GSA2_MCT) {
      block.reset(new bpBlock);
      for (int64_t i = seedIndex; i < extension; i++) block->insert(aux_gsa[i]);
      localThreadStore.push_back(block);
    } 
    seedIndex = extension++;
  }
  if (seedIndex == aux_gsa.size() - 1 && GSA2_MCT == 1) {
    shared_ptr<bpBlock> block;
    block.reset(new bpBlock);
    block->insert(aux_gsa[seedIndex]);
    localThreadStore.push_back(block);
  }
}


bool SNVIdentifier::extendBlock(int64_t seed_index, 
    bpBlock &block, bool orientation, int calibration) {
  // seed_index is the index of the unique suffix_t this function
  // was called with. 
  bool success_left{false}, success_right{false};
  int left_of_seed = 0, right_of_seed = 0;

  if (seed_index != 0){
    left_of_seed = computeLCP(gsa->suffix_at(gsa->sa_element(seed_index)),
                     gsa->suffix_at(gsa->sa_element(seed_index-1)));
  }
  if(seed_index != gsa->size()-1) {
    right_of_seed = computeLCP(gsa->suffix_at(gsa->sa_element(seed_index)), 
                     gsa->suffix_at(gsa->sa_element(seed_index+1)));
  }
  if (left_of_seed >= gsa->get_min_suf_size() && seed_index > 0){
    success_left = getSuffixesFromLeft(seed_index, block, orientation, calibration);
  }
  if (right_of_seed >= gsa->get_min_suf_size() && seed_index < (gsa->size() - 1)) {
    success_right = getSuffixesFromRight(seed_index, block, orientation, calibration);
  }
  return success_left || success_right;
}

bool SNVIdentifier::getSuffixesFromLeft(int64_t seed_index,
  bpBlock &block, bool orientation, int calibration) {
  bool success = false;
  int64_t left_arrow = seed_index-1;
  // While lexicographally adjacent suffixes share the same lcp value
  // they have the same branchpoint, thus they are in the same group,
  // so add them

  while( // lcps are same AND not out of bounds AND not already in group...
      left_arrow >= 0           
      && computeLCP(gsa->suffix_at(gsa->sa_element(left_arrow)),
                                    gsa->suffix_at(gsa->sa_element(seed_index))) >=
      gsa->get_min_suf_size()) {

    // ...add read pointed to by suffix to the block
    // however now add all reads as healthy. 
    // because of the way the set performs comparison, identical reads
    // now labled healthy will be rejected
    int64_t sa_pos = gsa->sa_element(left_arrow);
    read_tag next_read;
    next_read.read_id = gsa->read_id_of_suffix(sa_pos);
    next_read.orientation = orientation;
    if (orientation == RIGHT) {
      next_read.offset = gsa->offset(sa_pos) + calibration;
    } else { // orientation == LEFT
      next_read.offset = gsa->offset(sa_pos) - calibration;
    }

    if (gsa->tissuetype(sa_pos) == HEALTHY) {
      next_read.tissue_type = HEALTHY;
    }
    else {    // we need to switch the type of TUMOUR to SWITCHED
      next_read.tissue_type = SWITCHED;
    }

    // insert tag into block
    std::pair<bpBlock::iterator, bool> 
      insertion = block.insert(next_read);
    if (insertion.second == true) success = true;
    left_arrow--;
  }
  return success;
}

bool SNVIdentifier::getSuffixesFromRight(int64_t seed_index,
    bpBlock&block, bool orientation, int calibration) {
  bool success = false;
  int64_t right_arrow = seed_index+1;
  // While lexicographically adjacent suffixes share the same lcp val
  // they have the same branchpoint, thus they are in the same group
  // so add them

  while ( // lcps are the same AND not out of bounds AND not already in group...
      right_arrow < gsa->size()         // max LCP size is one less than SA
      && computeLCP(gsa->suffix_at(gsa->sa_element(right_arrow)), 
                                    gsa->suffix_at(gsa->sa_element(seed_index))) >=
      gsa->get_min_suf_size()) {
    int64_t sa_pos = gsa->sa_element(right_arrow);

    // ...add read pointed to by suffix to the block
    read_tag next_read;
    next_read.read_id = gsa->read_id_of_suffix(sa_pos);
    next_read.orientation = orientation;

    if (orientation == RIGHT) {
      next_read.offset = gsa->offset(sa_pos) + calibration;
    } else { // orientation == LEFT
      next_read.offset = gsa->offset(sa_pos) - calibration;
    }

    if (gsa->tissuetype(sa_pos) == HEALTHY) {
      next_read.tissue_type = HEALTHY;
    }
    else {    // we need to switch the type of TUMOUR to SWITCHED
      next_read.tissue_type = SWITCHED;
    }

    std::pair<bpBlock::iterator, bool> insertion = block.insert(next_read);
    if (insertion.second == true) success = true;
    right_arrow++;
  }
  return success;
}


int SNVIdentifier::minVal(int a, int b) {
  return (a > b) ? b : a;
}

bool SNVIdentifier::lexCompare(string const & l, string const & r,  unsigned int min_lr) {
  // return true if l < r
  // min_lr avoids redundant searches 
  string::const_iterator l_iter = l.begin() + min_lr;
  string::const_iterator r_iter = r.begin() + min_lr;
  for (; (l_iter != l.end() && r_iter != r.end()); l_iter++, r_iter++) {
    // lex compare character
    if (*l_iter < *r_iter) return true;
    if (*r_iter < *l_iter) return false;
    // equiv char so move to next...
  }

  // One is prefix of other, return the prefix as higher suffix
  return (l_iter == l.end()) && (r_iter != r.end());
}

int64_t SNVIdentifier::binarySearch(string const &query) {
  // Search bounds
  int64_t right{gsa->size() - 1};   // start at non-out of bounds
  int64_t left{0};
  int64_t mid;

  // prefix lengths with query
  unsigned int min_left_right;
  unsigned int lcp_left_query;
  unsigned int lcp_right_query;

  // find minimum prefix length of left and right bounds with query
  lcp_left_query = lcp(gsa->get_suffix_string(gsa->sa_element(left)), query, 0);
  lcp_right_query = lcp(gsa->get_suffix_string(gsa->sa_element(right)), query, 0);
  min_left_right = minVal(lcp_left_query, lcp_right_query);
  while (left <= right) {

    bool left_shift{false};
    mid = (left + right) / 2;
    if(lcp(gsa->get_suffix_string(gsa->sa_element(mid)), query, min_left_right) == query.size()) {
      // 30bp stretch covered. Arrived at genomic location. Return.
      return mid; 
    }
    if(lexCompare(gsa->get_suffix_string(gsa->sa_element(mid)), query, min_left_right)) {
      // then query lexicographically lower (indexed > mid) (higher ranked
      // characters) so move left bound towards right
      left = mid+1;
      left_shift = true;
    }
    else {
      // then query is  lexicographically higher (indexed < mid due to lower
      // ranking characters) than mid. Therefore, need to move right bound
      // towards left
       right = mid-1;
    }

    // only recompute the moved bound
    if (left_shift) {
      if (left >= 0 && left < gsa->size()) {
        lcp_left_query  = lcp(gsa->get_suffix_string(gsa->sa_element(left)), query, min_left_right);
      }
    }
    else {  // must be right_shift
      if (right >= 0 && right < gsa->size()) {
        lcp_right_query = lcp(gsa->get_suffix_string(gsa->sa_element(right)), query, min_left_right);
      }
    }
    min_left_right = minVal(lcp_left_query, lcp_right_query);
  }
  return -1; // no match
}

int SNVIdentifier::lcp(string const & l, string const & r, unsigned int mlr) {
  while (mlr < l.length() && mlr < r.length() && l[mlr] == r[mlr]) {
    ++mlr;
  }
  return mlr;
}


int SNVIdentifier::convertOffset(read_tag const& tag) {
  int64_t  sz = gsa->len(tag.read_id);
  return sz - (tag.offset + gsa->get_min_suf_size() + 1);
}

int64_t SNVIdentifier::getSize() {
  return BreakPointBlocks.size();
}


// end of file

/*
void SNVIdentifier::unifyBlocks(vector<bpBlock> & seedBlocks) {
//START(SNVIdentifier_unifyBlocks);
  std::sort(
    seedBlocks.begin(), seedBlocks.end(), 
    [](bpBlock const& a, bpBlock const& b) {
      return a.size() > b.size(); // sort ascending
    }
  );

  struct blockMergeStatus{
    unsigned int id;
    unsigned int mergeTo;
    bool merged;
    blockMergeStatus(unsigned int i, unsigned int mt, bool m): id(i),
    mergeTo(mt), merged(m){}
  };

  vector<blockMergeStatus> mTable;    // stores merge info
  for (int i = 0; i < seedBlocks.size(); i++) {
    mTable.push_back(
        blockMergeStatus(i, std::numeric_limits<unsigned int>::max(), false)
    );
  }

  unsigned int min = std::numeric_limits<unsigned int>::max();
  unsigned int max = 0;
  for (bpBlock const& b : seedBlocks) {      // find global min/max
    for (bpBlock::const_iterator block_it = b.begin();
        block_it != b.end();
        block_it++) {
      if (block_it->read_id < min) {
        min = block_it->read_id;
      }
      if (block_it->read_id > max) {
        max = block_it->read_id;
      }
    }
  }

  // begin merging
  bool mergeOccured = false;
  do {
    mergeOccured = false;   // reset after loop
    vector<vector<vector<blockMergeStatus>::iterator> >  tagArray(max+1 - min);

    // construct the tag array
    for (unsigned int mIdx = 0; mIdx < seedBlocks.size(); mIdx++) {
      if (mTable[mIdx].merged) continue;
      for (read_tag const& r : seedBlocks[mIdx]) {
        vector<blockMergeStatus>::iterator tag = mTable.begin() + mIdx;
        tagArray[r.read_id - min].push_back(tag);
      }
    }

    // merge blocks
    for (vector<vector<blockMergeStatus>::iterator> & queue : tagArray) {
      if (queue.size() < 2) continue;   // nothing to merge
      unsigned int mergeTo;
      if (queue[0]->merged) {
        mergeTo = queue[0]->mergeTo;
      } else {
        mergeTo = queue[0]->id;
      }

      for (unsigned int qIdx= queue.size()-1; qIdx > 0; qIdx--) {
        if (queue[qIdx]->merged) continue;
        mergeOccured = true;
        mergeBlocks(seedBlocks[queue[0]->id], seedBlocks[queue[qIdx]->id]);
        queue[qIdx]->merged = true;   // update mTable status post-merge
        queue[qIdx]->mergeTo = mergeTo;
      }
    }
  } while (mergeOccured);

  for (blockMergeStatus const& status : mTable) {
    if (!status.merged) { // blocks that never merged were merged to
      if (seedBlocks[status.id].size() < GSA2_MCT) continue;
      BreakPointBlocks.push_back(seedBlocks[status.id]);
    }
  }
//COMP(SNVIdentifier_unifyBlocks);
}
void SNVIdentifier::mergeBlocks(bpBlock & to, bpBlock & from) {
//START(SNVIdentifier_mergeBlocks);
  bpBlock::iterator f_it = from.begin();
  bpBlock::iterator common_read;
  // Find read in common between blocks.
  for(; (common_read = to.find(*f_it)) == to.end(); f_it++);

  // The orientation associated with each read, is its aligning orientation
  // relative to its current block. At this moment in time, the final
  // orientation (which we consider the orientation relative to the reference)
  // is unknown. Therefore, if the common read is in opposite orientations
  // in the from and to block, the orientations of each read in the
  // from block needs to be switched.
  if (common_read->orientation != f_it->orientation) {
    for(bpBlock::iterator it = from.begin();
        it != from.end(); it++) {
      if (it->orientation == LEFT) {
        it->orientation = RIGHT;
      } else {
        it->orientation = LEFT;
      }
    }
  }

  // correct for symetric blocks
  // adjustment will only be correctly calculated if reads are both
  // in the same orientation
  int common_read_offset{common_read->offset}, f_it_offset{f_it->offset};
  if (common_read->orientation == LEFT) {
    //int read_size = reads->getReadByIndex(common_read->read_id, common_read->tissue_type).size();
    //  common_read_offset = (read_size - (common_read->offset + reads->getMinSuffixSize() + 1));
    common_read_offset = convertOffset(*common_read);
  }
  if (f_it->orientation == LEFT) {
    //int read_size = reads->getReadByIndex(f_it->read_id, f_it->tissue_type).size();
    //  f_it_offset = (read_size - (f_it->offset + reads->getMinSuffixSize() + 1));
    f_it_offset = convertOffset(*f_it);
  }

  // Adjust offsets to to block.
  // Orientation check performs the correct adjustment in either
  // orientation.
  int adjustment = common_read_offset - f_it_offset;
  cout << adjustment << endl;

  for (bpBlock::iterator it = from.begin();
       it != from.end();
       it++) {
    if (it->orientation == RIGHT) {
      it->offset += adjustment;
    } else {
      it->offset -= adjustment;
    }
    to.insert(*it);
  }
//COMP(SNVIdentifier_mergeBlocks);
}
//void SNVIdentifier::transformBlock(int64_t* from, 
//     int64_t* to, vector< pair<int64_t, int64_t> > *bsa, string const& concat,
//     vector<read_tag> *block) {
//
//  auto orientation = [from, to, &concat] (int64_t *sa_pos) {
//    if (sa_pos < from || sa_pos > to) exit(1);
//    string::const_iterator it = concat.cbegin() + *sa_pos;
//    for (; it < concat.cend(); ++it) {
//      if (*it == TERM) return RIGHT;
//      else if (*it == '#') return LEFT;
//    }
//  };
//
//  for (; from < to; from++) {
//    pair<int64_t, int64_t> readConcat = binarySearch(*bsa,*from);
//    int offset = *from - readConcat.second;
//    int readSize = gsa->len(readConcat.first);
//    bool ori = orientation(from);
//    if (ori == LEFT) {
//      offset -= readSize;
//    }
//    if (readSize - offset <= gsa->get_min_suf_size()) continue;
//    if (ori == LEFT) {
//      offset = readSize - offset - gsa->get_min_suf_size() - 1;
//    }
//    read_tag t;
//    t.read_id = readConcat.first;
//    t.orientation = ori;
//    t.offset = offset;
//    t.tissue_type = TUMOUR;
//    block->push_back(t);
//  }
//}

void BreakPointBlocks::outputExtractedCancerReads(std::string const& filename) {
  ofstream ofHandle(filename.c_str());
  for (unsigned int read_idx : CancerExtraction) {
    ofHandle << "(" << read_idx << ",T)" << std::endl;
  }
  ofHandle.close();
}

void BreakPointBlocks::outputFromBPB(std::string const& filename) {
  ofstream ofHandle(filename.c_str());
  for (bpBlock const& b : BreakPointBlocks) {
    for (read_tag const& tag : b.block) {
      ofHandle << "(" << tag.read_id << ","
               << ((tag.tissue_type) ? "H" : "T")
               << ")" << std::endl;
    }
  }
  ofHandle.close();
}

BreakPointBlocks::~BreakPointBlocks() {
}
string BreakPointBlocks::buildQualityString(vector< vector<int> > const& freq_matrix,
    string const& cns, bool tissue) {
  // Function steps through each position of the string and
  // determines whether a position should be masked if:
  //  -- Cancer: if the number of bases contributing to the 
  //             consensus base is < CTR then mask. Or, if the 
  //             number of bases with frequency above the error threshold
  //             (ALLELIC_ERROR_THRESH) is > 1 then mask.
  //  -- Health: if the number of bases with frequency above the error threshold
  //             (ALLELIC_ERROR_THRESH) is > 1 then mask.
  
  string q_str("");
  int m_start{0};
  for (int i=0; i < freq_matrix[0].size(); i++) {
    if (freq_matrix[0][i] != -1) {
      m_start = i;
      break;
    }
  }

  //if (cns.size() == 0) {
  //  cout << "No consensus sequence! " << endl;
  //  exit(1);
  //}
  for (int pos=0; pos < cns.size(); pos++) {

    // Determine the number of bases above the error frequency
    // by the calulation:
    //      bases / total bases = freq.
    // bases is equivalent to the number of reads with a base of type given
    // base. 
    if (tissue == HEALTHY || tissue == TUMOUR) { // || TUMOUR testing

      double total_bases = 0;                     // get total_bases
      for (int base=0; base < 4; base++) {
        total_bases += freq_matrix[base][pos + m_start];
      }
      int n_bases_above_err_freq{0};
      for(int base=0; base < 4; base++) {
        if((freq_matrix[base][pos + m_start] / total_bases) > ALLELIC_FREQ_OF_ERROR) {
          n_bases_above_err_freq++;
        }
      }
      if (n_bases_above_err_freq > 1) { // then mask
        q_str += "L";
        continue;
      }
    }

    // Unique masking logic to cancer reads
    // The supporting evidence for the chosen consensus
    // position must be above 4
    if (tissue == TUMOUR) {
      int base{0};
      switch (cns[pos]) {
        case 'A': base = 0; break;
        case 'T': base = 1; break;
        case 'C': base = 2; break;
        case 'G': base = 3; break;
      }
      if (freq_matrix[base][pos + m_start] < GSA2_MCT) {
        q_str += "L";
        continue;
      }
    }

    q_str += "-";
  }
  return q_str;
}
   // SERIAL BUILD CNS ROUTINE
  //cout <<  "Number of cns pairs: " << consensus_pairs.size() << endl;
  //for (bpBlock &block : SeedBlocks) {
  //  consensus_pair pair;
  //  pair.left_ohang = pair.right_ohang = 0;
  //  generateConsensusSequence(TUMOUR, block, pair.mut_offset, pair.pair_id, pair.mutated, pair.mqual);

  //  if (block.block.size() > COVERAGE_UPPER_THRESHOLD) {
  //    block.block.clear();
  //    continue;
  //  }
  //  extractNonMutatedAlleles(block, pair);
  //  generateConsensusSequence(HEALTHY, block, pair.nmut_offset, pair.pair_id, pair.non_mutated, pair.nqual);
  //  if (block.block.size() > COVERAGE_UPPER_THRESHOLD) {
  //    block.clear();
  //    continue;
  //  }
  //  block.clear();
  //  trimHealthyConsensus(pair); // MUST trim healthy first
  //  trimCancerConsensus(pair);
  //  bool low_quality_block {false};
  //  maskLowQualityPositions(pair, low_quality_block);
  //  if (low_quality_block) {
  //    continue;
  //  }

  //  if (pair.mutated.empty() || pair.non_mutated.empty()) {
  //    continue;
  //  }
  //  consensus_pairs.push_back(pair);
  //}
  /// cout << "Pair id: " << pair.pair_id << endl;
  /// cout << "Tumour: " << endl;
  /// cout << pair.mutated << endl;
  /// cout << pair.mqual << endl;
  /// cout << "Block mutated offset: " << pair.mut_offset << endl;
  /// cout << "Healthy: " << endl;
  /// cout << pair.non_mutated << endl;
  /// cout << pair.nqual << endl;
  /// cout << "Block non_mut offset: " << pair.nmut_offset << endl;
  //} 
  //Cout << "n skipped: " << skipped << endl;
  //Cout << "DONE BUILDING PAIRS" << endl;
  //extractNonMutatedAlleles();
  //outputFromBPB("/data/ic711/point5.txt");


//Old version, does not correctly take into account CTR, and 
// also suffers from rare case bug
void BreakPointBlocks::extractionWorker(unsigned int seed_index, 
                                                        unsigned int to) {

  set<unsigned int> localThreadExtraction;
  double cancer_sequences = 0, healthy_sequences = 0;
  unsigned int extension;


  while (seed_index < to-1) {
   
    // Calc econt of suffixes sharing same genomic location (assumption: 
    // LCP >= 30). If econt below user specified value, extract group
    if(::computeLCP(SA->getElem(seed_index),
                                  SA->getElem(seed_index+1), *reads) >= 30) {

      // tally up the tissue types of seed suffix
      if (SA->getElem(seed_index).type == HEALTHY) {
        healthy_sequences++;
      } else {  // TUMOUR
        cancer_sequences++;
      }

      extension = seed_index+1;
      while ((::computeLCP(SA->getElem(seed_index), SA->getElem(extension), *reads) >= 30)
            ) {

        if (SA->getElem(extension).type == HEALTHY) {
          healthy_sequences++;
        } else { // TUMOUR 
          cancer_sequences++;
        }
        extension++;
        if(extension == SA->getSize()) {break;}
      }

      // check below econt
      if (cancer_sequences >= GSA1_MCT   && 
         ((healthy_sequences / cancer_sequences) <= econt )) { // permit group
        for (unsigned int i = seed_index; i < extension; i++) {
          if(SA->getElem(i).type == TUMOUR) {     // only extr. cancer reads
            localThreadExtraction.insert(SA->getElem(i).read_id);
            //{
            //  std::lock_guard<std::mutex> cout_guard(cout_lock);
            //  cout << reads->returnSuffix(SA->getElem(i)) << endl;
            //}
          }
        }
      }
      seed_index = extension; // jump scan to end of group
    }

    else { // if CTR == 1 {
      // In this case, there is no > 2 group of reads containing 
      // this 30bp stretch of DNA. However, the read is still unique
      // it could be unique due to a mutation, and a group could not
      // form due to MTT, mutation distribution or low coverage 
      // therefore, if it is a cancer read, extract it
      //if (SA->getElem(seed_index).type == TUMOUR) {
      //  localThreadExtraction.insert(seed_index);
      //}
      seed_index++;
    }

    healthy_sequences = 0; 
    cancer_sequences = 0; 

  }//end while

  // After extraction, threads need to load their data into global CancerExtractionSet


  std::lock_guard<std::mutex> lock(cancer_extraction_lock); // avoid thread interference
  for(unsigned int read_with_mutation : localThreadExtraction) {
      CancerExtraction.insert(read_with_mutation);
    }
}
void BreakPointBlocks::seedBreakPointBlocks() {

  string concat("");
  vector<pair<unsigned int, unsigned int>> binary_search_array;

  // load the first element
  set<unsigned int>::iterator it = CancerExtraction.begin();
  pair<unsigned int, unsigned int> first_read(*it,0);
  binary_search_array.push_back(first_read);    // first read starts at zero
  concat += reads->getReadByIndex(*it, TUMOUR);
  concat += reverseComplementString( 
      reads->getReadByIndex(*it, TUMOUR)
  ) + "$";

  cout << "Cancer Extraction size: " <<  CancerExtraction.size() << endl;
  unsigned int concat_idx = 0;
  it++;
  for (; it != CancerExtraction.end(); it++) {
    // build concat
    concat += reads->getReadByIndex(*it, TUMOUR);
    concat += reverseComplementString(reads->getReadByIndex(*it, TUMOUR)) + "$";
    // add bsa values
    concat_idx += reads->getReadByIndex(*std::prev(it), TUMOUR).size() * 2;
    pair<unsigned int, unsigned int> read_concat_pair (*it, concat_idx);
    binary_search_array.push_back(read_concat_pair);
  }

  // Build SA
  cout << "Building cancer specific sa" << endl;
  unsigned long long *radixSA = 
    Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();
  cout << "Size of cancer specific sa: " << concat.size() << endl;

  vector<read_tag> gsa;
  for (unsigned long long i=0; i < concat.size(); i++) {
    pair <unsigned int, unsigned int> read_concat_pair = 
      SA->binarySearch(binary_search_array, radixSA[i]);

    unsigned int offset = radixSA[i] - read_concat_pair.second;
    int read_size = reads->getReadByIndex(read_concat_pair.first, TUMOUR).size();

    bool orientation = RIGHT;
    if (offset >= read_size) {
      orientation = LEFT;
      offset -= read_size;
    }
    // remove suffixes that are too short
    if (read_size - offset <= reads->getMinSuffixSize())  continue;

    // read is stored in forward orientation, convert if LEFT
    if (orientation == LEFT) {
      offset = read_size - offset - reads->getMinSuffixSize() - 1; // dollar -1
    }

    read_tag tag;
    tag.read_id = read_concat_pair.first;
    tag.orientation = orientation;
    tag.offset = offset;
    tag.tissue_type = TUMOUR;
    gsa.push_back(tag);       // gsa should be built
  } 

  delete [] radixSA;
  //// PRINT gsa
  //for (read_tag const& tag : gsa) {
  //  std::cout << readTagToString(tag) 
  //            << ((tag.orientation) ? " -- R" : " -- L") << endl;
  //}

  extractGroups(gsa);
}

void BreakPointBlocks::maskLowQualityPositions(consensus_pair & pair, bool &
    low_quality) {
  int low_quality_count{0};
  for (int pos=0; pos < pair.mutated.size(); pos++) {
    if (pair.mqual[pos] != '-' || pair.nqual[pos + pair.left_ohang] != '-') {
      pair.mutated[pos] = pair.non_mutated[pos + pair.left_ohang];
    }
    if (pair.mqual[pos] == 'L' || pair.nqual[pos] == 'L') {
      low_quality_count++;
    }
  }
  if (low_quality_count > MAX_LOW_CONFIDENCE_POS) low_quality = true;
}
char BreakPointBlocks::revCompCharacter(char ch, bool rc) {
  if (!rc) return ch;
  switch (ch) {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'G': return 'C';
    case 'C': return 'G';

    default: {cout << "Hit Dollar" << endl; exit(1);}
  }
}

////// V2: Takes into account CTR and does not have rare case bug
//void BreakPointBlocks::extractGroups(vector<read_tag> const& gsa) {
//  // Adding from, to parameters, in order to make logic compatible
//  // with later multithreading modifications.
//  unsigned int to{gsa.size()}, seed_index{0};
//  unsigned int extension{seed_index + 1};
//  unsigned int block_id;
//
//  // < to is thread safe, != gsa.size() - 1 is rare case bound safe
//  while (seed_index < to && seed_index != gsa.size() - 1) {
//    bpBlock block;
//    // compute group size - avoid inserting here to avoid unnecessary mallocs
//    while (computeLCP(gsa[seed_index], gsa[extension]) >= reads->getMinSuffixSize()) {
//      extension++;
//      if (extension == gsa.size()) break;
//    }
//    // Make alloc if group at or above CTR
//    if (extension - seed_index >= GSA2_MCT) {
//      for (int i=seed_index; i < extension; i++) block.block.insert(gsa[i]); 
//    }
//    else {    // continue, discarding group
//      seed_index = extension++;
//      continue;
//    }
//  
//    // Load group into break point blocks
//    block.id = block_id;
//    //BlockSeeds.push_back(block);     // for loss of sensitivity by way of 2
//    SeedBlocks.push_back(block);
//    block_id++;
//    seed_index = extension++;
//  }
//
//  if (seed_index == gsa.size() - 1 && GSA2_MCT == 1) {
//    bpBlock block;
//    block.block.insert(gsa[seed_index]);
//    block.id = block_id;
//    //BlockSeeds.push_back(block);      // for loss of sensitivity by way of 2
//    SeedBlocks.push_back(block);
//  }
//}

//// V1: Does not take into account CTR, and rare case bug
//void BreakPointBlocks::extractGroups(vector<read_tag> &gsa) {
//
//  unsigned int seed_index{0};
//  unsigned int extension{0};
//  unsigned int block_id{0};
//  while (seed_index < gsa.size()-1) {
//   
//    // seed
//    if (computeLCP(gsa[seed_index], gsa[seed_index+1]) >= reads->getMinSuffixSize()) {
//      bpBlock block;
//      block.block.insert(gsa[seed_index]); // add seed to block
//      extension = seed_index+1;
//      // extend
//      while ((computeLCP(gsa[seed_index], gsa[extension])
//           >= reads->getMinSuffixSize())) {
//        block.block.insert(gsa[extension]);
//        extension++;
//        if (extension == gsa.size()) break;
//      }
//      if (block.block.size() < 4)  {
//        seed_index = extension;
//        continue;
//      }
//
//
//      block.id = block_id;
//      BreakPointBlocks.push_back(block);
//      block_id++;
//
//      seed_index = extension;
//    }
//    else {
//      seed_index++;
//    }
//  }
//  
//  // PRINT groups
//  for (bpBlock const& b : BreakPointBlocks) {
//    cout << "Block id: " << b.id << endl;
//    for (read_tag const& tag : b.block) {
//      cout << readTagToString(tag) << endl;
//    }
//  }
//}



void BreakPointBlocks::extractNonMutatedAlleles() {
//// NB: When we perform integer grouping, this will have be be performed
//// in forward and reverse orientation

  for (bpBlock &block : BreakPointBlocks) {
    read_tag tag = *block.block.begin();  // copy first element

    string read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
    read = read.substr(tag.offset, reads->getMinSuffixSize());
    string rev_read = reverseComplementString(read);

    long long int fwd_index = binarySearch(read);
    long long int rev_index = binarySearch(rev_read);

    if (fwd_index != -1) {
      extendBlock(fwd_index, block.block, tag.orientation, 0);
    }
    else {
      //std::cout << "Search failed" << std::endl;
    }

    if (rev_index != -1) {
      extendBlock(rev_index, block.block, !tag.orientation, 0);
    }
    else {
      //std::cout << "Search failed" << std::endl;
    }
  }
}







//void BreakPointBlocks::extractNonMutatedAlleles(bpBlock &block, consensus_pair
//    &pair) {
//
//  bool LEFT{false}, RIGHT{true};
//  // For each consensus sequence, where cns length permits, search three 30bp windows
//  // from indexes offset-30, offset, offset+30, assuming minimum suffix size
//  // is 30.
//  for (int i = pair.mut_offset - reads->getMinSuffixSize();
//       i <= pair.mut_offset + reads->getMinSuffixSize();
//       i += reads->getMinSuffixSize()) {
//    // substring in cns bounds
//    if (i < 0 || i > pair.mutated.size() - reads->getMinSuffixSize()) {
//      continue;
//    }
//
//    //START(SNVIdentifier_queryBuild);
//    string query = pair.mutated.substr(i, reads->getMinSuffixSize());
//    string rcquery = reverseComplementString(query);
//    END(queryBuild);
//    TIME(queryBuild);
//    PRINT(queryBuild);
//
//    //START(SNVIdentifier_binSearch);
//    long long int fwd_idx = binarySearch(query);
//    long long int rev_idx = binarySearch(rcquery);
//    END(binSearch);
//    TIME(binSearch);
//    PRINT(binSearch);
//
//    //START(SNVIdentifier_extendBlock);
//    if (fwd_idx != -1) {
//      extendBlock(fwd_idx, block.block, RIGHT, pair.mut_offset - i);
//      //for (read_tag const& tag : fwd_result) {
//      //  tag.offset += pair.mut_offset - i; // + (+30), + 0, + (-30)
//      //}
//    }
//    if (rev_idx != -1) {
//      extendBlock(rev_idx, block.block, LEFT, pair.mut_offset - i);
//      //for (read_tag const& tag : rev_result) {
//      //  tag.offset -= pair.mut_offset - i; // - (+30), - 0, - (-30)
//      //}
//    }
//    END(extendBlock);
//    TIME(extendBlock);
//    PRINT(extendBlock);
//  }
////  long long int fwd_idx = binarySearch(query);
////  long long int rev_idx = binarySearch(rcquery);
////  if (fwd_idx != -1) {
////    extendBlock(fwd_idx, block.block, RIGHT);
////  }
////  else {
////    cout << "forward search failed" << endl;
////  }
////  if (rev_idx != -1) {
////    extendBlock(rev_idx, block.block, LEFT);
////  }
////  else {
////    cout << "reverse search failed" << endl;
////  }
////
////  // if search failed, then re-search downstream and upstream 30bp
////  if (fwd_idx == -1 && rev_idx == -1) {
////    fail = true;
////    if (pair.mut_offset >= reads->getMinSuffixSize()) { // then enough dwnstr seq
////      int dwnstr = pair.mut_offset - reads->getMinSuffixSize();
////      query = pair.mutated.substr(dwnstr, reads->getMinSuffixSize());
////      rcquery = reverseComplementString(query);
////      fwd_idx = binarySearch(query);
////      rev_idx = binarySearch(rcquery);
////      set<read_tag, read_tag_compare> fwd_result, rev_result;
////      if (fwd_idx != -1) {
////        extendBlock(fwd_idx, fwd_result, RIGHT);
////      }
////      if (rev_idx != -1) {
////        extendBlock(rev_idx, rev_result, LEFT);
////      }
////      for (read_tag const& tag : fwd_result) {
////        tag.offset += reads->getMinSuffixSize();
////      }
////      for (read_tag const& tag : rev_result) {
////        tag.offset -= reads->getMinSuffixSize();
////      }
////      block.block.insert(fwd_result.begin(), fwd_result.end());
////      block.block.insert(rev_result.begin(), rev_result.end());
////    }
////
////    if (pair.mutated.size() - (pair.mut_offset + reads->getMinSuffixSize()) - 1 >=
////        reads->getMinSuffixSize()) { // then enough upstr seq
////      int upstr = pair.mut_offset + reads->getMinSuffixSize();
////      query = pair.mutated.substr(upstr, reads->getMinSuffixSize());
////      rcquery = reverseComplementString(query);
////      fwd_idx = binarySearch(query);
////      rev_idx = binarySearch(rcquery);
////      set<read_tag, read_tag_compare> fwd_result, rev_result;
////      if (fwd_idx != -1) {
////        extendBlock(fwd_idx, fwd_result, RIGHT);
////      }
////      if (rev_idx != -1) {
////        extendBlock(rev_idx, rev_result, LEFT);
////      }
////      for (read_tag const& tag : fwd_result) {
////        tag.offset -= reads->getMinSuffixSize();
////      }
////      for (read_tag const& tag : rev_result) {
////        tag.offset += reads->getMinSuffixSize();
////      }
////      block.block.insert(fwd_result.begin(), fwd_result.end());
////      block.block.insert(rev_result.begin(), rev_result.end());
////    }
////  }
//}
bool BreakPointBlocks::generateConsensusSequence(unsigned int block_idx, 
    int &cns_offset, bool tissue_type, unsigned int &pair_id, 
    string & cns, string & qual) {


  // DEBUG BOOL

  // all seqs get converted to RIGHT orientation, before consensus

  if (BreakPointBlocks[block_idx].size() > COVERAGE_UPPER_THRESHOLD) {
    return true;
  } 


  // select only one tissue type
  vector<read_tag> type_subset;
  for(read_tag tag : BreakPointBlocks[block_idx].block) {
    if (tissue_type == HEALTHY  && 
       (tag.tissue_type == HEALTHY || tag.tissue_type == SWITCHED)) {
      type_subset.push_back(tag);
    }
    else if(tissue_type == TUMOUR && tag.tissue_type == TUMOUR) {
      type_subset.push_back(tag);
    }
  }

  pair_id = BreakPointBlocks[block_idx].id;

  if(type_subset.size() == 0) {   // no seq. of tissue type, cannot be mapped
    return true;
   // DEBUG_BOOL = true;
  }

  // perform offset conversion, converting LEFT offsets to RIGHT
  int max_offset = 0;
  int min_offset = numeric_limits<int>::max();

  // LEFT reads are currently stored as RIGHT, with RIGHT offset
  // Need to convert offset to LEFT for correct alignment 
  // with the rest of the group
  for (read_tag &tag : type_subset) {
    if (tag.orientation == LEFT) {
      int read_size = reads->getReadByIndex(tag.read_id,
          tag.tissue_type).size();
      tag.offset = (read_size - (tag.offset + reads->getMinSuffixSize() + 1));
    }
  //for (read_tag &tag : type_subset) {
    if(tag.offset > max_offset) {   // also record max and min offset
      max_offset = tag.offset;
    }
    if(tag.offset < min_offset) {
      min_offset = tag.offset;
    }
  }

  // value did not change and will cause std::bad_alloc, so set to 0
  if(min_offset == numeric_limits<int>::max()) {
    min_offset = 0;
  }

  // decl alignment block
  vector<string> aligned_block;
  // initialize align_counter.
  vector< vector<int> > align_counter;
  for (int n_vectors=0; n_vectors < 4; n_vectors++) {
    vector<int> v(max_offset + reads->maxLengthRead() - min_offset, 0);
    align_counter.push_back(v);
  }


  for(read_tag tag : type_subset) {

    // LEFT reads are stored as RIGHT, need to convert sequence
    // to LEFT (reverse complement) as this is the orientation
    // in which the read contained a 30bp overlap with the
    // group
    string read;
    if(tag.orientation == LEFT) {
     read = reverseComplementString( 
         reads->getReadByIndex(tag.read_id, tag.tissue_type)
         );
    }
    else {
      read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
      read.pop_back();    // remove dollar symbol
    }

    for(int i=0; i < read.size(); i++) {
      switch(read[i]) {
        case 'A':
          align_counter[0][(max_offset - tag.offset) + i]++;
          break;

        case 'T':
          align_counter[1][(max_offset - tag.offset) + i]++;
          break;

        case 'C':
          align_counter[2][(max_offset - tag.offset) + i]++;
          break;

        case 'G':
          align_counter[3][(max_offset - tag.offset) + i]++;
          break;
      }
    }
    aligned_block.push_back(addGaps(max_offset - tag.offset) + read);
  }


  cns = ""; // seed empty consensus sequence

  int n_skipped_start_pos=0;
  bool hit_consensus_min = false;
  for (int pos=0; pos < align_counter[0].size(); pos++) {

    // only start when the number of reads aligned at that position
    /// is >= GSA2_MCT 

    int nreads=0;
    for (int base=0; base < 4; base++) {
      nreads += align_counter[base][pos];
    }
    if (nreads < GSA2_MCT) {
      if (!hit_consensus_min) {
        n_skipped_start_pos++;		// need to nibble off skipped positions from cns offset
      }
      invalidatePosition(align_counter, pos);   // invalidate rather than trim
      continue;                   // continue performs the trimming
    }
    else {
      hit_consensus_min = true;
    }


    int maxVal = 0, maxInd = 0;

    // find highest freq. base
    for(int base=0; base < 4; base++) {
      if (align_counter[base][pos] > maxVal) {
        maxVal = align_counter[base][pos];
        maxInd = base;
      }
    }

    switch(maxInd) {
      case 0:
        cns += "A";
        break;
      case 1:
        cns += "T";
        break;
      case 2:
        cns += "C";
        break;
      case 3:
        cns += "G";
        break;
    }
  }// end for
  
  // if empty string, return that we need to skip block
  if (cns == "") {
    return true;
  }
  //if (cns == "" && type_subset.size() < 4) DEBUG_BOOL = true;


  //////// DEBUG
  //std::cout << ((tissue_type) ? "Healthy" : "Cancer") << " sub-block below" <<
  //std::endl;
  //cout << "Block id: " << BreakPointBlocks[block_idx].id  << endl;
  //for (int i=0; i < aligned_block.size(); i++) { // SHOW ALIGNED BLOCK
  //  cout << aligned_block[i]  << ((type_subset[i].orientation == RIGHT) ? ", R" : ", L") 
  //       << ((type_subset[i].tissue_type % 2) ? ", (H," : 
  //           ((type_subset[i].tissue_type == SWITCHED) ? ", (S," : ", (T, ridx: ")) 
  //       << type_subset[i].read_id << ")" << endl;
  //}
  //cout << "CONSENSUS AND CNS LEN" <<  cns.size() << endl;
  //cout << cns << endl << endl;
  //cout << "QSTRING" << endl;
  //cout << buildQualityString(align_counter, cns,  tissue_type) << endl;


  ////for(vector<int> v : align_counter) {
  ////  for(int i : v) {
  ////    cout << i << ",";
  ////  }
  ////  cout << endl;
  ////}
  ////cout << endl;
  

  // set the offset value
  cns_offset = (max_offset - n_skipped_start_pos);

  // set the quality string
//  qual = buildQualityString(align_counter, cns, tissue_type);
  return false;
  //return DEBUG_BOOL;
}
void BreakPointBlocks::invalidatePosition(vector< vector<int> > &alignment_counter, int pos) {
  for (int base=0; base < 4; base++) {
    alignment_counter[base][pos] = -1;
  }
}
long long int BreakPointBlocks::backUpToFirstMatch(long long int bs_hit, string query) {
  while (bs_hit >= 0) {
    if (lcp(reads->returnSuffix(SA->getElem(bs_hit)), query, 0) != query.size()){
      return bs_hit+1;
    }
    else {
      bs_hit--;
    }
  }
  return bs_hit + 1; // 0
}

void BreakPointBlocks::printAlignedBlocks() {
  for (bpBlock &block : BreakPointBlocks) {
    printAlignedBlock(block);
  }
}
void BreakPointBlocks::buildConsensusPairs() {
  for (bpBlock block : BreakPointBlocks) {
    cns_pair pair = buildConsensusPair(block);
    gmp->consensus_pairs.push_back(pair);
  }
}

void BreakPointBlocks::buildConsensusPair(bpBlock & block) {
  consensus_pair pair;
  generateConsensus(pair.ccns, pair.cqual, block, TUMOUR);
  findHealthyReads(pair.ccns, pair.cqual, block);
  generateConsnsus(pair.hcns, pair.hqual, block, HEALTHY);
  trimCancerConsensus(pair);
  maskLowQualityPositions(pair);
  gmp->consensus_pairs.push_back(pair);
}
void BreakPointBlocks::printAlignedBlocks() {
  for (bpBlock const & block : BreakPointBlocks) {
  int max = std::numeric_limits<int>::min();
    for (read_tag tag : block.block) {
      if (tag.offset > max) {

        max = tag.offset;
      } 
    }

    cout << "block id: " << block.id << endl;
    cout << "block size " << block.block.size() << endl;

    for (read_tag tag : block.block) {
      string read;
      if (tag.orientation == LEFT) {
        int read_size = reads->getReadByIndex(tag.read_id,
            tag.tissue_type).size();
        tag.offset = (read_size - (tag.offset + reads->getMinSuffixSize() + 1));
        read = reverseComplementString(reads->getReadByIndex(tag.read_id,
              tag.tissue_type));
      } else {
        read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
        read.pop_back();
      }
      cout << addGaps(max - tag.offset) << read << " - " << tag.offset << 
           " - " << tag.orientation << endl;
    }
  }
}
void BreakPointBlocks::printBreakPointBlocks() {
  for(bpBlock block : BreakPointBlocks) {
    cout << "Block id: " << block.id << endl;
    cout << "block size" << block.size() << endl;
    for (read_tag tag: block.block) {
      cout << ((tag.tissue_type) ? "Healthy" : "Cancer") 
        << endl << "read_id: " << tag.read_id << endl
        << "offset: " << tag.offset << endl
        << ((tag.orientation) ? "RIGHT" : "LEFT" ) << endl;
      cout << reads->getReadByIndex(tag.read_id, tag.tissue_type) << endl;
    }

    cout << "End of block: -----------------------------------------" << endl;
  }
}
void BreakPointBlocks::unifyBreakPointBlocks(){ O(N^2)
  pair<bool, unsigned int> test(false, 100);
  pair<bool, unsigned int> test1(false, 54);
  pair<bool, unsigned int> test2(false, 88);

  (*BreakpointBlocks)[10].insert(test);
  (*BreakpointBlocks)[14].insert(test1);
  (*BreakpointBlocks)[25].insert(test2);
  for(set<pair<bool, unsigned int>> mutation_set : *BreakpointBlocks){
    for(pair<bool, unsigned int> read : mutation_set) {
      cout << read.first << " --- " << read.second  << endl;
    }
  }


  for(int i=0; i < BreakpointBlocks->size(); i++) { // compare each block
    for(int j=i; j < BreakpointBlocks->size(); j++) { // agains all other blocks

      if(i == j){ continue;}    // skip self comparisons

      std::set<pair<bool, unsigned int>>::iterator i_iter, j_iter, i_end, j_end;

      i_iter = (*BreakpointBlocks)[i].begin();    // point it to start of set
      i_end  = (*BreakpointBlocks)[i].end();

      j_iter = (*BreakpointBlocks)[j].begin();
      j_end  = (*BreakpointBlocks)[j].end();

      // quick compare.. if largest element is smaller than smallest in other
      // no overlap
      i_end--; j_end--;
      if( (*(i_end) < *j_iter)   // no elements can be same
          || 
          (*(j_end) < *i_iter)   // again no elems can be same
        ) { continue;}

      i_end++; j_end++;

      // make pass through set, comparing elems
      while ( (i_iter != i_end) && (j_iter != j_end) ) {

        if (*i_iter > *j_iter) {
          j_iter++;
        }
        else if (*j_iter > *i_iter) {
          i_iter++;
        }
        else { // *i_iter == *j_iter so unify groups

          // unify to i
          std::set<pair<bool, unsigned int> >::iterator unify;
          unify = (*BreakpointBlocks)[j].begin();
          for (;unify != (*BreakpointBlocks)[j].end(); unify++) {
            (*BreakpointBlocks)[i].insert(*unify);  // add unique elems j to i
          }

          // unify to j
          unify = (*BreakpointBlocks)[i].begin();
          for (;unify != (*BreakpointBlocks)[i].end(); unify++) {
            (*BreakpointBlocks)[j].insert(*unify);  // add unique elems i to j
          }
        }

        i_iter++; j_iter++;

      } // end while


    } // end j for 
  }   // end i for

}
void BreakPointBlocks::extractNonMutatedAllelesDeep() {
  for(set<read_tag, read_tag_compare> &block : BreakPointBlocks) {

    for(read_tag tag : block) {
      string read = reads->getReadByIndex(tag.read_id, tag.tissue_type);
      for(int i=0; i < reads->getReadByIndex(tag.read_id,
            tag.tissue_type).size()-30; i++) {

        string aligned_section = read.substr(i, 30);

        long long int index = binarySearch(aligned_section);
        if(index == -1) {
          cout << "Error, a previously identified sequence no longer exists" <<
            " in suffix array... program terminating. "  << endl;
          exit(1);
        }
        else {
          extendBlock(index, block, tag.orientation);
        }
      }
    }
  }
}









bool hitDollar(Suffix_t &lhs, Suffix_t &rhs) {

  // Already know that the 30 positions match, so start computation from 
  // 29th position in string
  string::iterator lhs_string_iter = returnStartIterator(&lhs, *reads) + 29;
  string::iterator rhs_string_iter = returnStartIterator(&rhs, *reads) + 29;

  string::iterator lhs_end =         returnEndIterator(&lhs, *reads); 
  string::iterator rhs_end =         returnEndIterator(&rhs, *reads);

  // extend further and compare
  while (*lhs_string_iter  == *rhs_string_iter) {
    lhs_string_iter++; rhs_string_iter++;
  }

  // hit a mismatch, check if the mismatch is a dollar symbol
  if(*lhs_string_iter == '$' || *rhs_string_iter == '$') {
    return true;
  }

  // else
  return false;
}




void BreakPointBlocks::mergeHash() {

  struct hashtag{
    unsigned long long int hash_id;
    unsigned int group;
    unsigned int read_id;
    bool orienation;
  };

  unsigned long long int gmin, gmax;
  gmin = std::numerical_limits<unsigned long long int>::max();
  gmax = 0;

  std::hash<std::string> h;
  vector<vector<hashtag>> hash_list;
  hashtag next_forward_hash;
  hashtag next_reverse_hash;
  for(int i=0; i < ComplementaryUnified->size(); i++){         // M
    set<unsigned int> group = (*ComplementaryUnified)[i];    
    vector<hashtag> group_hashes;

    for(unsigned int read : group) {                          // M*C

      string read_string = reads->getReadByIndex(read, TUMOUR);  
      for(int i=0; i < read_string.size()-30; i++) {            // M*C*n
        next_forward_hash.hash_id = h(read_string.substr(i, 30));
        next_forward_hash.group = i;
        next_forward_hash.orientation = true;
        next_forward_hash.read_id = read;
      }

      string reverse_read = reverseComplementString(read_string);
      for(int i=0; i < reverse_read.size()-30; i++) {
        next_reverse_hash.hash_id = h(reverse_read.substr(i, 30));
        next_reverse_hash.group = i;
        next_reverse_hash.orientation = false;
        next_reverse_hash.read_id = read;
      }

      group_hashes.push_back(next_forward_hash);
      group_hashes.push_back(next_reverse_hash);
    }

    hash_list.push_back(group_hashes);
  }


  // find global min max
  for(hashtag hashed_substring : hash_list) {      // M.C.n
    if(hashed_substring.hash_id < gmin) {
      gmin = hashed_substring.hash_id;
    }
    if(hashed_substring.hash_id > gmax) {
      gmax = hashed_substring.hash_id;
    }
  }

  struct group_unif_info {
    unsigned int group;
    unsigned int prev_unification_group;
    bool unified;
  };


  // make unification table
  vector<group_unif_info> unif_table;
  unif_table.reserve(CompelmentaryUnified->size());
  for(unsigned int i=0; i < ComplemenratyUnified->size(); i++) {
    group_unif_unfo info;
    info.group = i;
    info.unified = false;
    unif_table.push_back(info);
  }


  // Make tagArray

  typedef pair<vector<group_unif_info>::iterator, vector<hashtag>::iterator> 
  lut_hash_pointer;

  lut_hash_pointer next_ptr;
  vector<vector<lut_hash_pointer>> tagArray;    // tagArrau

  for(unsigned int i=0; i < ComplemenratyUnified->size(); i++) {
    for(int j=0; j < hash_list[i].size(); j++) {
      lut_hash_pointer.first = unif_table.begin() + i;
      lut_hash_pointer.second = hash_list[i].begin() + j;

      tagArray[(*lut_hash_pointer.second).hash_id] = lut_hash_pointer;
    }
  }



  // begin unif

  for(int i; i < (gmax - gmin); i++) {

    if(tagArray[i].size() ==  1) {continue;}

    bool forward = false, reverse = false;
    for(lut_hash_pointer tag: tagArray[i]) {
      if(tag.second->orientation = true) {
        forward = true;
      }
      if(tag.second->orientation = false) {
        reverse = true;
      }
    }

    // group not single orientation
    if (forward == false ||  reverse == false) {continue;}

    // all reads in should share 30bp in common, and if no collisions
    // counter should thus equal size()-1
    vector<unsigned int> collisions;

    // for every tag in hashpoint, try to merge to base
    for(int j=0; j< tagArray[i].size(); j++) {

        if(tagArray[i][0].second->orientation == 
           tagArray[i][j].second->orientation) {counter++; continue;}

        // collision check, i.e validate with actuall slow test
        // should be an exact match!! so change for seq matcher
        if(thirtyBasePairOverlap( tagArray[i][0].second->read_id, 
                                  tagArray[i][j].second->read_id) ) {
          // merge

          // if the base group has already been unified with group x, then
          // unify all in list, with group x
          int unify_upon_this_group;


          // check unif_table to see if the base group has been unified...
          if(tagArray[i][0].first->unified) {  // if so, unify to that group
            unify_upon_this_group = tagArray[i][0].first->prev_unification_group;  
          }
          else {  // if not, unify to the base
            unify_upon_this_group = tagArray[i][0].first->group; 
          }

          // begin merging
        
          // point to group to merge
          set<unsigned int>::iterator j_it = 
          (*ComplementaryUnfied)[tagArray[i][j].first->group].begin();

          set<unsigned int::iterator j_end = 
          (*ComplementaryUnfied)[tagArray[i][j].first->group].end();

          // merge with base group
          while (j_it != j_end) {
            (*CompelementaryUnfied)[unify_upon_this_group]->insert(*j_it);  
          }

        }
        else {
          collisions.push_back(j);
        }
      }
      for(int collision : collisions) {
        cout << "Collision with: " << tagArray[i][collision].second->read_id <<
        endl;
      }
    
  }
}
void BreakPointBlocks::discardSingleOrientationBlocks() {

  bool right_pair = false, left_pair = false;

  vector<set<pair<bool, unsigned int>>>::iterator block = BreakpointBlocks->begin();


  while(block != BreakpointBlocks->end()) {


    for (pair<bool, unsigned int> read : *block) {

      // check if the block contains two different read orientations
      if ((*reads->TumourMateOrder)[read.second]) {
        right_pair = true;
      }
      else {
        left_pair = true;
      }
    }

    // delete if group has one orientation
    if (right_pair == false || left_pair == false) {
      BreakpointBlocks->erase(block);
    }
    else {
      block++;
    }

    right_pair = false; left_pair = false;    // reset

  }
}
void BreakPointBlocks::makeHashGroups() {
  enum {LEFT, RIGHT};      // mate pair enum

  struct hashtag {
    unsigned int hash_id;
    unsigned int read_id;
    unsigned int offset;
    bool orientation;
    bool inblock;
  };


  // Make all the hashes...

  // this will contain all 100 hashes from each extracted read and its rev comp
  vector<hashtag> hash_list;
  hash_list.reserve(CancerExtraction->size() * 30 * 100);

  std::hash<std::string> h; // the glorious h() function

  // the forward and reverse hashes for each string
  hashtag fwd_hash;
  hashtag rev_hash;
  
  // for each read from each extracted group..
  for(set<unsigned int> mutation_site: *CancerExtraction) {
    for (unsigned int read_index : mutation_site) {  
    // .. hash all reads and rev comps by sliding a 30pb window

      string read = reads->getReadByIndex(read_index, TUMOUR);
      string rev_read = reverseComplementString(read);

      for(int i=0; i < read.size()-30; i++) {

        // iterate a hash over a 30bp window
        fwd_hash.hash_id = h(read.substr(i, 30)); 
        rev_hash.hash_id = h(rev_read.substr(i, 30));

        fwd_hash.read_id = read_index;
        rev_hash.read_id = read_index;

        fwd_hash.offset = i;
        rev_hash.offset = i;

        fwd_hash.inblock= false;
        rev_hash.inblock= false;

        fwd_hash.orientation = RIGHT;
        rev_hash.orientation = LEFT;    // call LEFT reverse

        hash_list.push_back(fwd_hash);
        hash_list.push_back(rev_hash);

      }
    }
  }


  hash_list.shrink_to_fit();      // neaten up vector


  // find global min and max hash value...
  unsigned int gmin, gmax;
  gmin = std::numeric_limits<unsigned int>::max();
  gmax = 0;

  for(hashtag hashed_substring : hash_list) {      // M.C.n
    if(hashed_substring.hash_id < gmin) {
      gmin = hashed_substring.hash_id;
    }
    if(hashed_substring.hash_id > gmax) {
      gmax = hashed_substring.hash_id;
    }
  }

  // using its hashvalue as an index, load every hashtag into an array
  // common 30bp seqs from different reads will have hashed to the same value

  cout << "make past here" << endl;
  cout << "gmax - gmin" << gmax - gmin << endl;
  vector<vector<hashtag>> commonReadsArray(gmax+1-gmin);
  cout << commonReadsArray.size();
  for (hashtag tag : hash_list) {
    cout << tag.hash_id << endl;
    commonReadsArray[tag.hash_id - gmin].push_back(tag);  // place tag in array
  }

  for(vector<hashtag> commonReads: commonReadsArray) {
    if (commonReads.size() < 2) { continue; } 


    set<set<unsigned int>> blocks;
    set<unsigned int> block;

    bool block_made = true;
    bool left = false;
    bool right = false;

    for(int i=0; i < commonReads.size(); i++) {
      if(commonReads[i].inblock) {continue;}

      block.insert(commonReads[i].read_id);
      commonReads[i].inblock = true;

      for(int j=0; j < commonReads.size(); j++) {
        string block_seed = 
        reads->getReadByIndex(commonReads[0].read_id,
        TUMOUR).substr(commonReads[0].offset, 30);

        string potential_member = 
        reads->getReadByIndex(commonReads[i].read_id,
        TUMOUR).substr(commonReads[i].offset, 30);


        // sequences should match, if not, its a collision, this is 
        // not a member, so continue
        if (!sequenceMatch(block_seed, potential_member)) { continue; }

        // group!!
        block_made = true;
        if(commonReads[i].orientation == RIGHT) {
          right = true;
        }
        else {
          left = true;
        }

        if(commonReads[j].orientation== RIGHT) {
          right = true;
        }
        else {
          left = true;
        }
        block.insert(commonReads[j].read_id);
        commonReads[j].inblock = true;
      }
  
      if(left && right && block_made) {
        blocks.insert(block);
      }
      
      // reset
      block.clear();
      block_made = left = right = false;
    }

    for (set<unsigned int> block : blocks) {
      if(block.size() >= 4) {
        HashedPreBlocks->push_back(block);
      }
    }
  }

}
bool BreakPointBlocks::sequenceMatch(string right, string left) {

  for(int i=0; i < right.size(); i++) {
    if(right[i] != left[i]) {return false;}
  }

  return true;
}



void BreakPointBlocks::unifyBreakpointBlocks() {
  unsigned int gmin, gmax;
  gmin = reads->n_tumour_reads;
  gmax = 0;

  // scan all blocks and identify global min and max
  for(set<unsigned int> block : *HashedPreBlocks) {
    set<unsigned int>::iterator it, end;
    it  = block.begin();
    end = block.end();

    for (; it != end; it++) {
      if (*it < gmin) {gmin = *it;}   // update global minimum
      if (*it > gmax) {gmax = *it;}   // update global max
    }
  }

  // unif_table inst.
  // unif_table keeps track of which groups have been unified to which others

  struct block_unif_info {    // the elements of the table
    unsigned int id;          // the id of the group
    unsigned int prev_unification_group;  // what group it was unified to
    bool unified;             // whether the group has been unified
  };

  vector<block_unif_info> unif_table;
  unif_table.reserve(HashedPreBlocks->size());
  //cout << "unif table size: " << CancerExtraction->size() << endl;

  for(unsigned int i=0; i < HashedPreBlocks->size(); i++) {
    block_unif_info b;
    b.id = i;
    b.unified = false;

    unif_table.push_back(b);          // intialize table
  }



  // make tagArray:
  // tagArray is a 2d array of the range gmin-gmax. 
  // The idea is that, each block id "tag" is added to this 2d array, 
  // such that, the read_id is used as an index. 
  // For example, given block_1 = {1, 2, 4, 8}. A tag to block_1
  // will be placed at positions 1, 2, 4, 8 in tag array. 
  // the tag is accompanied by a pointer to other blocks, and whether the block
  // has been merged. This allows effective merging

  // decl. tagArray as an array of vectors of pointers that will be used
  // to look up the unification status of a group in unif_table

  bool unification_occured = false;
  do {
    vector<vector<vector<block_unif_info>::iterator>> tagArray(gmax+1 -gmin);

    //cout << "unif_table status: " << endl;
    //for (block_unif_info block: unif_table) {
    //  cout << block.id << endl;
    //  cout << block.unified << endl;
    //  cout << block.prev_unification_group << endl;
    //}
    unification_occured = false;    // reset after loop

    // load read_id's
    for (unsigned int tag=0; tag < HashedPreBlocks->size(); tag++) {
      if (unif_table[tag].unified) {continue;}  // dont add info if merged

      // generate pointer to block
      vector<block_unif_info>::iterator ptr_to_block_info = 
        unif_table.begin() + tag;


      // set up iterators to next block, to read the read_id's
      set<unsigned int>::iterator read_id, end;
      read_id = (*HashedPreBlocks)[tag].begin();
      end = (*HashedPreBlocks)[tag].end();

      for (; read_id != end; read_id++) {   // for each read in block

        tagArray[*read_id- gmin].push_back(ptr_to_block_info);     
        // load the pointer to the group in
        // unif_table into the vector att he read_id'th index in tagArray
      }
    }
    //for(vector<vector<block_unif_info>::iterator> v : tagArray) {
    //  cout << "Tag elem size: " << v.size() << endl;
    //}

    // begin unification
    for (unsigned int i=0; i < (gmax+1-gmin); i++) {


      if (tagArray[i].size() < 2) {
          continue;
      } // nothing in list to unify

      else {

        // if the base group has already been unified with group x, then
        // unify all in list, with group x
        int unify_upon_this_group;

        // make syntax easier on the eye ;-)

        // check unif_table to see if the base group has been unified...
        if(tagArray[i][0]->unified) {  // if so, unify to that group
          unify_upon_this_group = tagArray[i][0]->prev_unification_group;  
        }
        else {  // if not, unify to the base
          unify_upon_this_group = tagArray[i][0]->id; 
        }

        // working backwards, unify the groups in block_tags, to 
        // unify_upon_this_group
        for (int k=tagArray[i].size()-1; k > 0; k--) {

          if (tagArray[i][k]->unified == true) {continue;} // already merged!

          unification_occured = true;   // going to unify groups


          set<unsigned int>::iterator it = 
            (*HashedPreBlocks)[tagArray[i][k]->id].begin(); // start of group

          set<unsigned int>::iterator end = 
            (*HashedPreBlocks)[tagArray[i][k]->id].end(); // end of group

          while(it != end) {
            (*HashedPreBlocks)[unify_upon_this_group].insert(*it); it++;
          }
          tagArray[i][k]->unified = true;
          tagArray[i][k]->prev_unification_group = unify_upon_this_group;
        }
      }
    }
  } while(unification_occured);


  // now, just extract merged groups and load to ComplementryUnified
  for (unsigned int i=0; i < unif_table.size(); i++) {
    if (unif_table[i].unified == false) {

      // pass the group into ComplementaryUnified
      FinalUnify->push_back((*HashedPreBlocks)[i]);
    }
  }

  // done with this data as processed into ComplementyUnfied
  //cout << "PRINTING NUMERICAL GROUPS" << endl;
  //printCancerExtractionGroups();
  delete CancerExtraction;
  CancerExtraction = nullptr;

  // finally, merging has completed. So, load the merged reads into 
  // ComplementaryGroups
  //for (int i=0; i < unif_table.size(); i++) {

  //  if(unif_table[i].unified == false) { // only add groups that were not unified

  //    // need to convert each read from each group into 
  //    // a tissue type, read_id pair, of a "block"

  //    set<unsigned int>::iterator it = (*CancerExtraction)[i].begin();
  //    set<unsigned int>::iterator end = (*CancerExtraction)[i].end();

  //    vector<pair<bool, unsigned int> > block;
  //    pair<bool, unsigned int> next_read;   
  //    while (it != end) {
  //      next_read.first = TUMOUR;
  //      next_read.second = *it;
  //      block.push_back(next_read);

  //      it++;
  //    }

  //    ComplementaryUnfied->push_back(block);
  //  }
  //}

}

void BreakPointBlocks::unifyComplementaryGroups() {
  unsigned int gmin, gmax;
  gmin = reads->n_tumour_reads;
  gmax = 0;

  // scan all blocks and identify global min and max
  for(set<unsigned int> block : *CancerExtraction) {
    set<unsigned int>::iterator it, end;
    it  = block.begin();
    end = block.end();

    for (; it != end; it++) {
      if (*it < gmin) {gmin = *it;}   // update global minimum
      if (*it > gmax) {gmax = *it;}   // update global max
    }
  }

  // unif_table inst.
  // unif_table keeps track of which groups have been unified to which others

  struct block_unif_info {    // the elements of the table
    unsigned int id;          // the id of the group
    unsigned int prev_unification_group;  // what group it was unified to
    bool unified;             // whether the group has been unified
  };

  vector<block_unif_info> unif_table;
  unif_table.reserve(CancerExtraction->size());
  //cout << "unif table size: " << CancerExtraction->size() << endl;

  for(unsigned int i=0; i < CancerExtraction->size(); i++) {
    block_unif_info b;
    b.id = i;
    b.unified = false;

    unif_table.push_back(b);          // intialize table
  }



  // make tagArray:
  // tagArray is a 2d array of the range gmin-gmax. 
  // The idea is that, each block id "tag" is added to this 2d array, 
  // such that, the read_id is used as an index. 
  // For example, given block_1 = {1, 2, 4, 8}. A tag to block_1
  // will be placed at positions 1, 2, 4, 8 in tag array. 
  // the tag is accompanied by a pointer to other blocks, and whether the block
  // has been merged. This allows effective merging

  // decl. tagArray as an array of vectors of pointers that will be used
  // to look up the unification status of a group in unif_table

  bool unification_occured = false;
  do {
    vector<vector<vector<block_unif_info>::iterator>> tagArray(gmax+1 -gmin);

    //cout << "unif_table status: " << endl;
    //for (block_unif_info block: unif_table) {
    //  cout << block.id << endl;
    //  cout << block.unified << endl;
    //  cout << block.prev_unification_group << endl;
    //}
    unification_occured = false;    // reset after loop

    // load read_id's
    for (unsigned int tag=0; tag < CancerExtraction->size(); tag++) {
      if (unif_table[tag].unified) {continue;}  // dont add info if merged

      // generate pointer to block
      vector<block_unif_info>::iterator ptr_to_block_info = 
        unif_table.begin() + tag;


      // set up iterators to next block, to read the read_id's
      set<unsigned int>::iterator read_id, end;
      read_id = (*CancerExtraction)[tag].begin();
      end = (*CancerExtraction)[tag].end();

      for (; read_id != end; read_id++) {   // for each read in block

        tagArray[*read_id- gmin].push_back(ptr_to_block_info);     
        // load the pointer to the group in
        // unif_table into the vector att he read_id'th index in tagArray
      }
    }
    //for(vector<vector<block_unif_info>::iterator> v : tagArray) {
    //  cout << "Tag elem size: " << v.size() << endl;
    //}

    // begin unification
    for (unsigned int i=0; i < (gmax+1-gmin); i++) {


      if (tagArray[i].size() < 2) {
          continue;
      } // nothing in list to unify

      else {

        // if the base group has already been unified with group x, then
        // unify all in list, with group x
        int unify_upon_this_group;

        // make syntax easier on the eye ;-)

        // check unif_table to see if the base group has been unified...
        if(tagArray[i][0]->unified) {  // if so, unify to that group
          unify_upon_this_group = tagArray[i][0]->prev_unification_group;  
        }
        else {  // if not, unify to the base
          unify_upon_this_group = tagArray[i][0]->id; 
        }

        // working backwards, unify the groups in block_tags, to 
        // unify_upon_this_group
        for (int k=tagArray[i].size()-1; k > 0; k--) {

          if (tagArray[i][k]->unified == true) {continue;} // already merged!

          unification_occured = true;   // going to unify groups


          set<unsigned int>::iterator it = 
            (*CancerExtraction)[tagArray[i][k]->id].begin(); // start of group

          set<unsigned int>::iterator end = 
            (*CancerExtraction)[tagArray[i][k]->id].end(); // end of group

          while(it != end) {
            (*CancerExtraction)[unify_upon_this_group].insert(*it); it++;
          }
          tagArray[i][k]->unified = true;
          tagArray[i][k]->prev_unification_group = unify_upon_this_group;
        }
      }
    }
  } while(unification_occured);


  // now, just extract merged groups and load to ComplementryUnified
  for (unsigned int i=0; i < unif_table.size(); i++) {
    if (unif_table[i].unified == false) {

      // pass the group into ComplementaryUnified
      ComplementaryUnfied->push_back((*CancerExtraction)[i]);
    }
  }

  // done with this data as processed into ComplementyUnfied
  //cout << "PRINTING NUMERICAL GROUPS" << endl;
  //printCancerExtractionGroups();
  delete CancerExtraction;
  CancerExtraction = nullptr;

  // finally, merging has completed. So, load the merged reads into 
  // ComplementaryGroups
  //for (int i=0; i < unif_table.size(); i++) {

  //  if(unif_table[i].unified == false) { // only add groups that were not unified

  //    // need to convert each read from each group into 
  //    // a tissue type, read_id pair, of a "block"

  //    set<unsigned int>::iterator it = (*CancerExtraction)[i].begin();
  //    set<unsigned int>::iterator end = (*CancerExtraction)[i].end();

  //    vector<pair<bool, unsigned int> > block;
  //    pair<bool, unsigned int> next_read;   
  //    while (it != end) {
  //      next_read.first = TUMOUR;
  //      next_read.second = *it;
  //      block.push_back(next_read);

  //      it++;
  //    }

  //    ComplementaryUnfied->push_back(block);
  //  }
  //}

}



void BreakPointBlocks::unifyReverseComplementaryReads() {

  for(unsigned int i = 0; i < ComplementaryUnfied->size(); i++ ) {

    // if the first read of the ith block of ComplementryUnified 
    // is oposite pair of the first read from the jth block of ComplementryUnified

    for(unsigned int j = i; j < ComplementaryUnfied->size(); j++) {
      if (i == j) {continue;} // dont self compare

      // if pair type of j and i is the same, we dont want to compare them
      std::set<unsigned int>::iterator i_it, j_it;
      i_it = (*ComplementaryUnfied)[i].begin();
      j_it = (*ComplementaryUnfied)[j].begin();

      if( (*reads->TumourMateOrder)[*i_it] ==   // skip if common mate 
          (*reads->TumourMateOrder)[*j_it]) {
        continue;
      }


      // groups are opposite pair type, so need rev comp comparison

      if(thirtyBasePairOverlap(*i_it, *j_it)) {

        // unify groups
        // unify to i
        for (;j_it != (*ComplementaryUnfied)[j].end(); j_it++) {
          (*ComplementaryUnfied)[i].insert(*j_it);  // add unique elems j to i
        }

        // unify to j
        for (;i_it != (*ComplementaryUnfied)[i].end(); i_it++) {
          (*ComplementaryUnfied)[j].insert(*i_it);  // add unique elems i to j
        }
      }
    }
  }


  // one compared reverse complementary sequences, make unique groups
  // by passing the groups through a set

  set<set<unsigned int>> filter_uniques;
  for(unsigned int i=0; i < ComplementaryUnfied->size(); i++) {
    filter_uniques.insert((*ComplementaryUnfied)[i]);
  }

  cout << "FILTER UNIQUES" << endl;
  for(set<unsigned int> s : filter_uniques) {
    for(unsigned int n : s) {
      cout << n << ", ";
    }
    cout << endl;
//  }




  // after passed through set groups unique, so place into vector

  // pass through each group
  for(set<unsigned int> group : filter_uniques) {

    set<pair<bool, unsigned int>> block;
    pair<bool, unsigned int> tumour_read;

    // load each read into next block
    for(unsigned int read_index : group) {
      tumour_read.first = TUMOUR;
      tumour_read.second = read_index;

      block.insert(tumour_read);
    }
    // load block into BreakPointBlocks 
    BreakpointBlocks->push_back(block);
  }
}


bool BreakPointBlocks::thirtyBasePairOverlap(unsigned int lhs, 
                           unsigned int rhs) {

  // string as rev comp of each other, by getting the rev comp of 
  // one, the two string should represent the same strand, and
  // will have a match of 30bp if they share the same genomic location

  string rhs_read_rev_comp = reverseComplementString(      // rev comp of right
                     reads->getReadByIndex(rhs, TUMOUR)
                     );
  string lhs_read = reads->getReadByIndex(lhs, TUMOUR);
  


  // window along 30bp along sequence looking for match
  for(unsigned int i=0; i < lhs_read.size() - 30; i++) {
    string lhs_sub = lhs_read.substr(i, 30);

    // if we find a 30bp match...
    if(rhs_read_rev_comp.find(lhs_sub) != std::string::npos) {
      return true;    // report match found
    }
  }


  return false;   // didnt find a match
}

int BreakPointBlocks::computeLCP(read_tag a, read_tag b) {
  string::const_iterator a_it, b_it, a_end, b_end;
  const string &astr = reads->getReadByIndex(a.read_id, TUMOUR);
  const string &bstr = reads->getReadByIndex(b.read_id, TUMOUR);

  int a_inc, b_inc;
  bool a_rc = false, b_rc = false;

  if (a.orientation == RIGHT) {
    a_it = astr.begin() + a.offset;
    a_end = astr.end();
    a_inc = 1;
  }
  else { 
    a_it = astr.begin() + a.offset + reads->getMinSuffixSize() - 1;
    a_end = astr.begin();
    a_inc = -1;
    a_rc = true;

  }

  if (b.orientation == RIGHT) {
    b_it = bstr.begin() + b.offset;
    b_end = bstr.end();
    b_inc = 1;
  }
  else {
    b_it = bstr.begin() + b.offset + reads->getMinSuffixSize() - 1;
    b_end = bstr.begin();
    b_inc = -1;
    b_rc = true;
  }

  for(; revCompCharacter(*a_it, a_rc) == revCompCharacter(*b_it, b_rc) &&
        b_it != b_end &&
        a_it != a_end; 
        a_it += a_inc, b_it += b_inc, lcp++);

}


int BreakPointBlocks::computeLCP(read_tag a, read_tag b) {
  string::const_iterator a_it, b_it, a_end, b_end;
  const string &astr = reads->getReadByIndex(a.read_id, TUMOUR);
  const string &bstr = reads->getReadByIndex(b.read_id, TUMOUR);
  int a_inc, b_inc;
  bool a_rc = false, b_rc = false;

  if (a.orientation == RIGHT) {
    a_it = astr.begin() + a.offset;
    a_end = astr.end();
    a_inc = 1;
  }
  else {
    a_it = astr.end() - 2 - a.offset; // 2 avoids dollar symbol
    a_end = astr.begin();
    a_inc = -1;
    a_rc = true;
  }

  if (b.orientation == RIGHT) {
    b_it = bstr.begin() + b.offset;
    b_end = bstr.end();
    b_inc = 1;
  }
  else {
    b_it = bstr.end() - 2 - b.offset; // 2 avoids dollar symbol
    b_end = bstr.begin();
    b_inc = -1;
    b_rc = true;
  }

  int lcp = 0;
  for(; revCompCharacter(*a_it, a_rc) == revCompCharacter(*b_it, b_rc) &&
        b_it != b_end &&
        a_it != a_end; 
        a_it += a_inc, b_it += b_inc, lcp++);

  return lcp;
}
void BreakPointBlocks::makeBreakPointBlocks(){
  if(CancerExtraction.size() == 0) {
    cout << "No mutations were identified " << endl;
    exit(1);
  }
 

  std::hash<std::string> hash;
  multimap<string, read_tag> mutation_grouper;

  for(unsigned int read_index : CancerExtraction) {
    string cancer_read = reads->getReadByIndex(read_index, TUMOUR);
    for(int i=0; i < cancer_read.size()-30; i++) {    // make 30bp fwd/rev

      pair<string, read_tag> fwd_key_val, rev_key_val;

      // forward and reverse pair of each 30pb window
      string sub_str_fwd = cancer_read.substr(i, 30);
      string sub_str_rev = reverseComplementString(sub_str_fwd);
      unsigned long long fwd_hash, rev_hash;
      fwd_hash = hash(sub_str_fwd);
      rev_hash = hash(sub_str_rev);
    

      // init forward tag
      read_tag fwd_tag;
      fwd_tag.read_id = read_index;
      fwd_tag.offset = i;
      fwd_tag.orientation = RIGHT;
      fwd_tag.tissue_type = TUMOUR;

      // init reverse tag
      read_tag rev_tag;
      rev_tag.read_id = read_index;
      rev_tag.offset = i;
      rev_tag.orientation = LEFT;
      rev_tag.tissue_type = TUMOUR;

      // make pairs and load into map
      fwd_key_val.first = sub_str_fwd;
      fwd_key_val.second = fwd_tag;

      rev_key_val.first = sub_str_rev;
      rev_key_val.second = rev_tag;

      mutation_grouper.insert(fwd_key_val);
      mutation_grouper.insert(rev_key_val);

    }
 }

  // now loaded into map, extract blocks
  bpBlock block;

  multimap<string, read_tag>::iterator it = mutation_grouper.begin();  

  unsigned int block_id = 0;
  while(it != std::prev(mutation_grouper.end())){
    bool left_mate = false, right_mate = false;

    if(it->first == std::next(it)->first) { // then same genomic location

      string seed = it->first;
      block.insert(it->second);           

      if(it->second.orientation) {
        right_mate = true;
      }
      else {
        left_mate = true;  
      }

      // extend group, all equal to seed
      it++;  
      while(seed == it->first) {
        block.insert(it->second);

        if(it->second.orientation) {
          right_mate = true;
        }
        else {
          left_mate = true;
        }

        it++;
        if(it == mutation_grouper.end()) {break;}
      }    
    
      // only extract double orientation groups with CTR 4
      if(left_mate & right_mate &&  block.size() >= 4) {
        block.id = block_id;
        BreakPointBlocks.push_back(block);
        block_id++;
      }

      block.clear();
    }

    else {
      it++;
    }

    if(it == mutation_grouper.end()) {break;}
  }


}
*/
