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
#include <list>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <iterator>

#include <zlib.h>
#include <omp.h>

#include "divsufsort64.h"
#include "util_funcs.h"
#include "SNVIdentifier.h"
#include "GenomeMapper.h"
#include "ssw_cpp.h"

using namespace std;

const char TERM = '$';
const int BASE33_CONVERSION = 33;



SNVIdentifier::SNVIdentifier(GSA &_gsa,
                                     string outpath,
                                     string const & basename,
                                     char mpq, int g1, int g2, int cut,
                                     int t, int mlcp, double e, double a):
                                     MIN_PHRED_QUAL(mpq), 
                                     PRI_MSS(g1), 
                                     AUX_MSS(g2), 
                                     COVERAGE_UPPER_THRESHOLD(cut),
                                     N_THREADS(t),
                                     MAX_SNPS(mlcp),
                                     ECONT(e), 
                                     ALLELIC_FREQ_OF_ERROR(a) {
  if (PRI_MSS > AUX_MSS) AUX_MSS = PRI_MSS;
  if (outpath[outpath.size()-1] != '/') outpath += "/";
  gsa = &_gsa;    
  cout << "Extracting reads from tumour-suffix enriched sections: pMSS = "
       << PRI_MSS << endl;
  extractCancerSpecificReads(); 

  cout << "Constructing auxiliary suffix array..." << endl;
  seedBreakPointBlocks();
  CancerExtraction.clear();
  cout << "Number of variant blocks extracted: " << SeedBlocks.size() 
       << endl;

  cout << endl << endl << "Consensus pair construction and filtering." << endl;
  cout << "Minimum phred score of included bases = " 
       << ((int) MIN_PHRED_QUAL) - BASE33_CONVERSION << endl;
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
    buildConsensusPairs(from, to, threadWork);
#pragma omp ordered
    {
      consensus_pairs.insert(consensus_pairs.end(), threadWork.begin(), threadWork.end());
    }
  }
  cout << "Number of filtered consensus pairs: " << consensus_pairs.size() << endl;

  string fastq_data;
  buildFastqAndTumourCNSData(consensus_pairs, fastq_data);
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
void SNVIdentifier::buildConsensusPairs(shared_ptr<bpBlock> * block,
    shared_ptr<bpBlock> * end,
    vector<consensus_pair> & threadWork) {
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
    // check prior to masking to remove pairs for which containsIndel
    // would be uneccessarily applied
    if (noSNV(pair)) {
      continue;
    }

    if (containsIndel(pair.mutated, pair.non_mutated)) {
      continue;
    }
    maskLowQualityPositions(pair);

    if (noSNV(pair)) {
      continue;
    }
    threadWork.push_back(pair);
  }
}

bool SNVIdentifier::singleIndel(string const & cigar) {
  int n_indels = 0;
  int mismatch_or_clip = 0;
  int pos = 0;

  // mismatch
  while ((pos = cigar.find('X', pos)) != string::npos) {
    pos++;
    mismatch_or_clip++;
  }
  pos = 0;
  // clip
  while ((pos = cigar.find('S', pos)) != string::npos) {
    pos++;
    mismatch_or_clip++;
  }
  pos = 0;
  // Count insertions
  while ((pos = cigar.find('I', pos)) != string::npos) {
    pos++;
    n_indels++;
  }
  // Count deletions
  pos = 0;
  while ((pos = cigar.find('D', pos)) != string::npos) {
    pos++;
    n_indels++;
  }
  return ((n_indels == 1 && mismatch_or_clip == 0) ? true : false);
}

bool SNVIdentifier::containsIndel(string const & tumour, string const & control) {
  // if (consecutiveMismatch()) return true;
  // requirement of SW library 
  int32_t maskLen = strlen(tumour.c_str())/2;
  maskLen = maskLen < 15 ? 15 : maskLen;
  
  /* Since we are interested only in identifying Indels, set
   * gap/extension penalty to zero. Accordingly, if an indel exists, the
   * gaped alignemnt will be the highest scoring.*/
  StripedSmithWaterman::Aligner aligner(1,1,0,0);
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;
  aligner.Align(tumour.c_str(), control.c_str(), control.size(), filter, &alignment, maskLen);
  if (singleIndel(alignment.cigar_string)) return true;
  else return false;
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
  if (numLowQuality > MAX_SNPS) return true;
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
  if (dir == true) {
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
    if (nreads < AUX_MSS) {
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
    if (idx < 0 || idx > gsa->len(r.read_id)-2) continue;
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
      if (freq_matrix[pos*4 + base] < AUX_MSS) {
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

    if (extension - seed_index == 1 && PRI_MSS   == 1 && gsa->tissuetype(seed_index) == TUMOUR) {
      threadExtr.insert(gsa->read_id_of_suffix(gsa->sa_element(seed_index)));           // extract read
    }
    else if (c_reads >= PRI_MSS  && (h_reads / (c_reads + h_reads)) <= ECONT)  {
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
    if (PRI_MSS  == 1 && gsa->tissuetype(gsa->sa_element(seed_index)) == TUMOUR) {
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

  int64_t concat_idx = 0;
  it++;
  for (; it != CancerExtraction.end(); it++) {
    concat += gsa->get_suffix_string(*it);
    concat += reverseComplementString(gsa->get_suffix_string(*it));
    concat += '#';
    concat_idx += gsa->len(*std::prev(it)) * 2;
    bsa.push_back(pair<int64_t, int64_t>(*it, concat_idx));
  }

  int64_t dSA_sz = concat.size();
  int64_t * dSA = (int64_t*) std::malloc(dSA_sz*sizeof(int64_t));
  divsufsort64((uint8_t*)const_cast<char*>(concat.c_str()), dSA,(int64_t)dSA_sz);

  // void function changes both dSA and dSA_sz
  remove_short_suffixes(dSA, dSA_sz, gsa->get_min_suf_size(), concat);

  omp_set_num_threads(N_THREADS);
  int64_t elementsPerThread = dSA_sz / N_THREADS;
  set<shared_ptr<bpBlock>, vBlockComp> result;
  cout << "Extracting variant blocks: aMSS = " << AUX_MSS << endl;
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
  if (omp_get_thread_num() == 0) {
    cout << "Removing redundant variant blocks" << endl;
  }
    result.insert(twork.begin(), twork.end());
  }
  SeedBlocks.insert(SeedBlocks.end(), result.begin(), result.end());
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
    if (ext - seed_idx >= AUX_MSS) {
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
  if (seed_idx == (dSA + dSA_sz - 1) && AUX_MSS == 1) {
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

//DEBUG
void SNVIdentifier::printExtractedCancerReads() {
  ofstream of("/data/ic711/extracted_cancer_read_ids.txt");
  for (unsigned int rid : CancerExtraction) of << rid << '\n';
  of << endl;
  of.close();
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

