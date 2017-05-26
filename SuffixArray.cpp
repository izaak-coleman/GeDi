// SuffixArray.cpp
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

// for use of radixSA
#include <sys/types.h>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <algorithm>
#include <thread>

#include "radix.h"
#include "Suffix_t.h"
#include "util_funcs.h"
#include "SuffixArray.h"
#include "Reads.h"
#include "benchmark.h"

using namespace std;

SuffixArray::SuffixArray(ReadPhredContainer &reads, int m, int n):
N_THREADS(n), MIN_SUFFIX(m) {
START(SuffixArray_SuffixArray);
  cout << "MIN SUFFIX " <<  MIN_SUFFIX << endl;
  this->reads = &reads;
  cout << "Starting parallelGenRadixSA:" << endl;
  parallelGenRadixSA(MIN_SUFFIX);
  cout << "GSA1 size: " << SA.size() << endl;
COMP(SuffixArray_SuffixArray);
}

void SuffixArray::parallelGenRadixSA(int min_suffix) {
START(SuffixArray_parallelGenRadixSA);
  vector<thread> BSA_and_SA;
  unsigned long long *radixSA;   // Pointer to suffix array
  unsigned int startOfTumour;    // Index of tumour reads in concat
  unsigned int radixSASize;

  // Construct binary search arrays and suffix array
  cout << "Making SA" << endl;
  vector<pair<unsigned int, unsigned int> > healthyBSA;
  vector<pair<unsigned int, unsigned int> > tumourBSA;
  BSA_and_SA.push_back(
      std::thread(&SuffixArray::buildBinarySearchArrays, this, 
        &healthyBSA, &tumourBSA)
  );
  BSA_and_SA.push_back(
      std::thread(&SuffixArray::generateParallelRadix, this, 
        &radixSA, &startOfTumour, &radixSASize)
  );
  for(auto & t : BSA_and_SA) t.join();
  // Construct Generalized Suffix Array
  cout << "radix sa size " << radixSASize << endl;
  cout << "Making GSA" << endl;
  vector<thread> workers;
  vector<vector<Suffix_t>> arrayBlocks(N_THREADS, vector<Suffix_t>());
  unsigned int elementsPerThread = (radixSASize/N_THREADS);
  unsigned int from{0}, to {elementsPerThread};
  for(unsigned int i = 0; i < N_THREADS; i++) {
    workers.push_back(
    std::thread(&SuffixArray::transformSuffixArrayBlock, this, &arrayBlocks[i], 
        &healthyBSA, &tumourBSA, radixSA, from, to, startOfTumour, min_suffix)
    );
    // Set up next threads block bounds
    from = to;
    if(i == N_THREADS - 2) {
      to = radixSASize;
    }
    else {
      to += elementsPerThread;
    }
  }
  for(auto &thread : workers) thread.join();

  // Load thread work into GSA
  delete radixSA;
  for(vector<Suffix_t> & b : arrayBlocks) {
    SA.insert(SA.end(), b.begin(), b.end());
    b.clear();
  }
  SA.shrink_to_fit();
COMP(SuffixArray_parallelGenRadixSA);
}

void SuffixArray::transformSuffixArrayBlock(vector<Suffix_t> *block, 
    vector<pair<unsigned int, unsigned int> > *healthyBSA,
    vector<pair<unsigned int, unsigned int> > *tumourBSA,
    unsigned long long *radixSA, unsigned int from, 
    unsigned int to, unsigned int startOfTumour, int min_suf) {
START(SuffixArray_transformSuffixArrayBlock);
  for(unsigned int i=from; i < to; i++) {
    // extract mapping, determining which read the suffix belongs to
    Suffix_t s;
    pair<unsigned int, unsigned int> rcPair;
    if (radixSA[i] < startOfTumour) { // HEALTHY
      rcPair = binarySearch(*healthyBSA,radixSA[i]);
      s.offset = radixSA[i] - rcPair.second;
      if (reads->getReadByIndex(rcPair.first, HEALTHY).size() - s.offset <=
          reads->getMinSuffixSize()) continue;
      s.type = HEALTHY;
    }
    else  { // TUMOUR
      rcPair = binarySearch(*tumourBSA, radixSA[i] - startOfTumour);
      s.offset = (radixSA[i] - startOfTumour) - rcPair.second;
      if (reads->getReadByIndex(rcPair.first, TUMOUR).size() - s.offset
          <= reads->getMinSuffixSize()) continue;
      s.type = TUMOUR;
    }
    s.read_id = rcPair.first;
    block->push_back(s);
  }
COMP(SuffixArray_transformSuffixArrayBlock);
}

void SuffixArray::buildBinarySearchArrays(
                  vector<pair<unsigned int, unsigned int> > *healthyBSA, 
                  vector<pair<unsigned int, unsigned int> > *tumourBSA) {
  generateBSA(*healthyBSA, HEALTHY);
  generateBSA(*tumourBSA, TUMOUR);
}

void SuffixArray::generateBSA(
     vector<pair<unsigned int, unsigned int>> &BSA, bool type) {
START(SuffixArray_generateBSA);
  unsigned int indexInConcat {0};
  BSA.reserve(reads->getSize(type));
  BSA.push_back(pair<unsigned int, unsigned int>(0,0));
  for(unsigned int i = 1; i < reads->getSize(type); i++) { 
    indexInConcat += reads->getReadByIndex(i-1, type).size();
    BSA.push_back(pair<unsigned int, unsigned int>(i,indexInConcat));
  }
COMP(SuffixArray_generateBSA);
}

void SuffixArray::generateParallelRadix(unsigned long long **radixSA, 
                                        unsigned int *startOfTumour,
                                        unsigned int *sizeOfRadixSA) {
START(SuffixArray_generateParallelRadix);
  // concatenate all reads
  string concat = concatenateReads(HEALTHY);
  *startOfTumour = concat.size();            // mark the end of healthy seqs 
  concat += concatenateReads(TUMOUR);
  *sizeOfRadixSA = concat.size();
  // Build SA
  *radixSA = Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();
COMP(SuffixArray_generateParallelRadix);
}

pair<unsigned int, unsigned int> 
SuffixArray::binarySearch(vector<pair<unsigned int, unsigned int> > &BSA, 
                  unsigned int suffixIndex) {
START(SuffixArray_binarySearch);
  unsigned int right = BSA.size();
  unsigned int left = 0;
  unsigned int mid;
  // binary search to home in on read
  while(left < right) {
    mid = left + ((right - left) / 2);
    if(suffixIndex == BSA[mid].second)  return BSA[mid];
    else if (suffixIndex > BSA[mid].second) left = mid+1;
    else right = mid;
  }
  left--;
  return BSA[left]; // left on the tuple
COMP(SuffixArray_binarySearch);
}

string SuffixArray::concatenateReads(bool type) {
START(SuffixArray_concatenateReads);
  string concat = "";
  for(int i=0; i < reads->getSize(type); i++) {
    concat += reads->getReadByIndex(i, type);
  }
  return concat;
COMP(SuffixArray_concatenateReads);
}

Suffix_t & SuffixArray::getElem(int index) {
  if (index >= SA.size() || index < 0) {
    cout << "getElem() out of bounds" << endl;
    exit(1);
  }
  return SA[index];
}

unsigned int SuffixArray::getSize() {
  return SA.size();
}

// End of file

/*
// debug ints for printing GSA
static const int READ_ID_IDX = 0;
static const int OFFSET_IDX  = 1;
static const int TYPE_IDX = 2;
static const int NUM_FIELDS = 3;


void SuffixArray::generalizedRadixSA(vector<Suffix_t> *TissueSA, bool type,
    uint8_t min_suf) {

  // concatenate all reads to one giant read
  string concat = concatenateReads(type);


  // generate the Binary Search Array used to transform output from 
  // suffix array (one string) to generalized suffix array
  // (multi string suffix array)
  vector<pair<unsigned int, unsigned int>> readToConcatMap;
  generateBSA(readToConcatMap, type);


  // Builds SA
  cout << "Generating " << ((type) ? "healthy" : "tumour" ) << " suffix" << 
    " array. " << endl;

  unsigned long long *SA;
  SA = Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();



  cout << "Generating " << ((type) ? "healthy" : "tumour" ) << 
    " general suffix array. " << endl;

  TissueSA->reserve(reads->getSize(type) * 50);


  // transform suffix array to generalized suffix array
  for(unsigned int i=0; i < concat.size(); i++) {

    // extract mapping, determining which read the suffix belongs to
    pair<unsigned int, unsigned int> read_concat_tup = 
                       binarySearch(readToConcatMap, SA[i]);



    // suffixes with length >= 30 only
    if((SA[i] - read_concat_tup.second) <=
        (reads->getReadByIndex(read_concat_tup.first, type).size() - min_suf)) {
      Suffix_t s;
      s.read_id = read_concat_tup.first;
      s.offset = SA[i] - read_concat_tup.second; 
      s.type = type;

      TissueSA->push_back(s);
    }
    else { // suffix was less than 30pb long so we dont want it  
      continue;
    }
  }

  TissueSA->shrink_to_fit();

  delete [] SA;   // done with local suffix array
}
void SuffixArray::constructTotalRadixSA(uint8_t min_suffix) {

  // local suffix arrays
  vector<Suffix_t> healthy_SA;
  vector<Suffix_t> tumour_SA;

  // buill local suffix arrays in parallel

  vector<thread> tissue_SA_threads;

  // thread for healthy 
  tissue_SA_threads.push_back(
      std::thread(&SuffixArray::generalizedRadixSA, this,
                  &healthy_SA, 
                  HEALTHY, min_suffix)
      );
 
  // thread for tumour
  tissue_SA_threads.push_back(
      std::thread(&SuffixArray::generalizedRadixSA, this,
        &tumour_SA, 
        TUMOUR, min_suffix)
      );

  // wait for threads to finish
  for(auto &thread : tissue_SA_threads) {
    thread.join();
  }



  SA.reserve(healthy_SA.size() + tumour_SA.size()); // make room

  bool end_of_healthy = false;
  bool end_of_tumour = false; // not for long...lets hope ;)

  unsigned int tind=0, hind=0;

  cout << "Merging tissue general suffix arrays to build final suffix " << 
    "array " << endl;


  while ((tind < tumour_SA.size()) || (hind < healthy_SA.size())) {

    // bound checks...
    if(hind == healthy_SA.size()) {
      end_of_healthy = true;
    }

    if(tind == tumour_SA.size()) {
      end_of_tumour = true;
    }


    // one has reached end, so add all of other
    if(end_of_healthy) {
      SA.push_back(tumour_SA[tind]);
      tind++;
    }

    else if(end_of_tumour) {
      SA.push_back(healthy_SA[hind]);
      hind++;
    }

    else {  // noone reached end so add based on lex order
      string::iterator h_start, h_end, t_start, t_end;

      t_start = reads->returnStartIterator(tumour_SA[tind]);
      t_end   = reads->returnEndIterator(tumour_SA[tind]);

      h_start = reads->returnStartIterator(healthy_SA[hind]);
      h_end   = reads->returnEndIterator(healthy_SA[hind]);


      if(lexicographical_compare(t_start, t_end, h_start, h_end)) {
        SA.push_back(tumour_SA[tind]);
        tind++;
      }
      else {
        SA.push_back(healthy_SA[hind]);
        hind++;
      }

    }

  }

  cout << "Done with suffix array build" << endl;

}

// ---------------- NAIVE MERGESORT FUNCS ----------------------------


void SuffixArray::sort(unsigned int from, unsigned int to, unsigned int level) {
  if ((to - from) <= 1) {return;} // dividing reaches single elem, hit base case

  // otherwise...keep dividing
  int mid = (from + to) / 2;
  sort(from, mid, level+1);
  sort(mid, to, level+1);

  string level_name;
  stringstream s;

  s << level;
  level_name = "level_";
  level_name += s.str();

  //MSTART(SuffixArray_level_name);
  merge(from, mid, to);
  //MEND(level_name);
  //MTIME(level_name);
  //MPRINT(level_name);
}

void SuffixArray::merge(unsigned int from, unsigned int mid, unsigned int to) {
  // make out of place copies of sa section 
  vector<Suffix_t> *left = copyOf(from, mid);
  vector<Suffix_t> *right = copyOf(mid, to);


  unsigned int left_ptr=0, right_ptr=0, sa_ptr=from;
  bool end_of_left = false;      // stop range bound errors 
  bool end_of_right = false;

  // add suffix_t's back to sa in lexicographical order
  while (sa_ptr < to) {

    // right finished... keep adding lefts elems
    if (!end_of_left && end_of_right) {
      SA[sa_ptr++] = (*left)[left_ptr++];
    }

    // left finished... keep adding rights elems
    else if (!end_of_right && end_of_left) {
      SA[sa_ptr++] = (*right)[right_ptr++];
    }

    // left lexiocographically before right element, so add left next
    else if (!end_of_left &&
             lexCompare((*left)[left_ptr], (*right)[right_ptr])
             ) {

      SA[sa_ptr++] = (*left)[left_ptr++];

    }
    else if(!end_of_right){   // right lexicographcially before left
      SA[sa_ptr++] = (*right)[right_ptr++];
    }
    else{
      cout << "merge sort error" << endl;

    }

    // check bounds
    if (left_ptr == left->size()) {
      end_of_left = true;
    }
    if (right_ptr == right->size()) {
      end_of_right = true;
    }
  }

  delete left;
  delete right;
}



vector<Suffix_t>* SuffixArray::copyOf(unsigned int from, unsigned int to){
  // make section copy of SA on heap 
  vector<Suffix_t> *SA_section_ptr;
  SA_section_ptr = new vector<Suffix_t>;
  SA_section_ptr->reserve(to-from);

  // load with section
  for (unsigned int i = from; i < to; i++) {
    SA_section_ptr->push_back(SA[i]);
  }
  SA_section_ptr->shrink_to_fit();    // Keep size down 
  return SA_section_ptr;
}

bool SuffixArray::lexCompare(Suffix_t &lhs, Suffix_t &rhs) {
  // Generate pointers to lhs and rhs suffixes in reads
  string::iterator lhs_iter  = reads->returnStartIterator(lhs);
  string::iterator lhs_end   = reads->returnEndIterator(lhs);
  string::iterator rhs_iter  = reads->returnStartIterator(rhs);
  string::iterator rhs_end   = reads->returnEndIterator(rhs);

  for( ; (lhs_iter != lhs_end && rhs_iter != rhs_end); lhs_iter++, rhs_iter++){
    // lex compare character
    if (*lhs_iter < *rhs_iter) { return true; }
    if (*rhs_iter < *lhs_iter) { return false; }
    // equiv char so move to next...
    }

    // One is prefix of other, return the prefix as higher suffix
    return (lhs_iter == lhs_end) && (rhs_iter != rhs_end);
}

void SuffixArray::lexMergeSort() {
   sort(0, SA.size(), 0);      // start recursive mergesort
}

void SuffixArray::loadUnsortedSuffixes(uint8_t min_suffix) {

    // Read length is ~100 bp, and stoping at 100 - min_suffix
    SA.reserve(    // reserve size for tumour + healthy arrays
        (reads->getSize(HEALTHY) + reads->getSize(TUMOUR)) * (100 - min_suffix)
        ); 

  // Loop through each read in HEALTHY set, and add suffixes
  for(unsigned int read_id = 0; read_id < reads->getSize(HEALTHY); read_id++){

    for(uint16_t offset = 0; 
        offset <= (reads->getReadByIndex(read_id, HEALTHY).size() - min_suffix)
        && offset < reads->getReadByIndex(read_id, HEALTHY).size(); // Safety bounds
        offset++) {

      // Construct Suffix_t with correct info 
      Suffix_t suf;
      suf.read_id    = read_id;
      suf.offset     = offset;
      suf.type = HEALTHY;
      // add to SA
      SA.push_back(suf);
    }
  }
  for(unsigned int read_id = 0; read_id < reads->getSize(TUMOUR); read_id++){

    for(uint16_t offset = 0; 
        offset <= (reads->getReadByIndex(read_id, TUMOUR).size() - min_suffix)
        && offset < reads->getReadByIndex(read_id, TUMOUR).size(); // Safety bounds
        offset++) {

      // Construct Suffix_t with correct info 
      Suffix_t suf;
      suf.read_id    = read_id;
      suf.offset     = offset;
      suf.type = TUMOUR;
      // add to SA
      SA.push_back(suf);
    }
  }

  // SA will now nolonger change size of the program. So make as
  // compact as possible 
  SA.shrink_to_fit();
}

void SuffixArray::buildGSAFile(vector<Suffix_t> &GSA, string filename) {
  if (filename.substr(filename.size() - EXT.size()) != EXT){
    filename += EXT;
  }
  ofstream gsa_file;
  gsa_file.open(filename);

  for(int i = 0; i < GSA.size(); i++) {
    string next_suf = "";
    next_suf += (to_string(GSA[i].read_id) + ",");
    next_suf += (to_string(GSA[i].offset)  + ",");
    next_suf += to_string(GSA[i].type);
    gsa_file << next_suf << endl;
  }
}

void SuffixArray::constructGSAFromFile(vector<Suffix_t> &GSA, string filename) {
  ifstream gsa_file;
  gsa_file.open(filename);
  string next_line;

  while (getline(gsa_file, next_line)) {
    vector<string> elements;
    split_string(next_line, ",", elements); // split.gsa line into suffix fields
    Suffix_t suf;

    if (elements.size() == NUM_FIELDS) {
      suf.read_id = stoi(elements[READ_ID_IDX]);    // load data
      suf.offset = stoi(elements[OFFSET_IDX]);
      suf.type = stoi(elements[TYPE_IDX]);
      GSA.push_back(suf);
    }
  }
}


void SuffixArray::printSuffixData() {
  for(Suffix_t s : SA) {
    cout << "read_id: " << s.read_id <<  " --- " 
         << "offset: "  << s.offset  <<  " --- "
         << std::boolalpha 
         << "tissue: " << ((s.type) ? "HEALTHY" : "TUMOUR") 
         << endl;
  }
}

void SuffixArray::printSuffixArray(std::string const& filename) { ofstream file(filename);
  for (Suffix_t & s : SA) {
    file << reads->returnSuffix(s) << endl;
  }
  file.close();
}
void SuffixArray::printSuffixes() {
  for(Suffix_t s : SA) {
    cout << reads->returnSuffix(s) << endl;
  }
}

void SuffixArray::printReadsInGSA(std::string const& filename) {
  ofstream ofHandle(filename.c_str());

  for (Suffix_t const& s : SA) {
    ofHandle << "(" << s.read_id << ","
             << ((s.type == HEALTHY) ? "H" : "T") 
             << ")" << std::endl;
  }
  ofHandle.close();
}

 */
