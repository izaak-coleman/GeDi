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
#include "string.h" // split_string()



#include "Suffix_t.h"
#include "util_funcs.h"
#include "SuffixArray.h"
#include "Reads.h"

#include "benchmark.h"

using namespace std;

// debug ints for printing GSA
static const int READ_ID_IDX = 0;
static const int OFFSET_IDX  = 1;
static const int TYPE_IDX = 2;
static const int NUM_FIELDS = 3;


static const string EXT = ".gsa";


SuffixArray::SuffixArray(ReadPhredContainer &reads, int min_suffix, 
                         int n_threads):
N_THREADS(n_threads),
MIN_SUFFIX(min_suffix) {
  cout << "MIN SUFFIX " << (int) min_suffix << endl;
  this->reads = &reads;      // store reads location
  cout << "Starting parallelGenRadixSA:" << endl;

  parallelGenRadixSA(min_suffix);
  
//  printReadsInGSA("/data/ic711/point2.txt");
}


SuffixArray::~SuffixArray() {
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

// PARALLELRADIXSACONSTRUCTION FUNCS


void SuffixArray::parallelGenRadixSA(int min_suffix) {

  vector<thread> BSA_and_SA;

  unsigned long long *radixSA;   // suffix array pointer
  unsigned int startOfTumour; // int marking start of tumour reads in concat
  unsigned int radixSASize;

  // Build binary search arrays and suffix array in parallel

  vector<pair<unsigned int, unsigned int> > healthyBSA;
  vector<pair<unsigned int, unsigned int> > tumourBSA;

  BSA_and_SA.push_back (
      std::thread(&SuffixArray::buildBinarySearchArrays, this, &healthyBSA, &tumourBSA)
      );

  BSA_and_SA.push_back (
      std::thread(&SuffixArray::generateParallelRadix, this, &radixSA, 
        &startOfTumour, &radixSASize)
      );

  for(auto &thread : BSA_and_SA) { 
    thread.join();
  }




  // begin parallel suffix array construction
  cout << "radix sa size " << radixSASize << endl;
  vector<thread> workers;
  vector<vector<Suffix_t>> array_blocks;
  unsigned int elements_per_thread = (radixSASize/N_THREADS);

  // initialze blocks
  for(int i=0; i < N_THREADS; i++) {
    vector<Suffix_t> init;
    init.reserve(elements_per_thread);      // make room
    array_blocks.push_back(init);
  }
  
  unsigned int from=0, to = elements_per_thread;
  for(unsigned int i=0; i < N_THREADS; i++) {
    // run worker thread
    workers.push_back(
    std::thread(&SuffixArray::transformSuffixArrayBlock, this, &array_blocks[i], 
        &healthyBSA, &tumourBSA, radixSA, from, to, startOfTumour, min_suffix)
    );
    // set up next worker thread
    from = to;
    if(i == N_THREADS - 2) {  // set up end chunk for final thread
      to = radixSASize;
    }
    else {
      to += elements_per_thread;
    }
  }

  // wait for workers
  for(auto &thread : workers) {
    thread.join();
  }

  delete radixSA;  // done with suffix array
  // Finally, load blocks into final SA in order
  for(int i=0; i < array_blocks.size(); i++) {
    for(int j=0; j < array_blocks[i].size(); j++) {
      SA.push_back(array_blocks[i][j]);
    }
    array_blocks[i].clear();
  }

  SA.shrink_to_fit();
}

void SuffixArray::transformSuffixArrayBlock(vector<Suffix_t> *block, 
    vector<pair<unsigned int, unsigned int> > *healthyBSA,
    vector<pair<unsigned int, unsigned int> > *tumourBSA,
    unsigned long long *radixSA, unsigned int from, 
    unsigned int to, unsigned int startOfTumour, int min_suf) {

  for(unsigned int i=from; i < to; i++) {
    // extract mapping, determining which read the suffix belongs to
    if(radixSA[i] < startOfTumour) { // HEALTHY
      pair<unsigned int, unsigned int > read_concat_tup = 
        binarySearch(*healthyBSA, radixSA[i]);
      Suffix_t s;
      s.offset = radixSA[i] - read_concat_tup.second;
      if (reads->getReadByIndex(read_concat_tup.first, HEALTHY).size() - s.offset
          <= reads->getMinSuffixSize()) { 
        continue;
      }
      else{
        s.read_id = read_concat_tup.first;
        s.type = HEALTHY;
        block->push_back(s);
      }

//      if((radixSA[i] - read_concat_tup.second) <=
//          (reads->getReadByIndex(read_concat_tup.first, HEALTHY).size() - min_suf)) {
//        s.read_id = read_concat_tup.first;
//        s.offset = radixSA[i] - read_concat_tup.second; 
//        s.type = HEALTHY;
//
//        block->push_back(s);
//      }
//      else { // suffix was less than 30pb long so we dont want it  
//        continue;
//      }
    }

    else {                      // TUMOUR Can logic be simplified??
      pair<unsigned int, unsigned int > read_concat_tup = 
        binarySearch(*tumourBSA, radixSA[i]-startOfTumour);

      Suffix_t s;
      s.offset = (radixSA[i] - startOfTumour) - read_concat_tup.second;
      if (reads->getReadByIndex(read_concat_tup.first, TUMOUR).size() - s.offset
          <= reads->getMinSuffixSize()) { 
        continue;
      }
      else{
        s.read_id = read_concat_tup.first;
        s.type = TUMOUR;
        block->push_back(s);
      }
    //  if(((radixSA[i] - startOfTumour) - read_concat_tup.second) <=
    //      (reads->getReadByIndex(read_concat_tup.first, TUMOUR).size() - min_suf)) {
    //    Suffix_t s;
    //    s.read_id = read_concat_tup.first;
    //    s.offset = (radixSA[i] - startOfTumour) - read_concat_tup.second; 
    //    s.type = TUMOUR;

    //    block->push_back(s);
    //  }
    //  else { // suffix was less than 30pb long so we dont want it  
    //    continue;
    //  }
    }
  }
}

void SuffixArray::buildBinarySearchArrays(
    vector<pair<unsigned int, unsigned int> > *healthyBSA, 
    vector<pair<unsigned int, unsigned int> > *tumourBSA) {
  generateBSA(*healthyBSA, HEALTHY);
  generateBSA(*tumourBSA, TUMOUR);

}

void SuffixArray::generateParallelRadix(unsigned long long **radixSA, 
                                        unsigned int *startOfTumour,
                                        unsigned int *sizeOfRadixSA) {
  // concatenate all reads
  string concat = concatenateReads(HEALTHY);
  *startOfTumour = concat.size();            // mark the end of healthy seqs 
  concat += concatenateReads(TUMOUR);
  *sizeOfRadixSA = concat.size();
  // Build SA
  *radixSA = Radix<unsigned long long>((uchar*) concat.c_str(), concat.size()).build();
}


// RADIXSAMERGECONSTRUCTOINFUNCS

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

pair<unsigned int, unsigned int> 
            SuffixArray::binarySearch(
                  vector<pair<unsigned int, unsigned int> > &BSA, 
                  unsigned int suffix_index) {

  unsigned int right = BSA.size();
  unsigned int left = 0;
  unsigned int mid;

  // binary search to home in on read
  while(left < right) {
    mid = left + ((right - left) / 2);

    if(suffix_index == BSA[mid].second) {
        return BSA[mid];
    }
    else if (suffix_index > BSA[mid].second) {
      left = mid+1;
    }
    else {
      right = mid;
    }
  }
  left--;
  return BSA[left]; // left should be on the seq
}

void SuffixArray::generateBSA(
     vector<pair<unsigned int, unsigned int>> &BSA, bool type) {

  
  unsigned int index_in_concat = 0;

  BSA.reserve(reads->getSize(type));
  pair<unsigned int, unsigned int> zeroPos(0,0);
  BSA.push_back(zeroPos);

  for(unsigned int i=1; i < reads->getSize(type); i++) { 
    pair<unsigned int, unsigned int> read_total;
    index_in_concat += reads->getReadByIndex(i-1, type).size();

    read_total.first = i;
    read_total.second = index_in_concat;    // start pos of read in concat

    BSA.push_back(read_total);
  }

}


string SuffixArray::concatenateReads(bool type) {
  string concat = "";

  for(int i=0; i < reads->getSize(type); i++) {
    concat += reads->getReadByIndex(i, type);
  }

  return concat;
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

  //MSTART(level_name);
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



// PUBLIC HELPER FUNCTIONS





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

void SuffixArray::printSuffixArray(std::string const& filename) {
  ofstream file(filename);
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
