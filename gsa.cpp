#include <iostream>
#include <fstream>
#include <string>
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm> 
#include <unistd.h>
#include <fcntl.h>
#include <zlib.h>
#include <omp.h>

#include "kseq.h"
#include "gsa.h"
#include "divsufsort64.h"
#include "util_funcs.h"

#define STR(str) #str
#define STRING(str) STR(str)

#define STR(str) #str
#define STRING(str) STR(str)


KSEQ_INIT(gzFile, gzread);    // initialize .gz parser

using namespace std;
static char const TERM_CHAR = '$';
static char const PHRED_DELIM = '\0';
static int64_t const MIN_SUF_LEN = 30;
static const          string HEALTHY_DATA = "C";
static const          string TUMOUR_DATA  = "T";
static const char     PHRED_20 = '5';
static const double   QUAL_THRESH = 0.1;
static const char     REMOVED_TOKEN = 'N';
static const string   ON = "on";

int64_t GSA::size() {
  return sa_sz;
}

int64_t GSA::get_min_suf_size() {
  return MIN_SUF_LEN;
}


GSA::GSA(string const& header_fname, int t, string const & ref, string const & em_filtering) {
  N_THREADS = t;
  max_read_len = 0;
  tsi = 0;
  sa_sz = 0;
  sa = nullptr;
  vector<string> h_fnames, t_fnames;
  read_header(header_fname, h_fnames, t_fnames);
  reportFilesDataset(h_fnames, HEALTHY);
  reportFilesDataset(t_fnames, TUMOUR);
  cout << endl << endl;
  cout << "Preprocessing." << endl;
  for (string const & f : h_fnames) {
    load_fq_data(f);
  }
  tsi = concat.size();
  if (em_filtering != ON) {
    for (string const & f : t_fnames) {
      load_fq_data(f);
    }
  }
  else {
    cout << "Running emfilter on: " << endl;
    for (string const & f : t_fnames) {
      cout << " - " + f + " ..." << endl;
      em_filter(ref, f);
    }
  }
  sa = (int64_t*) std::malloc(concat.size()*sizeof(int64_t));
  if (sa == nullptr) {
    cout << "Memory for suffix array allocation not available. Program terminating" << endl;
    exit(1);
  }
  cout << endl << endl << "SNV detection." << endl;
  cout << "Constructing primary suffix array..." << endl;
  sa_sz = concat.size();
  divsufsort64((uint8_t*)const_cast<char*>(concat.c_str()), sa, sa_sz);
  remove_short_suffixes(MIN_SUF_LEN);
}

void GSA::reportFilesDataset(vector<string> const & file_list, bool const tissue) {
  // Reports whether a fastq was loaded as tumour or control data.
  string type = (tissue == HEALTHY) ? "CONTROL" : "TUMOUR";
  for (string const & f : file_list) {
    cout << type + ": " + f << endl;
  }
}

void GSA::em_filter(string const & ref, string const & f) {
  // run bowtie2
  cout << "Bowtie2 output: " << endl;
  string pwd(STRING(PWD));
  string command_aln(pwd + "/bowtie2-2.3.4/bowtie2 -p 16 -x " + ref + " -U "
      + f + " -S " + f + ".sam");
  system(command_aln.c_str());
  string sf = f + ".sam";
  ifstream sam(sf.c_str()); // attach stream to samfile

  for (string line; getline(sam,line);) {
    if (line[0] == '@') continue;

    int64_t nm_index = line.find("NM:i:");
    if (nm_index != string::npos) {
      nm_index +=5;
      if (line[nm_index] == '0') continue; // read is exact match
    }
    // remaining reads are inexactly matched, or failed alignment
    vector<string> fields;
    split_string(line, "\t", fields);
    if (!fields[10].empty() && good_quality(fields[10])) {
      if (fields[1] == "16") {
        fields[9] = reverseComplementString(fields[9]);
        std::reverse(fields[10].begin(), fields[10].end());
      }
      add_fq_data(fields[9], fields[10]);
    }
  }
  sam.close();
  remove(sf.c_str()); // delete samfile
}


void GSA::load_fq_data(string const & fname) {
  gzFile data;
  data = gzopen(fname.c_str(), "r");
  kseq_t *fq = kseq_init(data);
  while (kseq_read(fq) >= 0) {
    if (fq->qual.l && good_quality(fq->qual.s)) {
      add_fq_data(fq->seq.s, fq->qual.s);
    }
  }
  kseq_destroy(fq);
  gzclose(data);
}

bool GSA::good_quality(string const & q) {
  double n_lowq_bases{0.0};
  for (int64_t i = 0; i < q.size(); i++) {
    if (q[i] < PHRED_20) ++n_lowq_bases; 
  }
  if ((n_lowq_bases / q.size()) > QUAL_THRESH) return false;
  else return true;
}

void GSA::add_fq_data(string const & seq, string const & qual) {
  vector<string> read_substrs;
  vector<string> phred_substrs;
  int left_arrow{0}, right_arrow{0};
  while((right_arrow = seq.find(REMOVED_TOKEN, left_arrow)) != string::npos) {
    if (left_arrow == right_arrow) {
      left_arrow++;
      continue;
    }
    else if ((right_arrow - left_arrow) < MIN_SUF_LEN) {
      left_arrow = right_arrow + 1;
      continue;
    }
    else {
      if ((right_arrow - left_arrow) > max_read_len) {
        max_read_len = right_arrow - left_arrow;
      }
      concat += seq.substr(left_arrow, right_arrow - left_arrow) + TERM_CHAR;
      phred  += qual.substr(left_arrow, right_arrow - left_arrow) + '\0';
      left_arrow = right_arrow + 1;
    }
  }
  if (seq.size() - left_arrow >= MIN_SUF_LEN)  {
    concat += seq.substr(left_arrow) + TERM_CHAR;
    phred  += qual.substr(left_arrow) + '\0';
  }
  if (seq.size() - left_arrow > max_read_len) {
    max_read_len = seq.size() - left_arrow;
  }
}

int64_t GSA::get_max_read_len() {
  return max_read_len;
}

void GSA::read_header(string const& header_fname,
                      vector<string> & h_fnames, vector<string> & t_fnames) {
  ifstream fin;
  fin.open(header_fname.c_str());
  string line;
  while (getline(fin, line)) {
    vector<string> fields;
    split_string(line, ",\t ", fields);
    if (fields[1] == HEALTHY_DATA) {
      h_fnames.push_back(fields[0]);
    } else if (fields[1] == TUMOUR_DATA) {
      t_fnames.push_back(fields[0]);
    } else {
      cout << "Invalid input found in " << header_fname 
           << ". Valid datatypes are either C or T." 
           << endl << "Program terminating." << endl;
      exit(1);
    }
  }
  fin.close();
}


void GSA::xorSwap(int64_t *x, int64_t *y) {
  if (x != y) {
    *x ^= *y;
    *y ^= *x;
    *x ^= *y;
  }
}

int64_t GSA::bubbleRemove(int64_t * const a, int64_t const sz, int64_t const invalid) {
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

void GSA::remove_short_suffixes(int64_t min_suffix_length) {
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
      if (len(*from) <= min_suffix_length) {
        *from = -1;
      }
      from++;
    }
  }
  sa_sz = bubbleRemove(sa, sa_sz, -1);
  int64_t* temp = (int64_t*) std::realloc(sa, sa_sz*sizeof(int64_t));
  if (temp == nullptr) {
    cout << "Not enough memory to construct suffix array. Program terminating."
      << endl;
    exit(1);
  } else {
    sa = temp;
  }
}

string GSA::get_suffix_string(int64_t i) {
  string::const_iterator it = concat.cbegin() + i;
  string s;
  while (*it != TERM_CHAR) {s += *it; it++;}
  return s + TERM_CHAR;
}

int64_t GSA::sa_element(int64_t pos) {
  if (pos < 0 || pos > sa_sz) return -1;
  return *(sa + pos);
}

bool GSA::tissuetype(int64_t const i) {
  if (i < 0 || i >= concat.size()) exit(1);
  return i < tsi;
}

string::const_iterator GSA::suffix_at(int64_t i) {
  if (i < 0 || i >= concat.size()) return concat.cend();
  return concat.cbegin() + i;
}

string::const_iterator GSA::phred_at(int64_t i) {
  if (i < 0 || i >= phred.size()) return phred.cend();
  return phred.cbegin() + i;
}

int64_t GSA::read_id_of_suffix(int64_t const i) {
  return i - offset(i);
}

int64_t GSA::offset(int64_t const i) {
  if (i < 0 || i >= concat.size()) return -1;
  int64_t offset{-1};
  for(string::const_iterator it = concat.cbegin() + i; 
      it > concat.cbegin() && *it != TERM_CHAR; --it) ++offset;
  return offset;
}

int64_t GSA::len(int64_t const i) {
  if (i < 0 || i >= concat.size()) return -1;
  int64_t len = 1;
  for (string::const_iterator it = concat.cbegin() + i; it < concat.cend() && *it != TERM_CHAR; ++it) ++len;
  return len;
}

GSA::~GSA() {
  std::free(sa);
}

// DEBUG FUNCS
void GSA::print_pos() {
  for (int64_t * it = sa; it < (sa+sa_sz); it++) {
    cout << *it << "\n";
  }
  cout << endl;
}

void GSA::print_gsa() {
  for (int64_t *it = sa; it < (sa+sa_sz); it++) {
    for(string::const_iterator s = suffix_at(*it); 
        s < concat.cend() && *s != TERM_CHAR; s++) {
      cout << *s;
    }
    cout << "$\n";
  }
}
void GSA::print_reads() {
  size_t base{0}, top{0};
  while ((top = concat.find(TERM_CHAR, base)) != string::npos) {
    cout << concat.substr(base, top-base+1) << endl;
    base = top+1;
  }
}
void GSA::print_phreds() {
  size_t base{0}, top{0};
  while ((top = phred.find(PHRED_DELIM, base)) != string::npos) {
    cout << phred.substr(base, top-base+1) << endl;
    base = top+1;
  }
}


void GSA::print_concat() {
  for(string::const_iterator it = concat.cbegin(); it < concat.cend(); it++) {
    cout << *it;
  }
  cout << endl;
}
