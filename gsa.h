#ifndef GSA_H
#define GSA_H

#include <cinttypes>
#include <string>
#include <vector>


class GSA {
  /* gsa takes the input healthy and tumour data, concatenates all the 
     healthy and tumour reads of the form conatenation = HT, where H, T
     are the concatenated healthy, tumour reads respectively. A suffix
     array is generated from the concatenation. By recording a sorted
     array of the positions of each read in the concatenation (binary
     search arrays), the read from which a given suffix derived can be 
     calculated. This enables : 
      - retrieval of a read that a given suffix in the concatenation derives
        from
      - calculation of the offset of a suffix from its read
      - calculation of the length of suffix
     Furthermore, the index at which concatenation transitions from H to T
     data is recorded (tsi: tumour start index), allowing the tissue type
     of a given suffix to be derived by equality comparison against tsi.
    */

  int N_THREADS;
  int64_t sa_sz; // size of suffix array
  int64_t tsi; // index in concat at which data transitions from H to T
  std::string concat; // = HT
  std::string phred;
  int64_t * sa;
  int64_t max_read_len;

  // Functions
private:
  void constructConcat(std::string const & fname);
  void load_fq_data(std::string const & fname);
  bool good_quality(std::string const & q);
  void add_fq_data(std::string const & fq, std::string const & q);
  void read_header(std::string const & header_fname, std::vector<std::string> & h_fnames, 
                   std::vector<std::string> & t_fnames);
  void split_string(std::string const & s, std::string const& tokens, std::vector<std::string>
      &split_strings);

  void remove_short_suffixes(int64_t min_suffix_length); 
  // removes suffixes of len() < min_suffix_length

  void xorSwap(int64_t *x, int64_t *y); // swap *x with *y

  int64_t bubbleRemove(int64_t * const a, int64_t const sz, int64_t const invalid);
  void em_filter(std::string const & ref, std::string const & f);

  void reportFilesDataset(std::vector<std::string> const & file_list, 
                          bool const tissue);
  // Reports to stdout whether a files was loaded as tumour or control data

public:
  GSA(std::string const & header_fname, int t, std::string const & ref, std::string const & em_filtering); // inits all class members
  ~GSA(); // free sa

  bool tissuetype(int64_t const i); // returns TUMOUR/HEALTHY
  std::string::const_iterator suffix_at(int64_t i);
  std::string::const_iterator phred_at (int64_t i);
  std::string::const_iterator read_of_suffix(int64_t const i);
  int64_t offset(int64_t const i);
  int64_t len(int64_t const i);
  int64_t get_max_read_len();
  int64_t get_min_suf_size();
  int64_t read_id_of_suffix(int64_t const i);
  int64_t sa_element(int64_t pos);
  int64_t size();
  std::string get_suffix_string(int64_t i);
  std::string get_phred_string(int64_t i);
  // debug funcs
  void print_gsa();
  void print_pos();
  void print_reads();
  void print_concat();
  void print_phreds();
};
#endif
