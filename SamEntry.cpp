#include <string>
#include <vector>
#include <map>
#include <boost/any.hpp> // used to make map with different Value types

#include "string.h"
#include "SamEntry.h"

using namespace std;


// COMPULSORY SAM FIELDS
const int SamEntry::HDR   = 0;
const int SamEntry::FLAG  = 1;
const int SamEntry::RNAME = 2;
const int SamEntry::POS   = 3;
const int SamEntry::MAPQ  = 4;
const int SamEntry::CIGAR = 5;
const int SamEntry::RNEXT = 6;
const int SamEntry::PNEXT = 7;
const int SamEntry::TLEN  = 8;
const int SamEntry::SEQ   = 9;
const int SamEntry::QUAL  = 10;

// BWA OPT FIELDS
const int SamEntry::NM = 11;
const int SamEntry::MD = 12;
const int SamEntry::AS = 13;
const int SamEntry::BC = 14;
const int SamEntry::X0 = 15;
const int SamEntry::X1 = 16;
const int SamEntry::XN = 17;
const int SamEntry::XM = 18;
const int SamEntry::XO = 19;
const int SamEntry::XG = 20;
const int SamEntry::XT = 21;
const int SamEntry::XA = 22;
const int SamEntry::XS = 23;
const int SamEntry::XF = 24;
const int SamEntry::XE = 25;

// ICSMuFin FIELDS
const int SamEntry::LEFT_OHANG = 26;
const int SamEntry::RIGHT_OHANG = 27;
const int SamEntry::BLOCK_ID = 28;

SamEntry::SamEntry(string const& entry) {
  vector<string> split_result;
  split_string(entry, "\t", split_result);

try {
    string header = split_result[0];  // load header for later parsing

    // Load sam fields, then erase from field list
    for (int i = 1; i < 11; i++) {
      if (i == FLAG || i == POS || i == MAPQ) {
        fields[i] = stoi(split_result[i]);
      }
      else {
        fields[i] = split_result[i];
      }
    }

    // parse the fastq header field
    vector<string> header_subs;
    // split off the cancer cns
    split_string(header, "[", header_subs);
    fields[HDR] = header_subs[0];
  
    // cut of the end "]" char
    header = header_subs[1].substr(0, header_subs[1].length()-1);
  
    // split up the remaining fields
    header_subs.clear();
    split_string(header, ";", header_subs);
    fields[LEFT_OHANG] = stoi(header_subs[0]);
    fields[RIGHT_OHANG] = stoi(header_subs[1]);
    fields[BLOCK_ID] = stoi(header_subs[2]);

    if (fields.size() > 11) {
      // extract bwa optional fields
      split_result.erase(split_result.begin(), split_result.begin() + 11);
  
      fields[NM] = startsWith("NM", split_result);
      fields[MD] = startsWith("MD", split_result);
      fields[AS] = startsWith("AS", split_result);
      fields[BC] = startsWith("BC", split_result);
      fields[X0] = startsWith("X0", split_result);
      fields[X1] = startsWith("X1", split_result);
      fields[XN] = startsWith("XN", split_result);
      fields[XM] = startsWith("XM", split_result);
      fields[XO] = startsWith("XO", split_result);
      fields[XG] = startsWith("XG", split_result);
      fields[XT] = startsWith("XT", split_result);
      fields[XA] = startsWith("XA", split_result);
      fields[XS] = startsWith("XS", split_result);
      fields[XF] = startsWith("XF", split_result);
      fields[XE] = startsWith("XE", split_result);
    }
  }
catch(...) {
  cout << "Exception!!" << endl;
}
}

string SamEntry::startsWith(string const& tok, vector<string> const& fields) {


  for (string const& s : fields) {
    if (tok == s.substr(0, 2)) return s;
  }
  return "\0";
}

void SamEntry::snv_push_back(int v) {
  SNVLocations.push_back(v);
}

int SamEntry::snvLocSize() { return SNVLocations.size();}
int SamEntry::snvLocation(int idx) {return SNVLocations[idx];}
void SamEntry::setSNVLocation(int idx, int val) {SNVLocations[idx] = val;}
bool SamEntry::containsIndel() {
  string cigar = get<string>(SamEntry::CIGAR);
  if (cigar.find('I') || cigar.find('D')) {
    return true;
  }
  return false;
}

