// cigarParser.h
// Author: Izaak Coleman

#ifndef CIGAR_PARSER
#define CIGAR_PARSER
#include <string> 

class CigarParser {
  /* Provides array access to cigar elements as the main functional interface. 
     Two functions in the main functional interface:
     CigarParser::length_at(i) -> returns number of bases assigned to operation i.
     CigarParser::operation_at(i) -> returns char representing operation i.
     e.g cp = CigarParser("80M12I8M");
     cout << cp.operation_at(1) << endl;
     cout << cp.length_at(1) << endl;
     
     will print to stdout:
     I
     12
  */

private: 
  std::string CIGAR; // The CIGAR STRING
  int64_t sz; // # operations in CIGAR
public:

  CigarParser(std::string const& cig);
  /*
     Takes CIGAR as input.
  */

  CigarParser(std::string const & cig, bool const reverse);
  /*
     Takes CIGAR as input. If reverse == true, then the CIGAR string
     is saved in reverse; this is useful when analysing reverse complements
     of aligned fastqs. 
  */

  int64_t length_at(int64_t const i) const;
  /* 
     Returns number of bases assigned to operation i.
  */

  char operation_at(int64_t const i) const;
  /*
     Returns char representing operation i.
  */

  int64_t size() const;
  /*
     Returns number of operations in CIGAR;
  */
};
#endif
