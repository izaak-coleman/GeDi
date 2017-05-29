/*
main.cpp
Author: Izaak Coleman
*/


#include <iostream>
#include <vector>
#include <string>
#include "boost/program_options.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

// benchmarking
#include <sys/time.h>
#include <sys/resource.h>
// testing
#include <typeinfo>

#include "SuffixArray.h"
#include "SNVIdentifier.h"
#include "util_funcs.h"
#include "Reads.h"
#include "GenomeMapper.h"
#include "benchmark.h"


using namespace std;
// Return values
static const int ERROR_IN_COMMAND_LINE = 1; 
static const int SUCCESS = 0; 
static const int ERROR_UNHANDLED_EXCEPTION = 2;

static const int CALIB = 4;
static const int BASE33_CONVERSION = 33;
static const int PHRED_LBOUND = 0;
static const int PHRED_UBOUND = 42;

// Default option values
static const int    GSA1_MCT               = 1;
static const int    GSA2_MCT               = 4; 
static const int    MAX_LOW_CONFIDENCE_POS = 10;
static const int    ECONT                  = 0;
static const int    MIN_PHRED_QUAL         = 22; 
static const int    MIN_MAPQ               = 42;
static const double ALLELE_FREQ_OF_ERR     = 0.1; 
static const string OUTPUT_PATH            = "./"; 


int main(int argc, char** argv) 
{ 
  try { 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help", "Print options description.\n") 
      ("gsa1_mct,1",  po::value<int>()->default_value(GSA1_MCT), 
       "Minimum group size extracted from the first GSA. Integer.\n")

      ("gsa2_mct,2",  po::value<int>()->default_value(GSA2_MCT), 
       "Minimum group size extracted from the second GSA. Integer.\n")

      ("min_phred,h", po::value<int>()->default_value(MIN_PHRED_QUAL),
       "Minimum allowed phred score of any character that contributes to a consensus sequence. Base-33 phred score. Integer ranged [0-42]\n")

      ("max_allele_freq_of_error,f", po::value<double>()->default_value(ALLELE_FREQ_OF_ERR), 
       "Maximum allelic frequency of a base within an aligned block that is considered an error frequency. Real number ranged [0-1].\n")
      
      ("max_low_confidence_positions,l",
       po::value<int>()->default_value(MAX_LOW_CONFIDENCE_POS),
       "Maximum number of low confidence positions allowed per block before the block is discarded. Integer.\n")

      ("min_mapq,q", po::value<int>()->default_value(MIN_MAPQ),
       "Defines the minimum Bowtie2 MAPQ score of the aligned healthy consensus sequence of a consensus pair that permits the consensus pair to undergo further analysis. Integer ranged [0-42].\n")

      ("expected_contamination,e", po::value<double>()->default_value(ECONT),
       "Proportion of cancer data set expected to contain healthy derived reads. Real number ranged [0-1].\n")

      ("output_path,p", po::value<string>()->default_value(OUTPUT_PATH), 
       "Path specifying write location of <output_prefix>.SNV_results <output_prefix>.fastq and <output_prefix>.out\n")

      ("expected_coverage,v", po::value<int>()->required(), 
       "Expected coverage of input datasets. Integer. Required.\n")

      ("n_threads,t", po::value<int>()->required(), 
       "Number of threads. Integer. Required.\n")

      ("chromosome,c", po::value<string>()->required(), 
       "Target chromosome for SNV calling. Only sam entries with an RNAME = <chromosome> will be extracted from the sam file. Required.")

      ("input_files,i", po::value<string>()->required(), 
       "Path and name of file containing the input file list. Required.\n")

      ("bt2-idx,x", po::value<string>()->required(), 
       "Basename of the index for the reference genome. Specified value should be identical to Bowtie2's -x option. Required.\n")

      ("output_basename,o", po::value<string>()->required(), 
       "Basename for output SNV call file (SNV_results), consensus fastq and sam alignment. Required.\n");

 
    po::variables_map vm; 
    try { 
      po::store(po::parse_command_line(argc, argv, desc),vm); // can throw 
      if (vm.count("help")) { 
        cout.precision(1);
        std::cout << "***********************************************************************" << std::endl
                  << "* Generalized Suffix Array based Direct Comparison (GeDi) SNV caller. *" << std::endl
                  << "***********************************************************************" << std::endl
                  << std::endl << std::endl;
        std::cout << "usage: [-1 1_arg] [-2 2_arg] [-h h_arg] [-f f_arg] [-e e_arg]"
                  << " [-p p_arg] -v v_arg -t t_arg -c c_arg -i i_arg -x x_arg -o o_arg" 
                  << std::endl;
        std::cout << desc 
                  << std::endl; 
        return SUCCESS; 
      } 
      // Throw if required options are missing.
      po::notify(vm);

      // Check supplied options are valid.
      if (vm["min_phred"].as<int>() < PHRED_LBOUND || vm["min_phred"].as<int>() > PHRED_UBOUND) {
        std::cerr << "ERROR: " 
                  << "--min_phred must take integer value in range [0-42]."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["gsa1_mct"].as<int>() < 1) {
        std::cerr << "ERROR: " 
                  << "--gsa1_mct must be at least 1."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["gsa2_mct"].as<int>() < 1) {
        std::cerr << "ERROR: " 
                  << "--gsa2_mct must be at least 1."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["max_allele_freq_of_error"].as<double>() < 0 ||
          vm["max_allele_freq_of_error"].as<double>() > 1) {
        std::cerr << "ERROR: " 
                  << "--max_allele_freq_of_error must a real number in range [0-1]."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["expected_contamination"].as<double>() < 0 ||
          vm["expected_contamination"].as<double>() > 1) {
        std::cerr << "ERROR: " 
                  << "--expected_contamination must a real number in range [0-1]."

                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["expected_coverage"].as<int>() < 0) {
        std::cerr << "ERROR: " 
                  << "--expected_coverage must be at least 0."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["n_threads"].as<int>() < 1) {
        std::cerr << "ERROR: " 
                  << "--n_threads must be at least 1."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["max_low_confidence_positions"].as<int>() < 0) {
        std::cerr << "ERROR: " 
                  << "--max_low_confidence_positions must be at least 0."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["min_mapq"].as<int>() < PHRED_LBOUND || 
          vm["min_mapq"].as<int>() > PHRED_UBOUND) {
        std::cerr << "ERROR: " 
                  << "--min_mapq must be an integer value in range [0-42]."
                  << std::endl << std::endl
                  << "Refer to --help for input desciption." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }

      struct stat info_p, info_i;
      if (stat(vm["output_path"].as<string>().c_str(), &info_p) != 0) {
        std::cerr << "ERROR: "
                  << "Cannot access " << vm["output_path"].as<string>()
                  << std::endl << "Please check path exists / is correct."
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      else if (!(info_p.st_mode & S_IFDIR)) {
        std::cerr << "ERROR: "
                  << vm["output_path"].as<string>() 
                  << " does not exist. Please make required directories"
                  << std::endl << "before running GeDi." 
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }

      if (stat(vm["input_files"].as<string>().c_str(), &info_i) != 0) {
        std::cerr << "ERROR: "
                  << "Cannot access " << vm["input_files"].as<string>()
                  << std::endl << "Please check path to file / filename."
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      else if (!(info_p.st_mode & S_IFDIR)) {
        std::cerr << "ERROR: "
                  << vm["input_files"].as<string>() 
                  << " cannot be found. Please check path to file/filename."
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
    } 
    catch(po::error& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << desc << std::endl
                << "Program terminating." << std::endl;
      return ERROR_IN_COMMAND_LINE; 
    } 
    // Run GeDi
    ReadPhredContainer  reads(vm["input_files"].as<string>());

    SuffixArray SA(reads, reads.getMinSuffixSize(), vm["n_threads"].as<int>());

    SNVIdentifier snvId(SA, reads, 
                         vm["min_phred"].as<int>()+BASE33_CONVERSION,
                         vm["gsa1_mct"].as<int>(),
                         vm["gsa2_mct"].as<int>(),
                         vm["expected_coverage"].as<int>()*CALIB,
                         vm["n_threads"].as<int>(),
                         vm["max_low_confidence_positions"].as<int>(),
                         vm["expected_contamination"].as<double>(),
                         vm["max_allele_freq_of_error"].as<double>());

    GenomeMapper mapper(snvId, reads,
                        vm["output_path"].as<string>(),
                        vm["output_basename"].as<string>(),
                        vm["chromosome"].as<string>(),
                        vm["bt2-idx"].as<string>(),
                        vm["min_mapq"].as<int>());
    return SUCCESS;
  } 
  catch(std::exception& e) 
  { 
    std::cerr << "Unhandled exception reached main: " 
              << e.what() << ", application will now exit" << std::endl; 
    return ERROR_UNHANDLED_EXCEPTION; 
  } 
  return SUCCESS; 
}
