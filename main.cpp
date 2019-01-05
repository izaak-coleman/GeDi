/*
main.cpp
Author: Izaak Coleman
*/


#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "boost/program_options.hpp"
#include <sys/stat.h>

// benchmarking
#include <sys/resource.h>
#include <chrono>

#include "SNVIdentifier.h"
#include "GenomeMapper.h"
#include "gsa.h"


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
static const int    PRI_MSS                = 2;
static const int    AUX_MSS                = 4; 
static const int    MAX_SNPS               = 5;
static const int    ECONT                  = 0;
static const int    MIN_PHRED_QUAL         = 35; 
static const int    MIN_MAPQ               = 0;
static const double ALLELE_FREQ_OF_ERR     = 0.1; 
static const string OUTPUT_PATH            = "./"; 
static const string CHR                    = "";
static const string EMFIL                  = "on";


int main(int argc, char** argv) 
{ 
  try { 
    namespace po = boost::program_options; 
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help", "Print options description.") 
      ("primary_mss",  po::value<int>()->default_value(PRI_MSS), 
       "Minimum suffix size of sections extracted from the primary " 
       "suffix array. Integer ranged [1, maxint].\n")

      ("emfilter", po::value<string>()->default_value(EMFIL),
       "Apply emfilter to input tumour data.(on|off)\n")

      ("auxiliary_mss",  po::value<int>()->default_value(AUX_MSS), 
       "Minimum suffix size of sections extracted from the auxiliary GSA. " 
       "Integer ranged [1, maxint].\n")

      ("min_phred,h", po::value<int>()->default_value(MIN_PHRED_QUAL),
       "Minimum phred score of base contributing to a phred-filtered " 
       "consensus sequence Base-33 phred score. Integer ranged [0,42]\n")

      ("max_allele_freq_of_error,f", po::value<double>()->default_value(ALLELE_FREQ_OF_ERR), 
       "All candidate variants with an allele frequency equal to or below "
        "this value will be considered sequencing errors and discarded. " 
        "Real number ranged [0,1]\n")
      
      ("max_SNPs,l",
       po::value<int>()->default_value(MAX_SNPS),
       "Max number of SNPs in consensus sequence allowed before multi-locus " 
       "filter will be triggered, resulting in the consensus pair being " 
       "discarded. Integer [0, maxint]. \n")

      ("min_mapq,q", po::value<int>()->default_value(MIN_MAPQ),
       "Aligned control consensus sequences with a Bowtie2 given MAPQ score " 
       "below this value will be discarded. Integer ranged [0,42] \n.")

      ("expected_contamination,e", po::value<double>()->default_value(ECONT),
       "Proportion of cancer dataset expected to contain healthy-tissue " 
       "derived reads. Real number ranged [0,1].\n")

      ("chromosome,c", po::value<string>()->default_value(CHR), 
       "For targeted sequencing data analysis. Only control consensus "  
       "consensus sequenceswith RNAME identical to the given value "  
       "will be considered for further analysis. Leaving this unspecified "  
       "removes constraint, allowing GeDi to analyse WGS datasets. \n")

      ("expected_coverage,v", po::value<int>()->required(), 
       "Expected coverage of input datasets. Integer range [0,maxint]. Required.\n")

      ("n_threads,t", po::value<int>()->required(), 
       "Number of threads. Integer. Required ranged [1,maxint].\n")


      ("input_files,i", po::value<string>()->required(), 
       "Path and name of file containing list of input fastq files. Required.\n")

      ("bt2-idx,x", po::value<string>()->required(), 
       "Basename of Bowtie2 reference genome index. Specified value " 
       "should be identical to Bowtie2's -x option if running Bowtie2 " 
       "from GeDi. Required.\n")

      ("output_basename,o", po::value<string>()->required(), 
       "Basename for SNV_results, fastq.gz, sam files output by GeDi. Required.\n");

 
    po::variables_map vm; 
    try { 
      po::store(po::parse_command_line(argc, argv, desc),vm); 
      if (vm.count("help")) { 
        cout.precision(1);
        std::cout << "***********************************************************************" << std::endl
                  << "* Generalized Suffix Array based Direct Comparison (GeDi) SNV caller. *" << std::endl
                  << "***********************************************************************" << std::endl
                  << std::endl << std::endl;
        std::cout << "usage: [optional parameters] -v coverage -t n_threads -i fastq_list -x bowtie_index -o output_file_basename"
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
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["emfilter"].as<string>() != "on" &&
          vm["emfilter"].as<string>() != "off") {
        std::cerr << "ERROR: " 
                  << "--emfilter must take value 'on' or 'off'"
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["primary_mss"].as<int>() < 1) {
        std::cerr << "ERROR: " 
                  << "--primary_mss must be at least 1."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["auxiliary_mss"].as<int>() < 1) {
        std::cerr << "ERROR: " 
                  << "--auxiliary_mss must be at least 1."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["max_allele_freq_of_error"].as<double>() < 0 ||
          vm["max_allele_freq_of_error"].as<double>() > 1) {
        std::cerr << "ERROR: " 
                  << "--max_allele_freq_of_error must a real number in range [0-1]."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["expected_contamination"].as<double>() < 0 ||
          vm["expected_contamination"].as<double>() > 1) {
        std::cerr << "ERROR: " 
                  << "--expected_contamination must a real number in range [0-1]."

                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["expected_coverage"].as<int>() < 0) {
        std::cerr << "ERROR: " 
                  << "--expected_coverage must be at least 0."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["n_threads"].as<int>() < 1) {
        std::cerr << "ERROR: " 
                  << "--n_threads must be at least 1."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["max_SNPs"].as<int>() < 0) {
        std::cerr << "ERROR: " 
                  << "--max_SNPs must be at least 0."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      if (vm["min_mapq"].as<int>() < PHRED_LBOUND || 
          vm["min_mapq"].as<int>() > PHRED_UBOUND) {
        std::cerr << "ERROR: " 
                  << "--min_mapq must be an integer value in range [0-42]."
                  << std::endl << std::endl
                  << "Read --help for required and default parameters." << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }

      struct stat info_p, info_i;
      string output_dir = get_dir_from_filename(vm["output_basename"].as<string>());
      if (stat(output_dir.c_str(), &info_p) != 0) {
        std::cerr << "ERROR: "
                  << "Cannot access " << output_dir
                  << std::endl << "Please check path is correct."
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      else if (!(info_p.st_mode & S_IFDIR)) {
        std::cerr << "ERROR: "
                  << vm["output_dir"].as<string>() 
                  << " cannot be found. Please check directory exists."
                  << std::endl << "before running GeDi." 
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }

      if (stat(vm["input_files"].as<string>().c_str(), &info_i) != 0) {
        std::cerr << "ERROR: "
                  << "Cannot access " << vm["input_files"].as<string>()
                  << std::endl << "Please check filename is correct."
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
      else if (!(info_p.st_mode & S_IFDIR)) {
        std::cerr << "ERROR: "
                  << vm["input_files"].as<string>() 
                  << " cannot be found. Please check filename is correct."
                  << std::endl
                  << "Program terminating." << std::endl;
        return ERROR_IN_COMMAND_LINE;
      }
    } 
    catch(po::error& e) 
    { 
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl; 
      std::cerr << "Read --help for required and default parameters." << std::endl
                << "Program terminating." << std::endl;
      return ERROR_IN_COMMAND_LINE; 
    } 

    // Run GeDi
    struct rusage gedi_rss; getrusage(RUSAGE_SELF, &gedi_rss);
    auto gedi_runtime_start = std::chrono::steady_clock::now();
    cout << "Running GeDi with " << vm["n_threads"].as<int>() << " threads." << endl;
    GSA gsa(vm["input_files"].as<string>(), vm["n_threads"].as<int>(),
            vm["bt2-idx"].as<string>(), vm["emfilter"].as<string>());
    SNVIdentifier snvId(gsa, 
                         vm["output_basename"].as<string>(),
                         vm["min_phred"].as<int>()+BASE33_CONVERSION,
                         vm["primary_mss"].as<int>(),
                         vm["auxiliary_mss"].as<int>(),
                         vm["expected_coverage"].as<int>()*CALIB,
                         vm["n_threads"].as<int>(),
                         vm["max_SNPs"].as<int>(),
                         vm["expected_contamination"].as<double>(),
                         vm["max_allele_freq_of_error"].as<double>());

    GenomeMapper mapper(snvId, 
                        vm["output_basename"].as<string>(),
                        vm["chromosome"].as<string>(),
                        vm["bt2-idx"].as<string>(),
                        vm["min_mapq"].as<int>());
    auto gedi_runtime_end = std::chrono::steady_clock::now();
    double gedi_runtime =
      std::chrono::duration_cast<std::chrono::milliseconds>(gedi_runtime_end -
          gedi_runtime_start).count();
    gedi_runtime = gedi_runtime / 60000; // minutes
    getrusage(RUSAGE_SELF, &gedi_rss);
    cout << std::fixed;
    cout << std::setprecision(2);
    cout << endl << endl;
    cout << "GeDi terminated successfully." << endl
         << "Runtime (mins): " << gedi_runtime << endl
         << "Memory usage (bytes) : " << gedi_rss.ru_maxrss << endl;
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
