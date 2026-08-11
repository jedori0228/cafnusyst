// std
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// fhiclcpp
#include "fhiclcpp/ParameterSet.h"
// systematicstools
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"
// nusystematics
#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/response_helper.hh"
// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
// ROOT
#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"

// cafnusyst
#include "cafnusyst/interface/WeightUpdater.h"

template class caf::Proxy<caf::StandardRecord>;

namespace cliopts {
  std::string fclname = "";
  std::string input_filename = "";
  std::string output_filename = "";
  std::string envvar = "FHICL_FILE_PATH";
  std::string fhicl_key = "generated_systematic_provider_configuration";
  size_t NMax = std::numeric_limits<size_t>::max();
  bool NMaxSet = false;
  size_t NSkip = 0;
  bool DoDebug = false;
  bool DoMonitor = false;
  bool WeightsOnly = false;
  std::string sr_branch = "rec";
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help          : Show this message.\n"
               "\t-c <config.fcl>    : fhicl file to read.\n"
               "\t-k <list key>      : fhicl key to look for parameter headers,\n"
               "\t                     "
               "\"generated_systematic_provider_configuration\"\n"
               "\t                     by default.\n"
               "\t-i <inputlist.txt> : List of input CAF files\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-s <NSkip>       : Number of events to skip.\n"
               "\t-o <out.root>    : File to write validation canvases to.\n"
               "\t--weights-only   : Emit a slim record with only the weights populated \n"
               "\t                   (for use as a friend of the input CAF).\n"
               "\t                   Incompatible with -N.\n"
               "\t--sr-branch <n>  : StandardRecord branch name in the output,\n"
               "\t                   \"rec\" by default.\n"
               "\t--debug          : Run debugging mode.\n"
               "\t--monitor        : Evaluate and print resource usage (wall/cpu time,\n"
               "\t                   peak RSS) for configuring the response_helper,\n"
               "\t                   creating the GlobalTree, and evaluating+saving\n"
               "\t                   reweights.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-c") {
      cliopts::fclname = argv[++opt];
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::fhicl_key = argv[++opt];
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::input_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = systtools::str2T<size_t>(argv[++opt]);
      cliopts::NMaxSet = true;
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::output_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "--weights-only") {
      cliopts::WeightsOnly = true;
      ++opt;
    } else if (std::string(argv[opt]) == "--sr-branch") {
      cliopts::sr_branch = argv[++opt];
    } else if (std::string(argv[opt]) == "--debug") {
      cliopts::DoDebug = true;
      ++opt;
    } else if (std::string(argv[opt]) == "--monitor") {
      cliopts::DoMonitor = true;
      ++opt;
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_laconic.xml"); // quiet mode

  HandleOpts(argc, argv);
  if (!cliopts::fclname.size()) {
    std::cout << "[ERROR]: Expected to be passed a -c option." << std::endl;
    SayUsage(argv);
    return 1;
  }
  if (!cliopts::input_filename.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }
  if (cliopts::WeightsOnly && cliopts::NMaxSet) {
    std::cout << "[ERROR]: -N (event cap) is incompatible with --weights-only: "
                 "a truncated output cannot stay entry-aligned with the parent CAF."
              << std::endl;
    return 1;
  }

  std::ifstream inputFile(cliopts::input_filename);
  if(!inputFile.is_open()){
    printf("[ERROR] %s does not exist\n", cliopts::input_filename.c_str());
    return 1;
  }

  std::string filePath;
  cafnusyst::WeightUpdater wu(
    "",
    "cafTree", cliopts::sr_branch,
    "globalTree", "global",
    "genieEvt", "genie_record"
  );
  wu.fWeightsOnly = cliopts::WeightsOnly;
  if(cliopts::DoDebug) wu.DoDebug = true;
  if(cliopts::DoMonitor) wu.DoMonitor = true;
  wu.SetOutputFileName(cliopts::output_filename);
  wu.SetNMaxCAFEventsToProcess(cliopts::NMax);
  wu.SetResponseHelper(cliopts::fclname);

  // Loop over input files
  while (std::getline(inputFile, filePath)) {
    printf("[Input] %s\n", filePath.c_str());
    wu.ProcessFile(filePath.c_str());
  }
  inputFile.close();

  wu.Save();

}
