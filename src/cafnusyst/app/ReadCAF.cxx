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
#include "cafnusyst/interface/LoopRecord.h"

template class caf::Proxy<caf::StandardRecord>;

namespace cliopts {
  std::string input_filename = "";
  std::string output_filename = "";
  std::string envvar = "FHICL_FILE_PATH";
  std::string fhicl_key = "generated_systematic_provider_configuration";
  size_t NMax = std::numeric_limits<size_t>::max();
  size_t NSkip = 0;
  bool DoDebug = false;
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
               "\t--debug          : Run debugging mode.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::fhicl_key = argv[++opt];
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::input_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::output_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "--debug") {
      cliopts::DoDebug = true;
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
  if (!cliopts::input_filename.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

  std::ifstream inputFile(cliopts::input_filename);
  if(!inputFile.is_open()){
    printf("[ERROR] %s does not exist\n", cliopts::input_filename);
    return 1;
  }

  std::string filePath;
  cafnusyst::LoopRecord LR(
    "",
    "cafTree", "rec",
    "globalTree", "global",
    "genieEvt", "genie_record"
  ); 
  LR.SetNMaxCAFEventsToProcess(cliopts::NMax);
  if(cliopts::DoDebug) LR.DoDebug = true;

  printf("Saving outputs to: %s\n", cliopts::output_filename.c_str());
  TFile *fOutputFile = new TFile(cliopts::output_filename.c_str(), "RECREATE");
  fOutputFile->cd();
  LR.nuTreeHelper->Init();

  // Loop over input files
  while (std::getline(inputFile, filePath)) {
    printf("[Input] %s\n", filePath.c_str());
    LR.ProcessFile(filePath.c_str());
  }
  inputFile.close();

  fOutputFile->cd();
  LR.nuTreeHelper->GetTree()->Write("NuTree");

  fOutputFile->Close();

}
