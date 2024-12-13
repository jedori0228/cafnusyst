// std
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// ROOT
#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"
// systematicstools
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"

// anaobj
#ifdef USE_SBNCAF
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/Flat/FlatRecord.h"
#endif
#ifdef USE_DUNECAF
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"
#endif


namespace cliopts {
  std::string input_filename = "";
  size_t NMax = std::numeric_limits<size_t>::max();
  size_t NSkip = 0;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-s <NSkip>       : Number of events to skip.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::input_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = systtools::str2T<size_t>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {

  HandleOpts(argc, argv);
  if (!cliopts::input_filename.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

#ifdef USE_SBNCAF
  std::string fCAFTreeName = "recTree";
  std::string fSRName = "rec";
  std::string fGlobalTreeName = "globalTree";
  std::string fSRGlobalName = "global";
#endif
#ifdef USE_DUNECAF
  std::string fCAFTreeName = "cafmaker/cafTree";
  std::string fSRName = "rec";
  std::string fGlobalTreeName = "cafmaker/globalTree";
  std::string fSRGlobalName = "global";
#endif

  // Input file
  TFile *f_input = new TFile(cliopts::input_filename.c_str());

  // Global tree
  TTree *t_input_global = (TTree *)f_input->Get(fGlobalTreeName.c_str());
  caf::SRGlobal* srglobal = nullptr;
  t_input_global->SetBranchAddress(fSRGlobalName.c_str(), &srglobal);
  t_input_global->GetEntry(0);

#ifdef USE_SBNCAF
  printf("[ReadCAF] - Number of Parameter sets = %d\n", srglobal->wgts.size());
  for(unsigned int i = 0; i < srglobal->wgts.size(); ++i){
    const caf::SRWeightPSet& pset = srglobal->wgts[i];
    std::cout << "  " << i << ": " << pset.name << ", type " << pset.type << ", " << pset.nuniv << " universes, adjusted parameters:" << std::endl;
    for(const caf::SRWeightMapEntry& entry: pset.map){
      std::cout << "    " << entry.param.name << std::endl;
    }
  }
#endif
#ifdef USE_DUNECAF
  printf("[ReadCAF] - Number of Parameter sets = %d\n", srglobal->wgts.params.size());
  for(unsigned int i = 0; i < srglobal->wgts.params.size(); ++i){
    const caf::SRSystParamHeader& pset = srglobal->wgts.params[i];
    std::cout << "  " << i << ": " << pset.name << ", id = " << pset.id << ", nshifts = " << pset.nshifts << std::endl;
  }
#endif

  // CAF tree
  std::cout << "[ReadCAF] fCAFTreeName = " << fCAFTreeName << std::endl;
  TTree *fInputCAFTree = (TTree *)f_input->Get(fCAFTreeName.c_str());
  size_t ThisNCAFEvents = fInputCAFTree->GetEntries();
  size_t NToRead = std::min(ThisNCAFEvents, cliopts::NMax);
  printf("@@ Number of CAF events = %ld\n", ThisNCAFEvents);
  printf("@@ Number of CAF events to read = %ld\n", NToRead);
 
  const caf::CAFType caftype = caf::GetCAFType(fInputCAFTree);

  // - SRProxy to access record
  caf::StandardRecordProxy* srproxy = new caf::StandardRecordProxy(fInputCAFTree, fSRName.c_str());

  // Loop over CAFTree
  unsigned int NProcessedCAFEvents = 0;
  for (size_t cafev_it = 0; cafev_it < NToRead; ++cafev_it) {

    printf("[ReadCAF] * CAF entry = %ld\n", cafev_it);
    if( cliopts::NMax>0 ){
      // if set, check if we have reached the maximum
      if(NProcessedCAFEvents>=NToRead){
        printf("[ReadCAF] * Reached the maximum events to process, N_MAX = %ld\n", cliopts::NMax);
        break;
      }
    }

    fInputCAFTree->GetEntry(cafev_it);

    const size_t N_MC = srproxy->mc.nu.size();
    printf("[ReadCAF]   * N_MC = %ld\n", N_MC);

    // now loop over true neutrinos
    for(size_t i_nu=0; i_nu<N_MC; i_nu++){

      auto& nu = srproxy->mc.nu[i_nu];

#ifdef USE_DUNECAF
      printf("[ReadCAF]     * i_nu = %ld\n", i_nu);
      printf("[ReadCAF]       * nu.E() = %f\n", nu.E.GetValue());
      printf("[ReadCAF]       * nu.wgt.size() = %ld\n", nu.wgt.size());
      for(size_t i_wgt=0; i_wgt<nu.wgt.size(); i_wgt++){
        const auto& srmult = nu.wgt[i_wgt];
        printf("[ReadCAF]         * i_wgt = %ld, srmult.univ.size() = %ld\n", i_wgt, srmult.univ.size());
        for(size_t i_univ=0; i_univ<srmult.univ.size(); i_univ++){
          printf("[ReadCAF]           * i_univ = %ld, weight = %f\n", i_univ, srmult.univ[i_univ].GetValue());
        }
      }
#endif

    }

    NProcessedCAFEvents++;

  }

}
