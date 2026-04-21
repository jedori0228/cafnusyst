#include "LoopRecord.h"
#include "TROOT.h"
#include "TSystem.h"

namespace cafnusyst{

LoopRecord::LoopRecord(
  std::string basedirname,
  std::string caftreename,
  std::string srname,
  std::string globaltreename,
  std::string srglobalname,
  std::string genietreename,
  std::string genierecname){

  fBaseDirName = basedirname=="" ? "" : basedirname+"/";
  fCAFTreeName = caftreename;
  fSRName = srname;
  fGlobalTreeName = globaltreename;
  fSRGlobalName = srglobalname;
  fGENIETreeName = genietreename;
  fGENIERecName = genierecname;

  NProcessedCAFEvents = 0;
  NMaxCAFEventsToProcess = 0;
  GlobalGENIEEventCounter = 0;
  NProcessedFiles = 0;

  DoDebug = false;

  nuTreeHelper = new NuTreeHelper();

}

LoopRecord::~LoopRecord(){

}

void LoopRecord::SetNMaxCAFEventsToProcess(size_t nmax){
  NMaxCAFEventsToProcess = nmax;
}

void LoopRecord::ProcessFile(std::string inputfile){

  TFile *f_input = TFile::Open(inputfile.c_str());

  // CAF tree

  TTree *fInputCAFTree = (TTree *)f_input->Get( (fBaseDirName+fCAFTreeName).c_str());
  size_t ThisNCAFEvents = fInputCAFTree->GetEntries();

  if(DoDebug){
    printf("[LoopRecord::ProcessFile] ThisNCAFEvents = %ld\n", ThisNCAFEvents);
  }

  caf::SRGlobal* srglobal = nullptr;
  TTree *t_input_global = (TTree *)f_input->Get( (fBaseDirName+fGlobalTreeName).c_str() );
  if(t_input_global){
    t_input_global->SetBranchAddress(fSRGlobalName.c_str(), &srglobal);
    t_input_global->GetEntry(0);
  }

  // - GENIE tree
  TTree *fInputGENIETree = (TTree *)f_input->Get( (fBaseDirName+fGENIETreeName).c_str() );
  genie::NtpMCEventRecord *fInputGENIENtp = nullptr;
  fInputGENIETree->SetBranchAddress(fGENIERecName.c_str(), &fInputGENIENtp);
  size_t ThisNGENIEEvents = fInputGENIETree->GetEntries();

  // - SRProxy to access record
  caf::StandardRecordProxy* srproxy = new caf::StandardRecordProxy(fInputCAFTree, fSRName.c_str());

  // Loop over CAFTree
  for (size_t cafev_it = 0; cafev_it < ThisNCAFEvents; ++cafev_it) {

    if(DoDebug){
      printf("[LoopRecord::ProcessFile] * CAF entry = %ld\n", cafev_it);
    }
    // Check if NMaxCAFEventsToProcess is set
    if( NMaxCAFEventsToProcess>0 ){
      // if set, check if we have reached the maximum
      if(NProcessedCAFEvents>=NMaxCAFEventsToProcess ){
        printf("[LoopRecord::ProcessFile] * Reached the maximum events to process, N_MAX = %ld\n", NMaxCAFEventsToProcess);
        break;
      }
    }

    fInputCAFTree->GetEntry(cafev_it);

    const size_t N_MC = srproxy->mc.nu.size();
    if(DoDebug){
      printf("[LoopRecord::ProcessFile] - N_MC = %ld\n", N_MC);
    }

    // now loop over true neutrinos
    for(size_t i_nu=0; i_nu<N_MC; i_nu++){

      if(DoDebug){
        printf("[LoopRecord::ProcessFile]   - i_nu = %ld\n", i_nu);
      }

      auto& nu = srproxy->mc.nu[i_nu];

      nuTreeHelper->FillVariable(nu);

/*
      // TODO Find the matched GENIE EventRecord from this SRTrueInteraction
      size_t genieIdx = nu.genieIdx;
      fInputGENIETree->GetEntry(genieIdx);

      // Get genie event record
      genie::EventRecord const &GenieGHep = *fInputGENIENtp->event;

      for(size_t i_syst_dials=0; i_syst_dials<nu.syst_dials.size(); i_syst_dials++){
        printf("[LoopRecord::ProcessFile]     - Dial index =  %ld\n", i_syst_dials);
        for(size_t i_rw=0; i_rw<nu.syst_dials[i_syst_dials].weights.size(); i_rw++){
          const auto& this_rw = nu.syst_dials[i_syst_dials].weights[i_rw];
          printf("[LoopRecord::ProcessFile]       - Weight index = %ld, weight = %1.4f\n", i_rw, this_rw.GetValue());
        }
      }
*/

      GlobalGENIEEventCounter++;

    } // END nu loop

    if(DoDebug){
      printf("[LoopRecord::ProcessFile] => MC Loop DONE\n");
    }

    NProcessedCAFEvents++;

  } // END caf event loop

  NProcessedFiles++;

  printf("[LoopRecord::ProcessFile] -----------------\n");
  printf("[LoopRecord::ProcessFile] File done\n");

}

} // END namespace cafnusyst
