#pragma once

#include <iostream>
#include <string>
// ROOT
#include "TChain.h"
#include "TBranch.h"
// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"

// anaobj
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"

// nusystematics
#include "nusystematics/utility/response_helper.hh"
// cafnusyst
#include "cafnusyst/utility/Utilities.h"

namespace cafnusyst{

class WeightUpdater{

public:

  WeightUpdater(
    std::string basedirname,
    std::string caftreename,
    std::string srname,
    std::string globaltreename,
    std::string srglobalname,
    std::string genietreename,
    std::string genierecname
  );
  ~WeightUpdater();

  nusyst::response_helper* fRH;
  void SetResponseHelper(std::string fclname);

  std::string fBaseDirName;
  std::string fCAFTreeName;
  std::string fSRName;
  std::string fGlobalTreeName;
  std::string fSRGlobalName;
  std::string fGENIETreeName;
  std::string fGENIERecName;
  size_t NProcessedCAFEvents;
  size_t NMaxCAFEventsToProcess;
  void SetNMaxCAFEventsToProcess(size_t nmax);
  size_t GlobalGENIEEventCounter;
  void ProcessFile(std::string inputfile);
  size_t NProcessedFiles;

  TFile *fOutputFile;
  TTree *fOutputGlobalTree;
  TTree *fOutputWeightTree;
  std::uint32_t fSourceFileHash;
  unsigned int fGenieEventCounter;
  std::vector<std::vector<double>> fWeights;

  void SetOutputFileName(std::string FileName);
  void CreateMetadataTree(); // TODO
  void CreateGlobalTree(caf::SRGlobal* input_srglobal);
  unsigned int NExpectedWeights;

  void Save();

  bool DoDebug;

};

} // END namespace cafnusyst
