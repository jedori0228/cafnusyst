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
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"

// nusystematics
#include "nusystematics/utility/response_helper.hh"
// cafnusyst
#include "cafnusyst/utility/Utilities.h"
#include "cafnusyst/utility/ResourceMonitor.h"

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
  TTree *fOutputCAFTree;
  TTree *fOutputGENIETree;
  TTree *fOutputGlobalTree;
  caf::FlatStandardRecord* fOutputFlatSR;
  genie::NtpMCEventRecord *fOutputGENIENtp;

  void SetOutputFileName(std::string FileName);
  void CreateMetadataTree(); // TODO
  void CreateGlobalTree(caf::SRGlobal* input_srglobal);
  unsigned int NExpectedWeights;

  void SetOutputPOTHistName(std::string name);
  void SetOutputLivetimeHistName(std::string name);
  std::string fPOTHistName;
  std::string fLivetimeHistName;
  TH1D *fOutputPOT;
  TH1D *fOutputLivetime;
  bool AddPOTHist(TH1D *h_input);
  bool AddLivetimeHist(TH1D *h_input);

  void Save();

  bool CheckCAFToGENIEMatching;
  bool DoDebug;

  bool fWeightsOnly;

  // When true, evaluate and print resource usage (wall/cpu time, peak RSS)
  // for the configuring-response_helper, creating-GlobalTree, and
  // evaluating+saving-reweights steps. Off by default since getrusage()
  // calls add a small overhead of their own. Set this (e.g. from a
  // --monitor CLI option) before calling SetResponseHelper()/ProcessFile()
  // so it takes effect for all monitored steps.
  bool DoMonitor;

  // Resource usage (wall/cpu time, peak RSS) accumulated across every
  // "evaluate + save reweights" step (one per SRTrueInteraction). Reported
  // once, in Save(), rather than per-neutrino to avoid flooding stdout.
  // Disabled by default; ProcessFile() syncs its enabled state with
  // DoMonitor on every call.
  cafnusyst::ResourceAccumulator fReweightResourceAcc{"Evaluating+saving reweights (per nu)", false};

};

} // END namespace cafnusyst
