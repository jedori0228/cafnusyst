#include "WeightUpdater.h"
#include "TROOT.h"
#include "TSystem.h"

namespace cafnusyst{

WeightUpdater::WeightUpdater(
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

  fRH = nullptr;

  NProcessedCAFEvents = 0;
  NMaxCAFEventsToProcess = 0;
  GlobalGENIEEventCounter = 0;
  NProcessedFiles = 0;

  fOutputFile = nullptr;
  fOutputCAFTree = nullptr;
  fOutputGENIETree = nullptr;
  fOutputGlobalTree = nullptr;
  fOutputFlatSR = nullptr;
  fOutputGENIENtp = nullptr;

  NExpectedWeights = 0;

  fOutputPOT = nullptr;
  fOutputLivetime = nullptr;

  CheckCAFToGENIEMatching = false;
  DoDebug = false;

  fWeightsOnly = false;

}

WeightUpdater::~WeightUpdater(){

}

void WeightUpdater::SetResponseHelper(std::string fclname){
  fRH = new nusyst::response_helper(fclname);
}

void WeightUpdater::SetNMaxCAFEventsToProcess(size_t nmax){
  NMaxCAFEventsToProcess = nmax;
}

void WeightUpdater::ProcessFile(std::string inputfile){

  TFile *f_input = TFile::Open(inputfile.c_str());

  // CAF tree

  TTree *fInputCAFTree = (TTree *)f_input->Get( (fBaseDirName+fCAFTreeName).c_str());
  size_t ThisNCAFEvents = fInputCAFTree->GetEntries();

  if(DoDebug){
    printf("[WeightUpdater::ProcessFile] ThisNCAFEvents = %ld\n", ThisNCAFEvents);
  }

  const caf::CAFType caftype = caf::GetCAFType(fInputCAFTree);

  // Access StandardRecord if nested
  caf::StandardRecord* fSR = nullptr;

  if(caftype==caf::kNested){
    fInputCAFTree->SetBranchAddress(fSRName.c_str(), &fSR);

    assert(fInputCAFTree);

    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] * Input is nested caf\n");
    }

  }
  else if(caftype==caf::kFlat){

    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] * Input is flatcaf\n");
    }

    printf("[ERROR] Flatcaf is not supported yet\n");
    abort();

  }
  else{
    printf("[ERROR] Unknown caf type from\nFile: %s\nTree: %s\n", inputfile.c_str(), fCAFTreeName.c_str());
    abort();
  }

  // - Check Global
  if(!fOutputGlobalTree){

    caf::SRGlobal* srglobal = nullptr;

    TTree *t_input_global = (TTree *)f_input->Get( (fBaseDirName+fGlobalTreeName).c_str() );
    if(t_input_global){
      t_input_global->SetBranchAddress(fSRGlobalName.c_str(), &srglobal);
      t_input_global->GetEntry(0);
    }

    CreateGlobalTree(srglobal);

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
      printf("[WeightUpdater::ProcessFile] * CAF entry = %ld\n", cafev_it);
    }
    // Check if NMaxCAFEventsToProcess is set
    if( NMaxCAFEventsToProcess>0 ){
      // if set, check if we have reached the maximum
      if(NProcessedCAFEvents>=NMaxCAFEventsToProcess ){
        printf("[WeightUpdater::ProcessFile] * Reached the maximum events to process, N_MAX = %ld\n", NMaxCAFEventsToProcess);
        break;
      }
    }

    fInputCAFTree->GetEntry(cafev_it);

    //===========================
    // UPDATE NU
    //===========================

    const size_t N_MC = srproxy->mc.nu.size();
    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] - N_MC = %ld\n", N_MC);
    }

    // In weights-only mode, emit a slim record with only syst_dials populated. 
    // Size mc.nu to the input so index alignment is preserved.
    caf::StandardRecord outSR;
    if(fWeightsOnly){
      outSR.mc.nu.resize(N_MC);
    }

    // now loop over true neutrinos
    for(size_t i_nu=0; i_nu<N_MC; i_nu++){

      if(DoDebug){
        printf("[WeightUpdater::ProcessFile]   - i_nu = %ld\n", i_nu);
      }

      auto& nu = srproxy->mc.nu[i_nu];

      // Recomputed weights go to the full input record (full mode)
      // or the slim output record (weights-only mode).
      caf::SRTrueInteraction& outNu = fWeightsOnly ? outSR.mc.nu[i_nu] : fSR->mc.nu[i_nu];

      // TODO Find the matched GENIE EventRecord from this SRTrueInteraction
      size_t genieIdx = nu.genieIdx;
      fInputGENIETree->GetEntry(genieIdx);

      // Get genie event record
      genie::EventRecord const &GenieGHep = *fInputGENIENtp->event;

      // CAF-to-GENIE matching validation
      if(CheckCAFToGENIEMatching){
        genie::GHepParticle *ISLep = GenieGHep.Probe();
        TLorentzVector ISLepP4 = *ISLep->P4();
        double enu_from_genie = ISLepP4.E();
        if( std::fabs(nu.E.GetValue() - enu_from_genie) > 1E-5 ){
          printf("[WeightUpdater::ProcessFile]     - ENu from CAF and GENIE are not close; matching problem\n");
          printf("[WeightUpdater::ProcessFile]       ENu (CAF, GENIE) = (%f, %f), diff = %f > 1E-5\n", nu.E.GetValue(), enu_from_genie, std::fabs(nu.E.GetValue() - enu_from_genie));
          abort();
        }
        if(DoDebug){
          printf("[WeightUpdater::ProcessFile]     - ENu (Proxy) = %f\n", nu.E.GetValue());
          printf("[WeightUpdater::ProcessFile]     - ENu from GENIE = %f\n", enu_from_genie);
          printf("[WeightUpdater::ProcessFile]     - Running responses..\n");
        }
      }

      // Evaluate reweights
      systtools::event_unit_response_w_cv_t resp = fRH->GetEventVariationAndCVResponse(GenieGHep);
      if(resp.size() != NExpectedWeights){
        printf("[WeightUpdater::ProcessFile] resp.size() = %ld but NExpectedWeights = %ld\n", resp.size(), NExpectedWeights);
        abort();
      }

      if(DoDebug){
        printf("[WeightUpdater::ProcessFile]     - => done.\n");
        printf("[WeightUpdater::ProcessFile]     - Current number size of syst_dials = %ld\n", outNu.syst_dials.size());
        printf("[WeightUpdater::ProcessFile]     - Now updating weights\n");
      }

      // It is possible that we are processing multiple input CAFs
      // Then the genieIdx should be updated using the current number of EventRecord read.
      // genieIdx refers into the re-indexed GENIE tree, which is not written in
      // weights-only mode; there we leave it at its default (-1).
      if(!fWeightsOnly){
        outNu.genieIdx = GlobalGENIEEventCounter;
      }
      // Start from an empty syst_dials so weights from a previous spill do not linger.
      outNu.syst_dials.clear();
      for(const auto& v: resp){
        const systtools::paramId_t& pid = v.pid;
        const double& CVw = v.CV_response;
        const std::vector<double>& ws = v.responses;

        systtools::SystParamHeader const &sph = fRH->GetHeader( pid );
        if(DoDebug){
          printf("[WeightUpdater::ProcessFile]     - Param name = %s\n", sph.prettyName.c_str());
          printf("[WeightUpdater::ProcessFile]       - ParamID = %ld\n", pid);
        }
        if(sph.isResponselessParam){
          if(DoDebug){
            printf("[WeightUpdater::ProcessFile]       - Responsless parameter, so skipping\n");
          }
          continue;
        }

        // Upated record (caf::StandardRecord),
        // convert this into FlatRecord using flat::Flat::Fill(const T& x)
        outNu.syst_dials.emplace_back();
        for(const auto& w: ws){
          if(DoDebug){
            printf("[WeightUpdater::ProcessFile]       - w =  = %f\n", w);
          }
          outNu.syst_dials.back().weights.push_back(w);
        }

      } // END resp loop

      // Also fill output GENIE tree (not emitted in weights-only mode)
      if(!fWeightsOnly){
        fOutputGENIENtp->Fill(GlobalGENIEEventCounter, &GenieGHep);
        fOutputGENIETree->Fill();
        GlobalGENIEEventCounter++;
      }


    } // END nu loop

    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] => MC Loop DONE\n");
    }

    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] => Current fSR->mc.nu.size() = %ld\n", fSR->mc.nu.size());
    }

    fOutputFlatSR->Clear();
    fOutputFlatSR->Fill(fWeightsOnly ? outSR : *fSR);
    fOutputCAFTree->Fill();

    NProcessedCAFEvents++;

  } // END caf event loop

  NProcessedFiles++;

  printf("[WeightUpdater::ProcessFile] -----------------\n");
  printf("[WeightUpdater::ProcessFile] File done\n");

}

void WeightUpdater::SetOutputFileName(std::string FileName){

  fOutputFile = new TFile(FileName.c_str(), "RECREATE");

  fOutputCAFTree = new TTree(fCAFTreeName.c_str(), fCAFTreeName.c_str());
  fOutputFlatSR = new caf::FlatStandardRecord(fOutputCAFTree, fSRName.c_str(), "", 0);

  // The GENIE tree only exists to feed nusystematics on the input side; a
  // weights-only friend file does not carry it.
  if(!fWeightsOnly){
    fOutputGENIETree = new TTree(fGENIETreeName.c_str(), fGENIETreeName.c_str());
    fOutputGENIETree->Branch(fGENIERecName.c_str(), &fOutputGENIENtp);
  }

  CreateMetadataTree();

}

// TODO
void WeightUpdater::CreateMetadataTree(){

  fOutputFile->cd();
  fOutputFile->mkdir("metadata");
  fOutputFile->cd("metadata");

  TTree *fMetadataTree = new TTree("metatree", "metatree");

  std::string fMetadata_key;
  std::string fMetadata_value;

  fMetadataTree->Branch("key", &fMetadata_key);
  fMetadataTree->Branch("value", &fMetadata_value);

  fMetadataTree->Write();

  fOutputFile->cd();

}

void WeightUpdater::CreateGlobalTree(caf::SRGlobal* input_srglobal){

  if(!fRH){
    printf("[WeightUpdater::CreateGlobalTree] Response helper is not set. Run WeightUpdater::SetResponseHelper()\n");
    abort();
  }

  fOutputFile->cd();
  fOutputGlobalTree = new TTree(fGlobalTreeName.c_str(), fGlobalTreeName.c_str());
  caf::SRGlobal srglobal = caf::SRGlobal();
  fOutputGlobalTree->Branch(fSRGlobalName.c_str(), &srglobal);

  if(input_srglobal){
    // Copying from input SRGlobal
    printf("[WeightUpdater::CreateGlobalTree] @@ Copying input SRGlobal\n");

    printf("[WeightUpdater::CreateGlobalTree] - Number of Parameter sets = %d\n", input_srglobal->wgts.params.size());
    for(unsigned int i = 0; i < input_srglobal->wgts.params.size(); ++i){
      const caf::SRSystParamHeader& pset = input_srglobal->wgts.params[i];

      srglobal.wgts.params.push_back( pset );
    }

  }

  // Now adding new weights

  // make a map of responsless-response params
  std::map<systtools::paramId_t, std::vector<systtools::paramId_t>> map_resp_to_respless;
  for(systtools::paramId_t pid : fRH->GetParameters()) {
    systtools::SystParamHeader const &sph = fRH->GetHeader(pid);
    if(sph.isResponselessParam){
      auto it = map_resp_to_respless.find( sph.responseParamId );
      if( it != map_resp_to_respless.end() ){
        it->second.push_back( sph.systParamId );
      }
      else{
        map_resp_to_respless[sph.responseParamId] = {};
        map_resp_to_respless[sph.responseParamId].push_back( sph.systParamId );
      }
    }
  }

  for(systtools::paramId_t pid : fRH->GetParameters()) {
    systtools::SystParamHeader const &sph = fRH->GetHeader(pid);

    if(sph.isResponselessParam){
      if(DoDebug) {
        printf("[WeightUpdater::CreateGlobalTree] Responsless dial found: %s, thus skipping\n", sph.prettyName.c_str());
      }
      continue;
    }

    srglobal.wgts.params.emplace_back();

    NExpectedWeights++;

    // Find the IGENIESystProvider_tool(ISystProviderTool) for this pid
    int matched_idx_sp = -1;
    for(int idx_sp=0; idx_sp<fRH->GetSystProvider().size(); idx_sp++){
      if(fRH->GetSystProvider()[idx_sp]->ParamIsHandled(pid)){
        matched_idx_sp = idx_sp;
      }
    }
    if(matched_idx_sp<0){
      printf("[WeightUpdater::CreateGlobalTree] IGENIESystProvider_tool not found from pid = %d\n", int(pid));
      abort();
    }

    // Name
    srglobal.wgts.params.back().name = fRH->GetSystProvider()[matched_idx_sp]->GetFullyQualifiedName()+"_"+sph.prettyName;

    printf("[WeightUpdater::CreateGlobalTree] Adding %s to globalTree\n", srglobal.wgts.params.back().name.c_str());

    // Weight map entry (e.g., dep dials)
    auto it = map_resp_to_respless.find( sph.systParamId );

    if(it!=map_resp_to_respless.end()){

      for(const auto depdialid: it->second){
        const auto& sph_dep = fRH->GetHeader(depdialid);
        std::vector<double> paramVars_dep = sph_dep.isCorrection ? std::vector<double>(1, sph_dep.centralParamValue) : sph_dep.paramVariations;
        std::vector<float> widths_dep { paramVars_dep.begin(), paramVars_dep.end() };
      }

    }
    else{
      // single dial

      srglobal.wgts.params.back().id = int(pid);

      std::vector<double> paramVars = sph.isCorrection ? std::vector<double>(1, sph.centralParamValue) : sph.paramVariations;
      srglobal.wgts.params.back().vals.clear();
      for(const auto& v: paramVars) srglobal.wgts.params.back().vals.push_back(v);

    }

  } // END Loop pid

  fOutputGlobalTree->Fill();

}

void WeightUpdater::Save(){

  printf("[WeightUpdater::Save] Saving output\n");

  fOutputFile->cd();

  if(fBaseDirName!=""){
    fOutputFile->mkdir( fBaseDirName.substr(0, fBaseDirName.size() - 1).c_str());
  }
  TDirectory *OutTDir = fBaseDirName=="" ? fOutputFile : (TDirectory *)fOutputFile->Get(fBaseDirName.substr(0, fBaseDirName.size() - 1).c_str());

  fOutputGlobalTree->SetDirectory(OutTDir);
  fOutputCAFTree->SetDirectory(OutTDir);
  if(fOutputGENIETree) fOutputGENIETree->SetDirectory(OutTDir);

  fOutputFile->Write();
  fOutputFile->Close();

  printf("[WeightUpdater::Save] Done\n");

}

} // END namespace cafnusyst
