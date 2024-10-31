#include "WeightUpdater.h"
#include "TROOT.h"

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

  DoDebug = false;

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

  // Check POT and Livetime first;
  TH1D *hInputPOT = (TH1D *)f_input->Get( (fBaseDirName+fPOTHistName).c_str() );
  if( !AddPOTHist(hInputPOT) ){
    printf("[WeightUpdater::ProcessFile] Input file does not have POT histogram, skipping:\n");
    printf("[WeightUpdater::ProcessFile] - %s\n", inputfile.c_str());
  }
  TH1D *hInputLivetime = (TH1D *)f_input->Get( (fBaseDirName+fLivetimeHistName).c_str() );
  if( !AddLivetimeHist(hInputLivetime) ){
    printf("[WeightUpdater::ProcessFile] Input file does not have Livetime histogram, skipping:\n");
    printf("[WeightUpdater::ProcessFile] - %s\n", inputfile.c_str());
  }

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

    // now loop over true neutrinos
    for(size_t i_nu=0; i_nu<N_MC; i_nu++){

      if(DoDebug){
        printf("[WeightUpdater::ProcessFile]   - i_nu = %ld\n", i_nu);
      }

      auto& nu = srproxy->mc.nu[i_nu];

#ifdef USE_SBNCAF
      size_t genieIdx = nu.genie_evtrec_idx;

      // 1) Direct access
      //fInputGENIETree->GetEntry(genieIdx);
      // 2) Use GetEntryWithIndex
      std::uint32_t this_tunc_hash = srproxy->hdr.sourceNameHash;
      Long64_t GENIEEntryFromIndex = fInputGENIETree->GetEntryNumberWithIndex(this_tunc_hash, genieIdx);
      fInputGENIETree->GetEntry(GENIEEntryFromIndex);

      printf("[WeightUpdater::ProcessFile]     - (Truncated hash, GENIE Tree position) = (%u, %u) -> GetEntryNumberWithIndex = %ld\n", this_tunc_hash, genieIdx, GENIEEntryFromIndex);
#endif

#ifdef USE_DUNECAF
      size_t genieIdx = nu.genieIdx;
      fInputGENIETree->GetEntry(genieIdx);
#endif

      // Get genie event record
      genie::EventRecord const &GenieGHep = *fInputGENIENtp->event;

      // CAF-to-GENIE matching validation
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

      // Evaluate reweights
      systtools::event_unit_response_w_cv_t resp = fRH->GetEventVariationAndCVResponse(GenieGHep);
      if(resp.size() != NExpectedWeights){
        printf("[WeightUpdater::ProcessFile] resp.size() = %ld but NExpectedWeights = %ld\n", resp.size(), NExpectedWeights);
        abort();
      }

      if(DoDebug){
        printf("[WeightUpdater::ProcessFile]     - => done.\n");
        printf("[WeightUpdater::ProcessFile]     - Now updating weights\n");
      }


      // Upated fSR (caf::StandardRecord*),
      // convert this into FlatRecord using flat::Flat::Fill(const T& x)
#ifdef USE_SBNCAF
      fSR->mc.nu[i_nu].genie_evtrec_idx = GlobalGENIEEventCounter;
#endif
#ifdef USE_DUNECAF
      fSR->mc.nu[i_nu].genieIdx = GlobalGENIEEventCounter;
#endif

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

        // Upated fSR (caf::StandardRecord*),
        // convert this into FlatRecord using flat::Flat::Fill(const T& x)
        fSR->mc.nu[i_nu].wgt.emplace_back();
        for(const auto& w: ws){
          fSR->mc.nu[i_nu].wgt.back().univ.push_back(w);

          if(DoDebug){
            printf("[WeightUpdater::ProcessFile]       - w =  = %f\n", w);
          }

        }

      } // END resp loop

      // Also fill output GENIE tree
      fOutputGENIENtp->Fill(GlobalGENIEEventCounter, &GenieGHep);
      fOutputGENIETree->Fill();

      // Increase genie event counter
      GlobalGENIEEventCounter++;

    } // END nu loop

    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] => MC Loop DONE\n");
    }

#ifdef USE_SBNCAF
    //===========================
    // UPDATE SLC
    //===========================

    const size_t N_SLC = srproxy->slc.size();
    if(DoDebug){
      printf("[WeightUpdater::ProcessFile] - N_SLC = %ld\n", N_SLC);
    }

    for(size_t i_slc=0; i_slc<N_SLC; i_slc++){

      if(DoDebug){
        printf("[WeightUpdater::ProcessFile]   - i_slc = %ld\n", i_slc);
      }

      int this_slc_truth_index = srproxy->slc[i_slc].truth.index;
      if(DoDebug){
        printf("[WeightUpdater::ProcessFile]     - truth index = %d\n", this_slc_truth_index);
      }

      if(this_slc_truth_index>=0){
        fSR->slc[i_slc].truth = fSR->mc.nu[this_slc_truth_index];
      }

    }
#endif

    fOutputFlatSR->Clear();
    fOutputFlatSR->Fill(*fSR);
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

  fOutputGENIETree = new TTree(fGENIETreeName.c_str(), fGENIETreeName.c_str());
  fOutputGENIETree->Branch(fGENIERecName.c_str(), &fOutputGENIENtp);

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

#ifdef USE_SBNCAF
    printf("[WeightUpdater::CreateGlobalTree] - Number of Parameter sets = %d\n", input_srglobal->wgts.size());
    for(unsigned int i = 0; i < input_srglobal->wgts.size(); ++i){
      const caf::SRWeightPSet& pset = input_srglobal->wgts[i];
      std::cout << "  " << i << ": " << pset.name << ", type " << pset.type << ", " << pset.nuniv << " universes, adjusted parameters:" << std::endl;
      for(const caf::SRWeightMapEntry& entry: pset.map){
        std::cout << "    " << entry.param.name << std::endl;
      }

      srglobal.wgts.push_back( pset );
#endif
#ifdef USE_DUNECAF
    printf("[WeightUpdater::CreateGlobalTree] - Number of Parameter sets = %d\n", input_srglobal->wgts.params.size());
    for(unsigned int i = 0; i < input_srglobal->wgts.params.size(); ++i){
      const caf::SRSystParamHeader& pset = input_srglobal->wgts.params[i];
      std::cout << "  " << i << ": " << pset.name << ", id = " << pset.id << ", nshifts = " << pset.nshifts << std::endl;

      srglobal.wgts.params.push_back( pset );
#endif

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

#ifdef USE_SBNCAF
    srglobal.wgts.emplace_back();
#endif
#ifdef USE_DUNECAF
    srglobal.wgts.params.emplace_back();
#endif

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

#ifdef USE_SBNCAF
    // Name
    srglobal.wgts.back().name = fRH->GetSystProvider()[matched_idx_sp]->GetFullyQualifiedName()+"_"+sph.prettyName;
    srglobal.wgts.back().nuniv = sph.isCorrection ? 1 : sph.paramVariations.size();

    printf("[WeightUpdater::CreateGlobalTree] Adding %s to globalTree\n", srglobal.wgts.back().name.c_str());

    // TODO
    srglobal.wgts.back().type = caf::kMultiSim;
#endif
#ifdef USE_DUNECAF
    // Name
    srglobal.wgts.params.back().name = fRH->GetSystProvider()[matched_idx_sp]->GetFullyQualifiedName()+"_"+sph.prettyName;
    srglobal.wgts.params.back().nshifts = sph.isCorrection ? 1 : sph.paramVariations.size();

    printf("[WeightUpdater::CreateGlobalTree] Adding %s to globalTree\n", srglobal.wgts.params.back().name.c_str());
#endif

    // Weight map entry (e.g., dep dials)
    auto it = map_resp_to_respless.find( sph.systParamId );

    if(it!=map_resp_to_respless.end()){

      for(const auto depdialid: it->second){
        const auto& sph_dep = fRH->GetHeader(depdialid);
        std::vector<double> paramVars_dep = sph_dep.isCorrection ? std::vector<double>(1, sph_dep.centralParamValue) : sph_dep.paramVariations;
        std::vector<float> widths_dep { paramVars_dep.begin(), paramVars_dep.end() };

#ifdef USE_SBNCAF
        srglobal.wgts.back().map.emplace_back();

        srglobal.wgts.back().map.back().param.name = sph_dep.prettyName;
        srglobal.wgts.back().map.back().vals = widths_dep;
#endif

      }
    }
    else{
      // single dial

#ifdef USE_SBNCAF
      srglobal.wgts.back().map.emplace_back();

      srglobal.wgts.back().map.back().param.name = sph.prettyName;
      std::vector<double> paramVars = sph.isCorrection ? std::vector<double>(1, sph.centralParamValue) : sph.paramVariations;
      std::vector<float> widths { paramVars.begin(), paramVars.end() };
      srglobal.wgts.back().map.back().vals = widths;
#endif
#ifdef USE_DUNECAF
      srglobal.wgts.params.back().id = int(pid);
#endif

    }

  } // END Loop pid

  fOutputGlobalTree->Fill();

}

void WeightUpdater::SetOutputPOTHistName(std::string name){

  fPOTHistName = name;

  fOutputPOT = new TH1D(fPOTHistName.c_str(), fPOTHistName.c_str(), 1, 0, 1);

}
void WeightUpdater::SetOutputLivetimeHistName(std::string name){

  fLivetimeHistName = name;

  fOutputLivetime = new TH1D(fLivetimeHistName.c_str(), fLivetimeHistName.c_str(), 1, 0, 1);

}

bool WeightUpdater::AddPOTHist(TH1D *h_input){

  if(fOutputPOT && !h_input){
    return false;
  }
  if(!fOutputPOT){
    return true;
  }

  std::cout << "[WeightUpdater::AddPOTHist] POT = " << h_input->GetBinContent(1) << std::endl;

  fOutputPOT->Add(h_input);
  return true;

}

bool WeightUpdater::AddLivetimeHist(TH1D *h_input){
  
  if(fOutputLivetime && !h_input){
    return false;
  }
  if(!fOutputLivetime){
    return true;
  }
  
  fOutputLivetime->Add(h_input);
  return true;

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
  fOutputGENIETree->SetDirectory(OutTDir);

  fOutputFile->Write();

  fOutputFile->Close();

  printf("[WeightUpdater::Save] Done\n");

}

} // END namespace cafnusyst
