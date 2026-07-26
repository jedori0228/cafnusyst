#include "NuTreeHelper.h"

namespace cafnusyst{

NuTreeHelper::NuTreeHelper(){

  DoDebug = false;

}

NuTreeHelper::~NuTreeHelper(){

}

void NuTreeHelper::Init(){

  if(DoDebug) std::cout << "[NuTreeHelper::Init()] Called" << std::endl;

  if(DoDebug) std::cout << "[NuTreeHelper::Init()] Creating TTree" << std::endl;

  fTree = new TTree("NuTree", "NuTree");

  if(DoDebug) std::cout << "[NuTreeHelper::Init()] Define branches" << std::endl;

  fTree->Branch("Enu", &Enu, "Enu/D");
  fTree->Branch("rws", &rws);
  fTree->Branch("Mode", &Mode, "Mode/I");

  if(DoDebug) std::cout << "[NuTreeHelper::Init()] Done!" << std::endl;

}

void NuTreeHelper::Reset(){

  Enu = -999;
  rws.clear();
  Mode = -999;

}

void NuTreeHelper::FillVariable(const caf::Proxy<caf::SRTrueInteraction>& nu){

  Reset();

  if(DoDebug) std::cout << "[NuTreeHelper::FillVariable] Called" << std::endl;

  Enu = nu.E;

  for(const auto& syst_dial: nu.syst_dials){
    rws.push_back({});
    for(const auto& rw: syst_dial.weights){
      rws.back().push_back(rw);
    }
  }

  Mode = int(nu.mode);

  if(DoDebug) std::cout << "[NuTreeHelper::FillVariable] Now running fTree->Fill().." << std::endl;

  fTree->Fill();

  if(DoDebug) std::cout << "[NuTreeHelper::FillVariable] Done!" << std::endl;


}

} // END namespace cafnusyst
