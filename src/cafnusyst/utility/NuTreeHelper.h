#pragma once

#include <iostream>
#include <string>
// ROOT
#include "TTree.h"
#include "TBranch.h"

// anaobj
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"

namespace cafnusyst{

class NuTreeHelper{

public:

  NuTreeHelper();
  ~NuTreeHelper();

  void Init();
  void Reset();
  void FillVariable(const caf::Proxy<caf::SRTrueInteraction>& nu);

  TTree* GetTree(){ return fTree; }

  bool DoDebug;

private:
  TTree* fTree{nullptr};

  // Variables
  Double_t Enu;
  std::vector<std::vector<double>> rws;
  Int_t Mode;


};

} // END namespace cafnusyst
