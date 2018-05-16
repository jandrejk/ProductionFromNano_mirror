#ifndef HElTauhTreeFromNano_h
#define HElTauhTreeFromNano_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "HTTEvent.h"
#include <iostream>

// Header of the base class
#include "HTauTauTreeFromNanoBase.h"

class HElTauhTreeFromNano : public HTauTauTreeFromNanoBase {
 public :

  /////////////////////////////////////////////////
  /// ET final state specific
  bool diElecVeto();
  bool pairSelection(unsigned int index);
  /////////////////////////////////////////////////
  
  HElTauhTreeFromNano(TTree *tree=0, bool doSvFit=false, bool correctRecoil=false, std::vector<std::string> lumis = std::vector<std::string>(), std::string prefix="HTTET");
  virtual ~HElTauhTreeFromNano();
  
};

#endif

#ifdef HElTauhTreeFromNano_cxx
HElTauhTreeFromNano::HElTauhTreeFromNano(TTree *tree, bool doSvFit, bool correctRecoil, std::vector<std::string> lumis, std::string prefix) : HTauTauTreeFromNanoBase(tree, doSvFit, correctRecoil, lumis, prefix)
{}

HElTauhTreeFromNano::~HElTauhTreeFromNano()
{}

#endif // #ifdef HElTauhTreeFromNano_cxx
