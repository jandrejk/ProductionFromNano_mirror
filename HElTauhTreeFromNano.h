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
  bool electronSelection(unsigned int index);
  bool tauSelection(unsigned int index);
  bool diElectronVeto();
  bool pairSelection(unsigned int index);
  /////////////////////////////////////////////////
  
  HElTauhTreeFromNano(TTree *tree=0, bool doSvFit=false, bool correctRecoil=false, bool isMC_ = false, std::vector<edm::LuminosityBlockRange> lumiBlock = std::vector<edm::LuminosityBlockRange>(), std::string prefix="HTTET");
  virtual ~HElTauhTreeFromNano();
  
};

#endif

#ifdef HElTauhTreeFromNano_cxx
HElTauhTreeFromNano::HElTauhTreeFromNano(TTree *tree, bool doSvFit, bool correctRecoil, bool isMC_, std::vector<edm::LuminosityBlockRange> lumiBlocks, std::string prefix) : HTauTauTreeFromNanoBase(tree, doSvFit, correctRecoil, isMC_, lumiBlocks, prefix)
{}

HElTauhTreeFromNano::~HElTauhTreeFromNano()
{}

#endif // #ifdef HElTauhTreeFromNano_cxx
