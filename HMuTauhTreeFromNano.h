#ifndef HMuTauhTreeFromNano_h
#define HMuTauhTreeFromNano_h

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

class HMuTauhTreeFromNano : public HTauTauTreeFromNanoBase {
 public :

  /////////////////////////////////////////////////
  /// MT final state specific
  bool diMuonVeto();
  bool pairSelection(unsigned int index);
  bool muonSelection(unsigned int index);
  bool tauSelection(unsigned int index);
  /////////////////////////////////////////////////
  
  HMuTauhTreeFromNano(TTree *tree=0, bool doSvFit=false, bool correctRecoil=false, bool isMC_ = false, std::vector<edm::LuminosityBlockRange> lumiBlocks = std::vector<edm::LuminosityBlockRange>(), std::string prefix="HTTMT");
  virtual ~HMuTauhTreeFromNano();
  
};

#endif

#ifdef HMuTauhTreeFromNano_cxx
HMuTauhTreeFromNano::HMuTauhTreeFromNano(TTree *tree, bool doSvFit, bool correctRecoil, bool isMC_, std::vector<edm::LuminosityBlockRange> lumiBlocks, std::string prefix) : HTauTauTreeFromNanoBase(tree, doSvFit, correctRecoil, isMC_, lumiBlocks, prefix)
{}

HMuTauhTreeFromNano::~HMuTauhTreeFromNano()
{}

#endif // #ifdef HMuTauhTreeFromNano_cxx
