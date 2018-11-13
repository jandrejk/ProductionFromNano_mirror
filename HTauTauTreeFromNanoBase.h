/*****************************
* M. Bluj, NCBJ, Poland
* 6/12/2017
*****************************/

#ifndef HTauTauTreeFromNanoBase_h
#define HTauTauTreeFromNanoBase_h

/*Include header of base class providing interface to the NanoAOD Events tree 
  as produced by nanoEventsChain->MakeClass("NanoEventsSkeleton") */
#undef NanoEventsSkeleton_cxx //undefine to protect againist problems with multile implementation
#include "NanoEventsSkeleton.h"
#include "EventWriter.h"
// #include "ParameterConfig.cc"

//#include <TROOT.h>
//#include <TChain.h>
//#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"

#include "HTTEvent.h"
#include <string.h>
#include <vector>
#include <iostream>
#include "json.hpp"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"

using json = nlohmann::json;

class HTauTauTreeFromNanoBase : public NanoEventsSkeleton {

public :

  /////////////////////////////////////////////////
  /// Lorentz vector
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PolarLorentzVector;
  /// Trigger data struct
  struct TriggerData {
    std::string path_name;
    unsigned int leg1Id=0,         leg2Id=0;//0-undefined, 11-electron, 13-muon, 15-tau
    int          leg1BitMask=0,    leg2BitMask=0;//definition depends on Id, cf. PhysicsTools/NanoAOD/python/triggerObjects_cff.py
    float        leg1Pt=0.,        leg2Pt=0.;
    float        leg1L1Pt=0.,      leg2L1Pt=0.;
    float        leg1Eta=0.,       leg2Eta=0.;
    float        leg1OfflinePt=0., leg2OfflinePt=0.;
  };

  virtual void initHTTTree(const TTree *tree, std::string prefix="HTT");
  void debugWayPoint(std::string description, std::vector<double> dbls = {}, std::vector<int> ints = {}, vector<string> descr = {""});
  void fillEvent(unsigned int bestPairIndex = 9999);
  virtual bool buildPairs();
  virtual void fillPairs(unsigned int bestPairIndex);
  virtual void fillJets(unsigned int bestPairIndex);
  virtual void fillLeptons();
  virtual void fillGenLeptons();
  void applyMetRecoilCorrections(HTTPair &aPair);
  virtual bool thirdLeptonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, int leptonPdg, double dRmin=-1);
  virtual bool extraMuonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin=-1);
  virtual bool extraElectronVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin=-1);
  int muonSelection(HTTParticle aLepton);
  int electronSelection(HTTParticle aLepton);
  int tauSelection(HTTParticle aLepton);
  bool failsGlobalSelection();
  virtual bool pairSelection(unsigned int index);
  virtual unsigned int bestPair(std::vector<unsigned int> &pairIndexes);
  void computeSvFit(HTTPair &aPair);
  TLorentzVector runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
			      const TVector2 &aMET, const TMatrixD &covMET);
  bool jetSelection(unsigned int index, unsigned int bestPairIndex);
  int getGenMatch(unsigned int index, std::string colType="");
  int getGenMatch(TLorentzVector selObj);
  //  int getTriggerMatching(unsigned int index, bool checkBit=false, std::string colType="");
  int getTriggerMatching(unsigned int index, TLorentzVector p4_1, bool checkBit=false, std::string colType="");
  int getMetFilterBits();
  double getZPtReweight(const TLorentzVector &genBosonP4, bool doSUSY=false);
  bool isGoodToMatch(unsigned int ind);
  TLorentzVector getGenComponentP4(std::vector<unsigned int> &indexes, unsigned int iAbsCharge);
  bool eventInJson();

  template<typename T> T getBranchValue(const char *branchAddress, unsigned int index);
  Double_t getProperty(std::string name, unsigned int index, std::string colType="");
  Double_t getProperty(std::string name, unsigned int index, TLorentzVector obj, std::string colType="");
  std::vector<Double_t> getProperties(const std::vector<std::string> & propertiesList, unsigned int index, std::string colType="");
  std::vector<Double_t> getProperties(const std::vector<std::string> & propertiesList, unsigned int index, TLorentzVector obj, std::string colType="");

  Int_t getFilter(std::string name);
  std::vector<Int_t> getFilters(const std::vector<std::string> & propertiesList);

  void writeJECSourceHeader(const std::vector<string> &jecSources);
  void writePropertiesHeader(const std::vector<std::string> & propertiesList);
  void writeTriggersHeader(const std::vector<TriggerData> &triggerBits);
  void writeFiltersHeader(const std::vector<std::string> &filterBits);

  void initJecUnc(std::string correctionFile);
  double getJecUnc(unsigned int index, unsigned int isrc, bool up=true);
  std::map<string, double> getValuesAfterJecSplitting(unsigned int iJet);

  static bool compareLeptons(const HTTParticle& i, const HTTParticle& j);
  static bool comparePairs(const HTTPair& i, const HTTPair& j);
  //int isGenPartDaughterPdgId(int index, unsigned int aPdgId);
  //bool isGenPartDaughterIdx(int index, int mother);
  bool getDirectDaughterIndexes(std::vector<unsigned int> & indexes, unsigned int motherIndex, bool ignoreNeutrinos=true);
  unsigned int findFinalCopy(unsigned int index);
  unsigned int findFirstCopy(unsigned int index);
  int genTauDecayMode(std::vector<unsigned int> &daughterIndexes);
  bool findBosonP4(TLorentzVector &bosonP4, TLorentzVector &visBosonP4);
  bool findTopP4(TLorentzVector &topP4, TLorentzVector &antiTopP4);

  std::vector<HTTPair> httPairCollection, httPairs_;
  HTTJetCollection         httJetCollection;
  std::vector<HTTParticle> httLeptonCollection;
  std::vector<HTTParticle> httGenLeptonCollection;
  std::vector<TriggerData> triggerBits_;
  std::vector<std::string> filterBits_;
  TTree *t_TauCheck;
  TTree *httTree;

  std::unique_ptr<EventWriter> evtWriter;
  std::unique_ptr<TFile> httFile;
  std::unique_ptr<HTTEvent> httEvent;
  TH1D* puweights_histo;
  TH1F* hStats;
  TH2D* zptmass_histo;  

  std::unique_ptr<ClassicSVfit> svFitAlgo_;
  std::unique_ptr<RecoilCorrector> recoilCorrector_;
  std::unique_ptr<TFile> zPtReweightFile, zPtReweightSUSYFile, puweights;
  TLorentzVector p4SVFit, p4Leg1SVFit, p4Leg2SVFit;   

  std::vector<edm::LuminosityBlockRange> jsonVector;

  bool firstWarningOccurence_; // used to print warnings only at first occurnece in the event loop
  bool isMC;
  bool isSync;
  bool applyRecoil;
  int passMask_;
  unsigned int check_event_number;
  unsigned int bestPairIndex_;


  std::vector<std::string> leptonPropertiesList, genLeptonPropertiesList, jecUncertList;
  std::vector<JetCorrectionUncertainty*> jecUncerts;

  json Settings;

  template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
  }

  HTauTauTreeFromNanoBase(TTree *tree=0, std::vector<edm::LuminosityBlockRange> lumiBlocks = std::vector<edm::LuminosityBlockRange>() , string prefix="HTT");
  virtual ~HTauTauTreeFromNanoBase();
  virtual Int_t    Cut(Long64_t entry);
  virtual void     Loop(Long64_t nentries_max=-1, unsigned int sync_event=-1);
};

#endif

