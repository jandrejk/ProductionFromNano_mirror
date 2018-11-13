#ifndef WarsawAnalysis_HTTDataFormats_HTTEvent_h
#define WarsawAnalysis_HTTDataFormats_HTTEvent_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH2F.h"
#include "TMatrixDEigen.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooFunctor.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TBits.h"
#include <map>
#include <vector>
#include <bitset>
#include <iostream>
#include <string.h>

#include "PropertyEnum.h"
#include "JecUncEnum.h"
#include "TriggerEnum.h"
#include "FilterEnum.h"
#include "SelectionBitsEnum.h"
#include "AnalysisEnums.h"
#include "LeptonCuts.h"
#include "utils/BTagCalibration/interface/BTagCalibrationStandalone.h"
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTEvent{


 public:

  static std::map<string,PropertyEnum> usePropertyFor;

  enum sampleTypeEnum {DUMMY, MuonData, EleData, TauData , DY, DYLowM, WJets, TTbar, ST, Diboson, EWK, h, H, A};

  static const int ntauIds = 13;
  static const int againstMuIdOffset = 0;
  static const int againstEIdOffset = againstMuIdOffset+2;
  static const int mvaIsoIdOffset = againstEIdOffset+5;
  static const TString tauIDStrings[ntauIds];//implementation in cxx

  HTTEvent(){clear();}
  ~HTTEvent(){}
  ///Data member setters.

  void setRun(unsigned int x){runId = x;}

  void setEvent(unsigned long int x){eventId = x;}

  void setLS(unsigned long int x){lsId = x;}

  void setNPU(float x){nPU = x;}

  void setPUWeight(float x){puWeight = x;}

  void setNPV(unsigned int x){nPV = x;}

  void setRho(float x){rho = x;}

  void setMCatNLOWeight(float x){aMCatNLOweight = x;}

  void setMCWeight(float x){mcWeight = x;}

  void setStage1Cat(int x){ htxs_stage1cat = x; }

  void setXsec(float x){xsec = x;}

  void setGenNEventsWeight(float x){genNEvents = x;}

  void setTopP4(const TLorentzVector &p4, const TLorentzVector &antiP4) {topP4 = p4; antiTopP4 = antiP4; isSetTopP4 = true; }

  void setTopPtReWeight(float x){topPtReWeight = x;}

  void setTopPtReWeightR1(float x){topPtReWeightR1 = x;}

  void setGenBosonP4(const TLorentzVector &p4, const TLorentzVector &visP4) {bosP4 = p4; bosVisP4 = visP4; isSetGenBosonP4 = true; }

  void setZPtReWeight(float x){zPtReWeight = x;}

  void setZPtReWeightSUSY(float x){zPtReWeightSUSY = x;}

  void setLHE_Ht(float x){lheHt = x;}

  void setNNLO_ggH_weight(double x){NNLO_ggH_weight = x;}

  void setLHEnOutPartons(int x){lheNOutPartons = x;}

  void setSampleType(string sampletype);

  void setDecayModeMinus(int x){decayModeMinus = x;}

  void setDecayModePlus(int x){decayModePlus = x;}

  void setDecayModeBoson(int x){decayModeBoson = x;}

  void setGenPV(const TVector3 & aPV) {genPV = aPV;}

  void setAODPV(const TVector3 & aPV) {AODPV = aPV;}

  void setRefittedPV(const TVector3 & aPV) {refittedPV = aPV;}

  void setIsRefit(bool aBit){isRefit = aBit;};

  void setNTracksInRefit(const int & nTracks) {nTracksInRefit = nTracks;};

  void setSelectionBit(SelectionBitsEnum iBit, bool value = true) {selectionWord.SetBitNumber((int)iBit, value);}

  void setMET(const TVector2 &aVector) {met = aVector;}

  void setMET_uncorr(const TVector2 &aVector) {met_uncorr = aVector;}

  void setMETFilterDecision(unsigned int aMETFilterDecision) {metFilterDecision = aMETFilterDecision;}

  void setFilters(const std::vector<Int_t> & aFilters) { filters = aFilters;}

  ////////////////////////

  ///Reset class data members
  void clear();

  void clearSelectionWord() {selectionWord.ResetAllBits();}

  ///Data member getters.
  unsigned int getRunId() const {return runId;}

  unsigned long int getEventId() const {return eventId;}

  unsigned long int getLSId() const {return lsId;}

  float getNPU() const {return nPU;}

  float getPUWeight() const {return puWeight;}

  unsigned int getNPV() const {return nPV;}

  float getRho() const {return rho;}

  float getMCatNLOWeight() const {return aMCatNLOweight;}

  double getTopPtReWeight(bool run_1=true) const {return run_1 ? topPtReWeightR1 : topPtReWeight;}

  double getZPtReWeight() const { return zPtReWeight;}

  double getNNLO_ggH_weight() const {return NNLO_ggH_weight;}

  float getMCWeight() const {return mcWeight;}

  int getStage1Cat() const {return htxs_stage1cat;}

  float getXsec() const {return xsec;}

  float getGenNEventsWeight() const {return genNEvents;}

  float getLHE_Ht() const {return lheHt;}

  int getLHEnOutPartons() const {return lheNOutPartons;}

  sampleTypeEnum getSampleType() const {return sampleType;}

  int getDecayModeMinus() const {return decayModeMinus;}

  int getDecayModePlus() const {return decayModePlus;}

  int getDecayModeBoson() const {return decayModeBoson;}

  TLorentzVector getGenBosonP4(bool visP4=false) const { return visP4 ? bosVisP4 : bosP4 ; }

  TLorentzVector getTopP4(bool anti=false) const { return anti ? antiTopP4 : topP4 ; }

  TVector2 getMET() const {return met;}

  TVector2 getMET_uncorr() const {return met_uncorr;}

  const TVector3 & getGenPV() const {return genPV;}

  const TVector3 & getAODPV() const {return AODPV;}

  const TVector3 & getRefittedPV() const {return refittedPV;}

  bool getIsRefit() const {return isRefit;};

  int getNTracksInRefit() const {return nTracksInRefit;}

  bool checkSelectionBit(SelectionBitsEnum iBit) const {return selectionWord.TestBitNumber((unsigned int)iBit);}

  unsigned int getMETFilterDecision() const { return metFilterDecision;}

  int getFilter(FilterEnum index) const {return (unsigned int)index<filters.size()?  filters[(unsigned int)index]: -999;}



 private:
  ///Event run, run and LS number

  unsigned int runId;
  unsigned long int eventId, lsId;

  //Generator event weight
  float mcWeight;
  float xsec;
  float genNEvents;
 
  ///Weight used to modify the pt shape.
  float topPtReWeight, topPtReWeightR1;
  float zPtReWeight,   zPtReWeightSUSY;

  // WG1 NNLO ggH reweighting: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/SignalModelingTools#Gluon_Fusion_NNLOPS_reweighting
  double NNLO_ggH_weight;

  ///Ht value from LHE record.
  float lheHt;

  ///Number of outgoing partons from LHE record
  int lheNOutPartons;
  int htxs_stage1cat = 0;

  ///MCatNLO weight
  float aMCatNLOweight;

  ///Number of true PU vertices from MC
  float nPU;

  float puWeight;

  //Number of reocnstructed PV
  unsigned int nPV;

  ///PU energy density with FastJet, rho
  float rho;

  ///Type of the physics process or DATA
  sampleTypeEnum sampleType;

  ///Boson (H, Z, W) decay mode
  int decayModeBoson;

  ///Boson (H, Z, W) p4 and visible p4
  bool isSetGenBosonP4;
  TLorentzVector bosP4, bosVisP4;

  ///top and antitop p4
  bool isSetTopP4;
  TLorentzVector topP4, antiTopP4;

  ///Tau decay modes
  int decayModeMinus, decayModePlus;

  ///Primary Vertices recontructed with different methods
  //Generated PV position
  TVector3 genPV;

  //PV stored in miniAOD
  TVector3 AODPV;

  ///PV recontructed from PF candidates, refitted
  TVector3 refittedPV;

  ///Flag marking if refit was successfull
  bool isRefit;

  ///Number of tracks used in the refit
  int nTracksInRefit;

  ///Bit word coding event selection result
  TBits selectionWord;

  //MET vector, uncorrected
  TVector2 met_uncorr;

  //MET vector, recoil corr
  TVector2 met;

  //MET filter decision
  unsigned int metFilterDecision;

  std::vector<Int_t> filters;

};

class HTTParticle
{

  

  public:

    static HTTAnalysis::sysEffects corrType;

    HTTParticle(){clear();}

    ~HTTParticle(){}

    void clear();

    ///Data member setters.
    void setP4(const TLorentzVector &aP4);
    void setChargedP4(const TLorentzVector &aP4) { chargedP4 = aP4;}
    void setNeutralP4(const TLorentzVector &aP4) { neutralP4 = aP4;}

    void setPCA(const TVector3 &aV3) {pca = aV3;}
    void setPCARefitPV(const TVector3 &aV3) {pcaRefitPV = aV3;}
    void setPCAGenPV(const TVector3 &aV3) {pcaGenPV = aV3;}

    void setCutBitmask(int bitmask) {cutBitmask = bitmask; }
    void setProperties(const std::vector<Double_t> & aProperties) { properties = aProperties;}

    ///Data member getters.
    const TLorentzVector & getP4(HTTAnalysis::sysEffects defaultType=HTTAnalysis::NOMINAL) const;
    const TVector2 & getDeltaVector() const { return deltaVector; }
    const TLorentzVector & getChargedP4() const {return chargedP4;}
    const TLorentzVector getNeutralP4() const {return neutralP4;}

    const TVector3 & getPCA() const {return pca;}
    const TVector3 & getPCARefitPV() const {return pcaRefitPV;}
    const TVector3 & getPCAGenPV() const {return pcaGenPV;}

    int getPDGid() const {return getProperty(PropertyEnum::pdgId);}
    int getCharge() const {return getProperty(PropertyEnum::charge);}
    float getMT(TVector2 met, HTTAnalysis::sysEffects defaultType=HTTAnalysis::NOMINAL) const { return TMath::Sqrt( 2. * getP4(defaultType).Pt() * met.Mod() * (1. - TMath::Cos( getP4(defaultType).Phi()-met.Phi()))); }

    int getCutBitmask() {return cutBitmask;}
    bool isBaseline()         { return (cutBitmask & 0x1) > 0; }
    bool isDiLepton()         { return (cutBitmask & 0x2) > 0; }
    bool isExtraLepton()      { return (cutBitmask & 0x4) > 0; }
    bool isAdditionalLepton() { return (cutBitmask & 0x8) > 0; }
    bool isSemiLepTau()       { return (cutBitmask & 0x10) > 0; }
    bool isFullHadLeadTau()   { return (cutBitmask & 0x20) > 0; }
    bool isFullHadSubTau()    { return (cutBitmask & 0x40) > 0; }

    Double_t getProperty(PropertyEnum index) const {return (unsigned int)index<properties.size()?  properties[(unsigned int)index]: -999;}

    bool hasTriggerMatch(TriggerEnum index) const {return (unsigned int)getProperty(PropertyEnum::isGoodTriggerType) & (1<<(unsigned int)index)
                                                           && (unsigned int)getProperty(PropertyEnum::FilterFired) & (1<<(unsigned int)index);}

   private:

    ///Return four-momentum shifted with scale.
    ///Shift modifies three-momentum transverse part only, leaving mass constant.
    const TLorentzVector getShiftedP4(HTTAnalysis::sysEffects shift) const;

    ///Nominal (as recontructed) four-momentum
    TLorentzVector p4;

    ///Scaled four-momentum;
    mutable TLorentzVector currentP4;
    mutable TVector2 deltaVector; // Shift propagated to MET
    mutable HTTAnalysis::sysEffects lastSystEffect;

    ///Charged and neutral components four-momentum
    TLorentzVector chargedP4, neutralP4;

    ///Vectors from primary vertex to point of closest approach (PCA)
    ///calculated with respect to AOD vertex, refitted and generated vertex.
    TVector3 pca, pcaRefitPV, pcaGenPV;

    ///Vector of various particle properties.
    ///Index generated automatically during conversion from
    ///LLR ntuple format
    std::vector<Double_t> properties;

    int cutBitmask = 0;
 
};
///////////////////////////////////////////////////
///////////////////////////////////////////////////
class HTTPair
{
  public:

    HTTPair(){ clear();}

    ~HTTPair(){}

    void clear();

    ///Data member setters.
    void setP4(const TLorentzVector aP4, string uncert) {p4Vector[uncert] = aP4;}

    void setLeg1(const HTTParticle &aParticle, int idx=-1){leg1 = aParticle; indexLeg1=idx;}
    void setLeg2(const HTTParticle &aParticle, int idx=-1){leg2 = aParticle; indexLeg2=idx;}

    void setMET(const TVector2 &aVector, string uncert, bool shift = true){ met[uncert] = shift ? aVector  - ( leg1.getDeltaVector() + leg2.getDeltaVector() ) : aVector ; metCache = met[uncert]; };
    void setCurrentMETShift(string uncert);
    void setMETMatrix(float m00, float m01, float m10, float m11) {metMatrix.push_back(m00); metMatrix.push_back(m01); metMatrix.push_back(m10); metMatrix.push_back(m11);}

    ///Data member getters.
    enum P4Type{ Vis = 0, Simple = 1, SVFit = 2};
    const TLorentzVector getP4(P4Type type = P4Type::SVFit);

    const TVector2 & getMET() const {return metCache;}

    float getMTLeg1(HTTAnalysis::sysEffects defaultType = HTTAnalysis::NOMINAL) const {return leg1.getMT( getMET(), defaultType ); }
    float getMTLeg2(HTTAnalysis::sysEffects defaultType = HTTAnalysis::NOMINAL) const {return leg2.getMT( getMET(), defaultType ); }
    float getMTTOT(HTTAnalysis::sysEffects defaultType = HTTAnalysis::NOMINAL) const;

    const HTTParticle & getLeg1() const {return leg1;}
    const HTTParticle & getLeg2() const {return leg2;}

    float getMVis(HTTAnalysis::sysEffects defaultType = HTTAnalysis::NOMINAL) const {return  ( leg1.getP4(defaultType) + leg2.getP4(defaultType) ).M(); }
    float getPTVis(HTTAnalysis::sysEffects defaultType = HTTAnalysis::NOMINAL) const {return ( leg1.getP4(defaultType) + leg2.getP4(defaultType) ).Pt(); }

    int getIndexLeg1() {return indexLeg1;}
    int getIndexLeg2() {return indexLeg2;}

    HTTAnalysis::finalState getFinalState();
    bool isInLooseSR();

    std::vector<float> getMETMatrix() const {return metMatrix;}

  private:

    ///Return MET modified according to given systematic effect.
    ///The MET is corrected for accorging leptons corrections.
    ///The recoil correctino is not updated.
    const TVector2 & getSystScaleMET(HTTAnalysis::sysEffects defaultType=HTTAnalysis::NOMINAL) const;

    ///Return transverse mass caluculated according to the scale shifts.
    float getSystScaleMT(const HTTParticle &aParticle, HTTAnalysis::sysEffects defaultType=HTTAnalysis::NOMINAL) const;

    ///Nominal met as calculated from PF.
    ///Includes recoil corrections.
    mutable std::map<string,TVector2> met;

    ///Scaled four-momentum cache;
    mutable TVector2 metCache;
    mutable string lastMETShift;
    mutable float mtCache;

    ///Vectors holding p4 and MET for
    ///for various scale variances.
    std::map<string, TLorentzVector> p4Vector;

    //MVAMET covariance matrix in order 00,01,10,11
    std::vector<float> metMatrix;

    ///Pair legs
    HTTParticle leg1, leg2;
    int indexLeg1, indexLeg2;

};

class HTTJet
{
  public:
    HTTJet(){clear();}
    ~HTTJet(){}

    void clear();

    void setP4(const TLorentzVector &aP4) { p4 = aP4; currentP4 = aP4; }
    void SetPtEtaPhiM(float pt, float eta, float phi, float m){ p4.SetPtEtaPhiM(pt,eta,phi,m); currentP4 = getP4(); }

    void setProperties(const std::vector<Double_t> & aProperties) { properties = aProperties;}
    void setJecUncertSourceValue(string uncert, double value, bool up){ jecUncertSourceValues[uncert] = value; };
    void setJecUncertValues( std::map<string, double> aUncertainties ){ jecUncertSourceValues = aUncertainties; }
    void setUncertShift( string uncert, bool up ){ currentP4 = getP4( uncert, up ); }

    TLorentzVector & P4(){ return currentP4; }
    float Pt(){ return currentP4.Pt(); }
    float Eta(){ return currentP4.Eta(); }
    float Phi(){ return currentP4.Phi(); }
    float M(){ return currentP4.M(); } 

    Double_t getProperty(PropertyEnum index) const {return (unsigned int)index<properties.size()?  properties[(unsigned int)index]: -999;}

    // Fallback when jec uncerts get asymmetric
    // double getJecUncertValue(unsigned int index, bool up){ return up ? jecUncertSourceValuesUp[index] : jecUncertSourceValuesDown[index]; };
    
    double getJecUncertValue(string uncert, bool up){ return jecUncertSourceValues[uncert]; };


  private:

    const TLorentzVector & getP4(string uncert = "", bool up = true);

    TLorentzVector p4;
    TLorentzVector currentP4;

    std::vector<Double_t> properties;
    std::map<string,double> jecUncertSourceValues;

    // Fallback when jec uncerts get asymmetric
    // std::vector<double> jecUncertSourceValuesUp;
    // std::vector<double> jecUncertSourceValuesDown;
};

class HTTJetCollection
{
  public:
    HTTJetCollection(){clear();};
    ~HTTJetCollection(){}

    void initCollection(bool isMC, bool applyRecoil, bool isSync);
    void clear();

    void initForPromoteDemote();

    void addJet(HTTJet newJet){ jetCollection.push_back(newJet); };
    void setCurrentUncertShift( string uncert, bool up ){ fillCurrentCollections(uncert,up); }
    void fillCurrentCollections(string uncert = "", bool up = true);

    vector< pair< string, pair<string,bool> > > getNeededJECShifts(){ return jecShifts; }
    const TVector2 getTotalJetShift(string uncert, bool up);
    HTTJet & getJet(unsigned int index){ return jetCurrentCollection[index]; }

    void btagPromoteDemote(string mistagsys = "central", string btagsys = "central");
    HTTJet & getBtagJet(unsigned int index){ return btagCurrentCollection[index]; }

    int getNJets(double pt = 20);
    int getNBtag(){ return btagCurrentCollection.size(); };
    int getNJetInGap(double pt = 20);


    const TLorentzVector & getDiJetP4(){return dijet;}

  private:

    static bool sortJets(HTTJet j1, HTTJet j2);

    bool usePromoteDemote = false;
    BTagCalibration calib;
    BTagCalibrationReader reader;
    TFile *eff_file;
    TH2D *hb_eff;
    TH2D *hc_eff;
    TH2D *hoth_eff;

    // Names of JES uncertainty shifts.
    vector< pair< string, pair<string,bool> > > jecShifts;

    std::vector<HTTJet> jetCollection;
    std::vector<HTTJet> jetCurrentCollection;
    std::vector<HTTJet> btagCollection;
    std::vector<HTTJet> btagCurrentCollection;
    std::vector<HTTJet> antibtagCollection;
    TLorentzVector dijet;
    void setDijetP4();
    void fillCurrentCollection();
    

};

#endif
