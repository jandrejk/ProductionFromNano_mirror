#define HTauhTauhTreeFromNano_cxx
#include "HTauhTauhTreeFromNano.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauhTauhTreeFromNano::pairSelection(unsigned int iPair){

  ///Baseline+post sync selection as on
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#Baseline_tau_h_tau_h
  ///Indices for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
  ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

  if(httPairs_.empty()) return false;

  int pdgIdLeg1 = httPairs_[iPair].getLeg1().getPDGid();
  int pdgIdLeg2 = httPairs_[iPair].getLeg2().getPDGid();

  if( std::abs(pdgIdLeg1)!=15 || std::abs(pdgIdLeg2)!=15 ) return 0;

  unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
  unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
  //MB sort taus within the pair
  double pt_1 = httLeptonCollection[indexLeg1].getP4().Pt();
  double pt_2 = httLeptonCollection[indexLeg2].getP4().Pt();

  if(pt_2>pt_1){//tau with higher-Pt first
    unsigned int indexLegTmp = indexLeg1;
    indexLeg1 = indexLeg2;
    indexLeg2 = indexLegTmp;
  }

  TLorentzVector tau1P4 = httLeptonCollection[indexLeg1].getP4();
  TLorentzVector tau2P4 = httLeptonCollection[indexLeg2].getP4();

  bool tauBaselineSelection1 =  httLeptonCollection[indexLeg1].isFullHadLeadTau();
  bool tauBaselineSelection2 = httLeptonCollection[indexLeg2].isFullHadSubTau();

  bool baselinePair = tau1P4.DeltaR(tau2P4) > 0.5;

  bool boolAntiEleLeg1 = ( (int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idAntiEle) & 0x1) == 0x1;   //Vloose AntiEle Id
  bool boolAntiEleLeg2 = ( (int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idAntiEle) & 0x1) == 0x1;   //Vloose AntiEle Id
  bool boolAntiMuLeg1  = ( (int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idAntiMu)  & 0x1) == 0x1;   //Loose AntiMu Id
  bool boolAntiMuLeg2  = ( (int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idAntiMu)  & 0x1) == 0x1;   //Loose AntiMu Id
  bool boolAntiLeptonId = boolAntiEleLeg1 && boolAntiMuLeg1 && boolAntiEleLeg2 && boolAntiMuLeg2;


  bool boolExtraElectronVeto = thirdLeptonVeto(indexLeg1, indexLeg2, 11);
  bool boolExtraMuonVeto     = thirdLeptonVeto(indexLeg1, indexLeg2, 13);
  bool boolThirdLeptonVeto   = boolExtraMuonVeto || boolExtraElectronVeto;

  httEvent->clearSelectionWord();
  httEvent->setSelectionBit(SelectionBitsEnum::diMuonVeto,0); //only set explicitly for mutau
  httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,0); //only set explicitly for etau  
  httEvent->setSelectionBit(SelectionBitsEnum::diLeptonVeto, 0);
  httEvent->setSelectionBit(SelectionBitsEnum::antiLeptonId, boolAntiLeptonId);
  httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,boolExtraMuonVeto);
  httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto, boolExtraElectronVeto);
  httEvent->setSelectionBit(SelectionBitsEnum::thirdLeptonVeto, boolThirdLeptonVeto);

  return tauBaselineSelection1 && tauBaselineSelection2 && baselinePair
         && true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauhTauhTreeFromNano::bestPair(std::vector<unsigned int> &pairIndexes){

  unsigned int bestIndex = 9999;
  ///Pairs are already sorted during the ntuple creation?
  double iso_1=std::numeric_limits<double>::infinity(), iso_2=std::numeric_limits<double>::infinity(), pt_1=-1, pt_2=-1;
  if(pairIndexes.size()) {
    //return pairIndexes[0];//MB
    for(unsigned int ii=0;ii<2*pairIndexes.size();++ii){
      unsigned int i=(ii<pairIndexes.size()?ii:ii-pairIndexes.size());
      unsigned int iPair = pairIndexes[i];
      unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
      unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
      if(ii>=pairIndexes.size()){//invert legs
        indexLeg1 = httPairs_[iPair].getIndexLeg2();
        indexLeg2 = httPairs_[iPair].getIndexLeg1();
      }
      double pt_1_i = httLeptonCollection[indexLeg1].getP4().Pt();
      double pt_2_i = httLeptonCollection[indexLeg2].getP4().Pt();
      //MB: More isolated for MVAIso means higher value so inverted here to keep standard convention in comparison
      double iso_1_i = -httLeptonCollection[indexLeg1].getProperty(PropertyEnum::rawMVAoldDM2017v2);
      double iso_2_i = -httLeptonCollection[indexLeg2].getProperty(PropertyEnum::rawMVAoldDM2017v2);

      if(iso_1_i>iso_1) continue;
      if(iso_1_i==iso_1 && pt_1_i<pt_1) continue;
      if(iso_2_i>iso_2) continue;
      if(iso_2_i==iso_2 && pt_2_i<pt_2) continue;
      bestIndex = iPair;
      iso_1 = iso_1_i;
      iso_2 = iso_2_i;
      pt_1 = pt_1_i;
      pt_2 = pt_2_i;
    }
  }
  /*
  if(pairIndexes.size() && bestIndex!=pairIndexes[0]){
    unsigned int iPair = pairIndexes[0];
    unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
    unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
    double pt_1_i = httLeptonCollection[indexLeg1].getP4().Pt();
    double pt_2_i = httLeptonCollection[indexLeg2].getP4().Pt();
    std::cout<<"Pair sorting: "<<std::endl
             <<"best index = "<<bestIndex<<", index[0] = "<<pairIndexes[0]<<std::endl
             <<"\tiso1[best]="<<-iso_1<<", iso1[0]="<<httLeptonCollection[indexLeg1].getProperty(PropertyEnum::rawMVAoldDM)<<std::endl
             <<"\tpt1[best]="<<pt_1<<", pt1[0]="<<pt_1_i<<std::endl
             <<"\tiso2[best]="<<-iso_2<<", iso1[0]="<<httLeptonCollection[indexLeg2].getProperty(PropertyEnum::rawMVAoldDM)<<std::endl
             <<"\tpt2[best]="<<pt_2<<", pt1[0]="<<pt_2_i<<std::endl;
  }
  */

  return bestIndex;
};
/////////////////////////////////////////////////
/////////////////////////////////////////////////
