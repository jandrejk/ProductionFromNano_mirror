#define HElTauhTreeFromNano_cxx
#include "HElTauhTreeFromNano.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HElTauhTreeFromNano::pairSelection(unsigned int iPair){

  ///Baseline+post sync selection as on
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#Baseline_e_tau_h
  ///Indices for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
  ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc
  debugWayPoint("[pairSelection] ------ Begin -------");
  if(httPairs_.empty()) return false;

  //##################################################################
  
  int pdgIdLeg1 = std::abs(httPairs_[iPair].getLeg1().getPDGid());
  int pdgIdLeg2 = std::abs(httPairs_[iPair].getLeg2().getPDGid());

  debugWayPoint("[pairSelection] pdgid of pair",{},{(int)pdgIdLeg1,(int)pdgIdLeg2});

  unsigned int indexElecLeg = -1;
  if(pdgIdLeg1==11)      indexElecLeg = httPairs_[iPair].getIndexLeg1();
  else if(pdgIdLeg2==11) indexElecLeg = httPairs_[iPair].getIndexLeg2();
  else return 0;

  debugWayPoint("[pairSelection] Found Electron",{},{(int)indexElecLeg});

  unsigned int indexTauLeg = -1;
  if(pdgIdLeg1==15)      indexTauLeg = httPairs_[iPair].getIndexLeg1();
  else if(pdgIdLeg2==15) indexTauLeg = httPairs_[iPair].getIndexLeg2();
  else return 0;

  debugWayPoint("[pairSelection] Found Tau",{},{(int)indexTauLeg});

  TLorentzVector elecP4 = httLeptonCollection[indexElecLeg].getP4();
  TLorentzVector tauP4 =  httLeptonCollection[indexTauLeg].getP4();


    bool electronBaselineSelection = httLeptonCollection[indexElecLeg].isBaseline();

    debugWayPoint("[pairSelection] Electron",
                  {(double)elecP4.Pt(), (double)elecP4.Eta()},
                  {(int)electronBaselineSelection},
                  {"pt","eta","passCuts"});

    bool tauBaselineSelection = httLeptonCollection[indexTauLeg].isSemiLepTau();

    debugWayPoint("[pairSelection] Tau",
              {(double)tauP4.Pt(), (double)tauP4.Eta()},
              {(int)tauBaselineSelection},
              {"pt","eta","passCuts"});

    bool baselinePair = elecP4.DeltaR(tauP4) > 0.5;

    debugWayPoint("[pairSelection] Overlap",{(double)elecP4.DeltaR(tauP4)}, {(int)baselinePair},{"DR","noOL"});

    bool boolDiElectronVeto = diElectronVeto();
    bool boolDiMuonVeto = false;
    bool boolDiLeptonVeto = boolDiElectronVeto || boolDiMuonVeto;

    bool boolAntiEle = ( (int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiEle) & 0x8) == 0x8;   //Tight AntiEle Id
    bool boolAntiMu  = ( (int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiMu)  & 0x1) == 0x1;   //Loose AntiMu Id
    bool boolAntiLeptonId = boolAntiEle && boolAntiMu;

    bool boolExtraElectronVeto = thirdLeptonVeto(indexElecLeg, indexTauLeg, 11);
    bool boolExtraMuonVeto     = thirdLeptonVeto(indexElecLeg, indexTauLeg, 13);
    bool boolThirdLeptonVeto   = boolExtraMuonVeto || boolExtraElectronVeto;

    httEvent->clearSelectionWord();
    httEvent->setSelectionBit(SelectionBitsEnum::diMuonVeto,boolDiMuonVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,boolDiElectronVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::diLeptonVeto, boolDiLeptonVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,boolExtraMuonVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto, boolExtraElectronVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::thirdLeptonVeto, boolThirdLeptonVeto);



  return electronBaselineSelection && tauBaselineSelection && baselinePair
         && true;
}

bool HElTauhTreeFromNano::diElectronVeto(){

    std::vector<int> elecIndexes;
    for(unsigned int iLepton=0;iLepton<httLeptonCollection.size();++iLepton){

      if(std::abs(httLeptonCollection[iLepton].getPDGid())!=11) continue;
      TLorentzVector elecP4 = httLeptonCollection[iLepton].getP4();

        bool passLepton = httLeptonCollection[iLepton].isDiLepton();
        debugWayPoint("[diElectronVeto] dilep cut",{},
                      {(int)iLepton,(int)passLepton},
                      {"iLep","pass"});

        if(passLepton) elecIndexes.push_back(iLepton);
    }

    if(elecIndexes.size()<2) return false;
    else{
        for(auto &iE1 : elecIndexes)
        {
          for(auto &iE2 : elecIndexes)
          {
                if(iE1 == iE2) continue;

                TLorentzVector elec1P4 = httLeptonCollection[iE1].getP4();
                TLorentzVector elec2P4 = httLeptonCollection[iE2].getP4();
                int elec1Charge = (int)httLeptonCollection[iE1].getProperty(PropertyEnum::charge);
                int elec2Charge = (int)httLeptonCollection[iE2].getProperty(PropertyEnum::charge);

                debugWayPoint("[diElectronVeto] find dilep",
                              {(double)elec1P4.DeltaR(elec2P4) },
                              {(int)elec1Charge,(int)elec2Charge,(int)iE1,(int)iE2},
                              {"DR","q1","q2","iE1","iE2"});

                if(elec2Charge*elec1Charge==-1 && elec1P4.DeltaR(elec2P4)>0.15) return true;
            }
        }
    }
    return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
