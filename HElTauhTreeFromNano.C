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

  int tauIDmask = 0;
  int tauID = (int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiMu);
  tauID += (int)std::pow(2,HTTEvent::againstEIdOffset)*(int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiEle);
  tauID += (int)std::pow(2,HTTEvent::mvaIsoIdOffset)*(int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idMVAoldDM2017v2);

  for(unsigned int iBit=0;iBit<HTTEvent::ntauIds;iBit++){
    if(HTTEvent::tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
    if(HTTEvent::tauIDStrings[iBit]=="againstMuonLoose3") tauIDmask |= (1<<iBit);
    if(HTTEvent::tauIDStrings[iBit]=="againstElectronTightMVA6") tauIDmask |= (1<<iBit);
  }

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
    bool postSynchElectron = httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::pfRelIso03_all)<0.1;
    bool loosePostSynchElectron = httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::pfRelIso03_all)<0.3;
    bool postSynchTau = (tauID & tauIDmask) == tauIDmask;

    debugWayPoint("[pairSelection] Overlap",{(double)elecP4.DeltaR(tauP4)}, {(int)baselinePair},{"DR","noOL"});

    httEvent->clearSelectionWord();
    httEvent->setSelectionBit(SelectionBitsEnum::electronBaselineSelection,electronBaselineSelection);
    httEvent->setSelectionBit(SelectionBitsEnum::tauBaselineSelection,tauBaselineSelection);
    httEvent->setSelectionBit(SelectionBitsEnum::baselinePair,baselinePair);
    httEvent->setSelectionBit(SelectionBitsEnum::postSynchElectron,postSynchElectron);
    httEvent->setSelectionBit(SelectionBitsEnum::postSynchTau,postSynchTau);
    httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,diElectronVeto());
    httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,thirdLeptonVeto(indexElecLeg, indexTauLeg, 13));
    httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto,thirdLeptonVeto(indexElecLeg, indexTauLeg, 11));



  return electronBaselineSelection && tauBaselineSelection && baselinePair
    //&& postSynchTau && loosePostSynchMuon //comment out for sync
    //&& !diMuonVeto() && !thirdLeptonVeto(indexElecLeg, indexTauLeg, 13) && !thirdLeptonVeto(indexElecLeg, indexTauLeg, 11) //comment out for sync
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
