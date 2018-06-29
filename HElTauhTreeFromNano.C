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

  //##################################################################
  if (event==check_event_number)  cout << "pS1 " << httPairs_.empty() << endl;
  //##################################################################

  if(httPairs_.empty()) return false;

  //##################################################################
  if (event==check_event_number)
  {
      cout << "pS2 ";
      cout << std::abs(httPairs_[iPair].getLeg1().getPDGid()) << "  ";
      cout << std::abs(httPairs_[iPair].getLeg2().getPDGid()) << endl;
  }
  //##################################################################
  
  int pdgIdLeg1 = std::abs(httPairs_[iPair].getLeg1().getPDGid());
  int pdgIdLeg2 = std::abs(httPairs_[iPair].getLeg2().getPDGid());

  unsigned int indexElecLeg = -1;
  if(pdgIdLeg1==11)      indexElecLeg = httPairs_[iPair].getIndexLeg1();
  else if(pdgIdLeg2==11) indexElecLeg = httPairs_[iPair].getIndexLeg2();
  else return 0;

  if (event==check_event_number) cout << "pS3 " << endl;

  unsigned int indexTauLeg = -1;
  if(pdgIdLeg1==15)      indexTauLeg = httPairs_[iPair].getIndexLeg1();
  else if(pdgIdLeg2==15) indexTauLeg = httPairs_[iPair].getIndexLeg2();
  else return 0;

  if (event==check_event_number) cout << "pS4 " << endl;

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

  //  if (event==check_event_number) cout << "pS4 " << endl;

    // bool electronBaselineSelection =  electronSelection(indexElecLeg);
    bool electronBaselineSelection = httLeptonCollection[indexElecLeg].isBaseline();

    // bool tauBaselineSelection = tauSelection(indexTauLeg);
    bool tauBaselineSelection = tauP4.Pt()> LeptonCuts::Baseline.Tau.SemiLep.pt
                                 && std::abs(tauP4.Eta())<LeptonCuts::Baseline.Tau.SemiLep.eta;



    if (event==check_event_number){

        cout << "pS4a " << elecP4.Pt() << " ";
        cout << elecP4.Eta() << " ";
        cout << std::abs(httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::dz)) << " ";
        cout << std::abs(httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::dxy))<< " ";
        cout << (int)std::abs(httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::mvaFall17Iso_WP80))<< endl;

        cout << "pS4b " << tauP4.Pt() << " ";
        cout << tauP4.Eta() << " ";
        cout << httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idDecayMode) << " ";
        cout << std::abs(httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::dz)) << " ";
        cout << (int)std::abs(httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::charge))<<  endl;
    }


    bool baselinePair = elecP4.DeltaR(tauP4) > 0.5;
    bool postSynchElectron = httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::pfRelIso03_all)<0.1;
    bool loosePostSynchElectron = httLeptonCollection[indexElecLeg].getProperty(PropertyEnum::pfRelIso03_all)<0.3;
    bool postSynchTau = (tauID & tauIDmask) == tauIDmask;

      ///SUSY synch selection
      //electronBaselineSelection &= elecP4.Pt()>23 && std::abs(elecP4.Eta())<2.1;
      //tauBaselineSelection &= tauP4.Pt()>30 && std::abs(tauP4.Eta())<2.3;
      ///////////////////////

    httEvent->clearSelectionWord();
    httEvent->setSelectionBit(SelectionBitsEnum::electronBaselineSelection,electronBaselineSelection);
    httEvent->setSelectionBit(SelectionBitsEnum::tauBaselineSelection,tauBaselineSelection);
    httEvent->setSelectionBit(SelectionBitsEnum::baselinePair,baselinePair);
    httEvent->setSelectionBit(SelectionBitsEnum::postSynchElectron,postSynchElectron);
    httEvent->setSelectionBit(SelectionBitsEnum::postSynchTau,postSynchTau);
    httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,diElectronVeto());
    httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,thirdLeptonVeto(indexElecLeg, indexTauLeg, 13));
    httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto,thirdLeptonVeto(indexElecLeg, indexTauLeg, 11));

    if (event==check_event_number) cout << "pS5 " << electronBaselineSelection << " " << tauBaselineSelection << " " << baselinePair << endl;

  return electronBaselineSelection && tauBaselineSelection && baselinePair
    //&& postSynchTau && loosePostSynchMuon //comment out for sync
    //&& !diMuonVeto() && !thirdLeptonVeto(indexElecLeg, indexTauLeg, 13) && !thirdLeptonVeto(indexElecLeg, indexTauLeg, 11) //comment out for sync
    && true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HElTauhTreeFromNano::electronSelection(unsigned int index){

  TLorentzVector aP4 = httLeptonCollection[index].getP4();


    return aP4.Pt()>LeptonCuts::Baseline.Electron.pt 
           && std::abs(aP4.Eta())<=LeptonCuts::Baseline.Electron.eta
           && std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dz)) <0.2
           && std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dxy))<0.045
           && httLeptonCollection[index].getProperty(PropertyEnum::convVeto)>0.5
           && httLeptonCollection[index].getProperty(PropertyEnum::lostHits)<1.5 //0 or 1
           && httLeptonCollection[index].getProperty(HTTEvent::usePropertyFor.at("electronIDWP80") )>0.5;        

}

bool HElTauhTreeFromNano::tauSelection(unsigned int index){

    TLorentzVector aP4 = httLeptonCollection[index].getP4();

    return  aP4.Pt()> 20
            && std::abs(aP4.Eta())<2.3 
            && httLeptonCollection[index].getProperty(PropertyEnum::idDecayMode)>0.5 
            && std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dz))<0.2 
            && (int)std::abs(httLeptonCollection[index].getProperty(PropertyEnum::charge))==1;

}



bool HElTauhTreeFromNano::diElectronVeto(){

    std::vector<int> elecIndexes;
    for(unsigned int iLepton=0;iLepton<httLeptonCollection.size();++iLepton){

      if(std::abs(httLeptonCollection[iLepton].getPDGid())!=11) continue;
      TLorentzVector elecP4 = httLeptonCollection[iLepton].getP4();

        // bool passLepton = elecP4.Pt()> LeptonCuts::Di.Electron.pt
        //      && std::abs(elecP4.Eta())< LeptonCuts::Di.Electron.eta
        //      && std::abs(httLeptonCollection[iLepton].getProperty(PropertyEnum::dz))<0.2
        //      && std::abs(httLeptonCollection[iLepton].getProperty(PropertyEnum::dxy))<0.045
        //      && httLeptonCollection[iLepton].getProperty( HTTEvent::usePropertyFor["electronIsolation"] )<0.3
        //      && true;
        bool passLepton = httLeptonCollection[iLepton].isDiLepton();
        if(passLepton) elecIndexes.push_back(iLepton);
    }

    if(elecIndexes.size()<2) return false;
    else{
        for(unsigned int iElec1=0;iElec1<elecIndexes.size()-1;++iElec1){
            for(unsigned int iElec2=iElec1+1;iElec2<elecIndexes.size();++iElec2){

                TLorentzVector elec1P4 = httLeptonCollection[iElec1].getP4();
                TLorentzVector elec2P4 = httLeptonCollection[iElec2].getP4();
                int elec1Charge = (int)httLeptonCollection[iElec1].getProperty(PropertyEnum::charge);
                int elec2Charge = (int)httLeptonCollection[iElec2].getProperty(PropertyEnum::charge);

                if(elec2Charge*elec1Charge==-1 && elec1P4.DeltaR(elec2P4)>0.15) return true;
            }
        }
    }
    return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
