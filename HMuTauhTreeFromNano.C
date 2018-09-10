#define HMuTauhTreeFromNano_cxx
#include "HMuTauhTreeFromNano.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HMuTauhTreeFromNano::pairSelection(unsigned int iPair)
{

    ///Baseline+post sync selection as on
    ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#Baseline_mu_tau_h
    ///Indices for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
    ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

    debugWayPoint("[pairSelection] ------ Begin -------");

    if(httPairs_.empty()) return false;

    int pdgIdLeg1 = std::abs(httPairs_[iPair].getLeg1().getPDGid());
    int pdgIdLeg2 = std::abs(httPairs_[iPair].getLeg2().getPDGid());

    debugWayPoint("[pairSelection] pdgid of pair",{},{(int)pdgIdLeg1,(int)pdgIdLeg2});

    unsigned int indexMuonLeg = -1;
    if(pdgIdLeg1==13)      indexMuonLeg = httPairs_[iPair].getIndexLeg1();
    else if(pdgIdLeg2==13) indexMuonLeg = httPairs_[iPair].getIndexLeg2();
    else return 0;

    debugWayPoint("[pairSelection] Found Muon",{},{(int)indexMuonLeg});

    unsigned int indexTauLeg = -1;
    if(pdgIdLeg1==15)      indexTauLeg = httPairs_[iPair].getIndexLeg1();
    else if(pdgIdLeg2==15) indexTauLeg = httPairs_[iPair].getIndexLeg2();
    else return 0;

    debugWayPoint("[pairSelection] Found Tau",{},{(int)indexTauLeg});

    TLorentzVector muonP4 = httLeptonCollection[indexMuonLeg].getP4();
    TLorentzVector tauP4 = httLeptonCollection[indexTauLeg].getP4();

    bool muonBaselineSelection = httLeptonCollection[indexMuonLeg].isBaseline();

    debugWayPoint("[pairSelection] Muon",
                  {(double)muonP4.Pt(), (double)muonP4.Eta()},
                  {(int)muonBaselineSelection},
                  {"pt","eta","passCuts"});

    bool tauBaselineSelection = httLeptonCollection[indexTauLeg].isSemiLepTau();

    debugWayPoint("[pairSelection] Tau",
              {(double)tauP4.Pt(), (double)tauP4.Eta()},
              {(int)tauBaselineSelection},
              {"pt","eta","passCuts"});


    bool baselinePair = muonP4.DeltaR(tauP4) > 0.5;

    debugWayPoint("[pairSelection] Overlap",{(double)muonP4.DeltaR(tauP4)}, {(int)baselinePair},{"DR","noOL"});

    bool boolDiElectronVeto = false;
    bool boolDiMuonVeto = diMuonVeto();
    bool boolDiLeptonVeto = boolDiElectronVeto || boolDiMuonVeto;

    bool boolAntiEle = ( (int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiEle) & 0x1) == 0x1;   //Vloose AntiEle Id
    bool boolAntiMu  = ( (int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiMu)  & 0x2) == 0x2;   //Tight AntiMu Id
    bool boolAntiLeptonId = boolAntiEle && boolAntiMu;

    bool boolExtraElectronVeto = thirdLeptonVeto(indexMuonLeg, indexTauLeg, 11);
    bool boolExtraMuonVeto     = thirdLeptonVeto(indexMuonLeg, indexTauLeg, 13);
    bool boolThirdLeptonVeto   = boolExtraMuonVeto || boolExtraElectronVeto;

    httEvent->clearSelectionWord();
    httEvent->setSelectionBit(SelectionBitsEnum::diMuonVeto,boolDiMuonVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,    boolDiElectronVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::diLeptonVeto,      boolDiLeptonVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::antiLeptonId,      boolAntiLeptonId);
    httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,     boolExtraMuonVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto, boolExtraElectronVeto);
    httEvent->setSelectionBit(SelectionBitsEnum::thirdLeptonVeto,   boolThirdLeptonVeto);

    

    return muonBaselineSelection && tauBaselineSelection && baselinePair
           && true;
}

bool HMuTauhTreeFromNano::diMuonVeto(){

    std::vector<int> muonIndexes;
    for(unsigned int iLepton=0;iLepton<httLeptonCollection.size();++iLepton)
    {
        if(std::abs(httLeptonCollection[iLepton].getPDGid())!=13) continue;
        if(httLeptonCollection[iLepton].isDiLepton() ) muonIndexes.push_back(iLepton);
    }

    if(muonIndexes.size()<2) return false;
    else{
        for(auto &iMuon1 : muonIndexes)
        {
            for(auto &iMuon2 : muonIndexes)
            {
                if(iMuon1 == iMuon2) continue;

                TLorentzVector muon1P4 = httLeptonCollection[iMuon1].getP4();
                TLorentzVector muon2P4 = httLeptonCollection[iMuon2].getP4();
                int muon1Charge = (int)httLeptonCollection[iMuon1].getProperty(PropertyEnum::charge);
                int muon2Charge = (int)httLeptonCollection[iMuon2].getProperty(PropertyEnum::charge);

                float deltaR = muon1P4.DeltaR(muon2P4);
                if(muon2Charge*muon1Charge==-1 && deltaR>0.15) return true;
            }
        }
    }
    return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
