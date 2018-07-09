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

    int tauIDmask = 0;
    int tauID = (int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiMu);
    tauID += (int)std::pow(2,HTTEvent::againstEIdOffset)*(int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idAntiEle);
    tauID += (int)std::pow(2,HTTEvent::mvaIsoIdOffset)*(int)httLeptonCollection[indexTauLeg].getProperty(PropertyEnum::idMVAoldDM2017v2);

    for(unsigned int iBit=0;iBit<HTTEvent::ntauIds;iBit++)
    {
        if(HTTEvent::tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT") tauIDmask |= (1<<iBit);
        if(HTTEvent::tauIDStrings[iBit]=="againstMuonTight3") tauIDmask |= (1<<iBit);
        if(HTTEvent::tauIDStrings[iBit]=="againstElectronVLooseMVA6") tauIDmask |= (1<<iBit);
    }

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
    bool postSynchMuon = httLeptonCollection[indexMuonLeg].getProperty(PropertyEnum::pfRelIso04_all)<0.15;
    bool loosePostSynchMuon = httLeptonCollection[indexMuonLeg].getProperty(PropertyEnum::pfRelIso04_all)<0.3;
    bool postSynchTau = (tauID & tauIDmask) == tauIDmask;

    debugWayPoint("[pairSelection] Overlap",{(double)muonP4.DeltaR(tauP4)}, {(int)baselinePair},{"DR","noOL"});

    httEvent->clearSelectionWord();
    httEvent->setSelectionBit(SelectionBitsEnum::muonBaselineSelection,muonBaselineSelection);
    httEvent->setSelectionBit(SelectionBitsEnum::tauBaselineSelection,tauBaselineSelection);
    httEvent->setSelectionBit(SelectionBitsEnum::baselinePair,baselinePair);
    httEvent->setSelectionBit(SelectionBitsEnum::postSynchMuon,postSynchMuon);
    httEvent->setSelectionBit(SelectionBitsEnum::postSynchTau,postSynchTau);
    httEvent->setSelectionBit(SelectionBitsEnum::diMuonVeto,diMuonVeto());
    httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,thirdLeptonVeto(indexMuonLeg, indexTauLeg, 13));
    httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto,thirdLeptonVeto(indexMuonLeg, indexTauLeg, 11));


    return muonBaselineSelection && tauBaselineSelection && baselinePair
          //&& postSynchTau && loosePostSynchMuon //comment out for sync
          //&& !diMuonVeto() && !thirdLeptonVeto(indexMuonLeg, indexTauLeg, 13) && !thirdLeptonVeto(indexMuonLeg, indexTauLeg, 11) //comment out for sync
          && true;
}

bool HMuTauhTreeFromNano::diMuonVeto(){

    std::vector<int> muonIndexes;
    for(unsigned int iLepton=0;iLepton<httLeptonCollection.size();++iLepton)
    {

        if(std::abs(httLeptonCollection[iLepton].getPDGid())!=13) continue;
        TLorentzVector muonP4 = httLeptonCollection[iLepton].getP4();

        // bool passLepton = muonP4.Pt()> LeptonCuts::Di.Muon.pt
        //                   && std::abs(muonP4.Eta())< LeptonCuts::Di.Muon.pt
        //                   && std::abs(httLeptonCollection[iLepton].getProperty(PropertyEnum::dz))<0.2
        //                   && std::abs(httLeptonCollection[iLepton].getProperty(PropertyEnum::dxy))<0.045
        //                   && httLeptonCollection[iLepton].getProperty( HTTEvent::usePropertyFor["muonIsolation"] )<0.3
        //                   //FIXME && ((daughters_typeOfMuon->at(iLepton) & ((1<<0) + (1<<1) + (1<<2))) == ((1<<0) + (1<<1) + (1<<2))) //0=PF, 1=Global, 2=Tracker, muons in Nano are loose, i.e. PF&(Global|Tracker)
        //                   && true;

        bool passLepton = httLeptonCollection[iLepton].isDiLepton();

        if(passLepton) muonIndexes.push_back(iLepton);
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
