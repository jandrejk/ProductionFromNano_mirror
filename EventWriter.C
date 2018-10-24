#include "EventWriter.h"

const float DEF = -10.;
const float DEFWEIGHT = 1.0;

const int gen_el_map[24]={ 6, 1,6,6,6,6, 6,6,6,6,6, 6,6,6,6,3, 6,6,6,6,6, 6,6,6 }; 
const int gen_mu_map[24]={ 6, 2,6,6,6,6, 6,6,6,6,6, 6,6,6,6,4, 6,6,6,6,6, 6,6,6 }; 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fill(HTTEvent *ev, HTTJetCollection *jets, std::vector<HTTParticle> leptons, HTTPair *pair){

    // Make sure that current Collections are filled
    jets->setCurrentUncertShift( "", true );
    pair->setCurrentMETShift("");

    channel = pair->getFinalState();

    leg1 = pair->getLeg1();
    leg2 = pair->getLeg2();

    leg1P4=leg1.getP4();
    leg2P4=leg2.getP4();

    gen_match_1=leg1.getProperty(PropertyEnum::mc_match);
    gen_match_2=leg2.getProperty(PropertyEnum::mc_match);

    pdg1=std::abs(leg1.getProperty(PropertyEnum::pdgId));
    pdg2=std::abs(leg2.getProperty(PropertyEnum::pdgId));

    run_syncro=ev->getRunId();
    lumi_syncro=ev->getLSId();
    evt_syncro=ev->getEventId();
    //  entry=DEF; //filled in main loop
    //  fileEntry=DEF; //filled in main loop

    unsigned genFlav1=leg1.getProperty(PropertyEnum::genPartFlav);
    unsigned genFlav2=leg2.getProperty(PropertyEnum::genPartFlav);

    npv=ev->getNPV();
    npvGood=DEF;
    npu=ev->getNPU();
    rho=ev->getRho();

    fillLeg1Branches();
    fillLeg2Branches();
    fillJetBranches(jets);
    fillPairBranches(pair, jets);
    fillAdditionalLeptons( leptons, pair );

    gen_top_pt_1=DEF;
    gen_top_pt_2=DEF;
    genJets=DEF;
    uncorrmet=ev->getMET_uncorr().Mod();
    //////////////////////////////////////////////////////////////////  

    failBadGlobalMuonTagger=DEF;
    failCloneGlobalMuonTagger=DEF;
    Flag_duplicateMuons=DEF;
    Flag_noBadMuons=DEF;


    Flag_goodVertices =                       ev->getFilter(FilterEnum::Flag_goodVertices);
    Flag_globalTightHalo2016Filter =          ev->getFilter(FilterEnum::Flag_globalTightHalo2016Filter);
    Flag_HBHENoiseFilter =                    ev->getFilter(FilterEnum::Flag_HBHENoiseFilter);
    Flag_HBHENoiseIsoFilter =                 ev->getFilter(FilterEnum::Flag_HBHENoiseIsoFilter);
    Flag_EcalDeadCellTriggerPrimitiveFilter = ev->getFilter(FilterEnum::Flag_EcalDeadCellTriggerPrimitiveFilter);
    Flag_BadPFMuonFilter =                    ev->getFilter(FilterEnum::Flag_BadPFMuonFilter);
    Flag_BadChargedCandidateFilter =          ev->getFilter(FilterEnum::Flag_BadChargedCandidateFilter);
    Flag_eeBadScFilter =                      ev->getFilter(FilterEnum::Flag_eeBadScFilter);
    Flag_ecalBadCalibFilter =                 ev->getFilter(FilterEnum::Flag_ecalBadCalibFilter);
    Flag_METFilters =                         ev->getFilter(FilterEnum::Flag_METFilters);    

    flagMETFilter = Flag_goodVertices 
                    && Flag_globalTightHalo2016Filter 
                    && Flag_HBHENoiseFilter
                    && Flag_HBHENoiseIsoFilter
                    && Flag_EcalDeadCellTriggerPrimitiveFilter
                    && Flag_BadPFMuonFilter
                    && Flag_BadChargedCandidateFilter
                    && Flag_eeBadScFilter
                    && Flag_ecalBadCalibFilter;



    //////////////////////////////////////////////////////////////////  

    extramuon_veto=ev->checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
    extraelec_veto=ev->checkSelectionBit(SelectionBitsEnum::extraElectronVeto);
    diMuonVeto= ev->checkSelectionBit(SelectionBitsEnum::diMuonVeto);
    diElectronVeto= ev->checkSelectionBit(SelectionBitsEnum::diElectronVeto);
    dilepton_veto = ev->checkSelectionBit(SelectionBitsEnum::diLeptonVeto);

    passesDiMuonVeto=!diMuonVeto;           // silly to pass a veto...
    passesDiElectronVeto=!diElectronVeto;   // kept for backward compatibility

    passesThirdLepVeto=!ev->checkSelectionBit(SelectionBitsEnum::thirdLeptonVeto); 
    passesTauLepVetos = ev->checkSelectionBit(SelectionBitsEnum::antiLeptonId);    


    //////////////////////////////////////////////////////////////////
    if (channel == HTTAnalysis::TauTau) passesIsoCuts=byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 
                                                      && byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2;//tautau
    else passesIsoCuts=byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2; //etau/mutau

    if (pdg1==11 || pdg1==13) passesLepIsoCuts=(iso_1<0.1);


    //////////////////////////////////////////////////////////////////

    //this is quite slow, calling the function for each trigger item...
    if ( channel == HTTAnalysis::MuTau ) //mu-tau
    {
        trg_singletau_leading= false;
        trg_singletau_trailing= leg2.hasTriggerMatch(TriggerEnum::HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);

        trg_singlemuon_27 =    leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu27) && pt_1 > 28;
        trg_singlemuon_24 =    leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu24) && pt_1 > 25;

        trg_crossmuon_mu20tau27=  leg1.hasTriggerMatch(TriggerEnum::HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1) 
                                  && leg2.hasTriggerMatch(TriggerEnum::HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1)
                                  && pt_1 > 21 && pt_2 > 32;

    }else if ( channel == HTTAnalysis::EleTau )
    {
        trg_singletau_leading= false;
        trg_singletau_trailing= leg2.hasTriggerMatch(TriggerEnum::HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);

        trg_singleelectron_35 =    leg1.hasTriggerMatch(TriggerEnum::HLT_Ele35_WPTight_Gsf)  && pt_1 > 36;
        trg_singleelectron_32 =    leg1.hasTriggerMatch(TriggerEnum::HLT_Ele32_WPTight_Gsf)  && pt_1 > 33;
        trg_singleelectron_27 =    leg1.hasTriggerMatch(TriggerEnum::HLT_Ele27_WPTight_Gsf)  && pt_1 > 28;

        trg_crossele_ele24tau30 = leg1.hasTriggerMatch(TriggerEnum::HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1) 
                                  && leg2.hasTriggerMatch(TriggerEnum::HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1)
                                  && pt_1 > 25 && pt_2 > 35;

    } else if ( channel == HTTAnalysis::TauTau )
    {
        trg_singletau_leading= leg1.hasTriggerMatch(TriggerEnum::HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);
        trg_singletau_trailing= leg2.hasTriggerMatch(TriggerEnum::HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1);

        trg_doubletau_40_tightiso=          ( leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg) && leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg) ) && pt_1 > 45 && pt_2 > 45;
        trg_doubletau_40_mediso_tightid=    ( leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg) && leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg) ) && pt_1 > 45 && pt_2 > 45;
        trg_doubletau_35_tightiso_tightid=  ( leg1.hasTriggerMatch(TriggerEnum::HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg) && leg2.hasTriggerMatch(TriggerEnum::HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg) ) && pt_1 > 40 && pt_2 > 40;
    }

    trg_muonelectron=DEF; //fires HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL



    genPt_1=DEF;
    genPt_2=DEF;

    TLorentzVector ll=ev->getGenBosonP4(false);
    TLorentzVector llvis=ev->getGenBosonP4(true);
    gen_Mll=ll.M();
    gen_ll_px=ll.Px();
    gen_ll_py=ll.Py();
    gen_ll_pz=ll.Pz();
    gen_vis_Mll=llvis.M();
    gen_ll_vis_px=llvis.Px();
    gen_ll_vis_py=llvis.Py();
    gen_ll_vis_pz=llvis.Pz();

    topPtReweightWeightRun2=ev->getTopPtReWeight(false);
    topPtReweightWeightRun1=ev->getTopPtReWeight(true);
    zPtReweightWeight=ev->getZPtReWeight();

    zpt_weight_nom=DEFWEIGHT;
    zpt_weight_esup=DEFWEIGHT;
    zpt_weight_esdown=DEFWEIGHT;
    zpt_weight_ttup=DEFWEIGHT;
    zpt_weight_ttdown=DEFWEIGHT;
    zpt_weight_statpt0up=DEFWEIGHT;
    zpt_weight_statpt0down=DEFWEIGHT;
    zpt_weight_statpt40up=DEFWEIGHT;
    zpt_weight_statpt40down=DEFWEIGHT;
    zpt_weight_statpt80up=DEFWEIGHT;
    zpt_weight_statpt80down=DEFWEIGHT;

    eleTauFakeRateWeight=DEFWEIGHT;
    muTauFakeRateWeight=DEFWEIGHT;
    antilep_tauscaling=DEFWEIGHT;

    if(isMC)
    {
        xsec = ev->getXsec();
        genWeight = ev->getMCWeight();
        genNEventsWeight = ev->getGenNEventsWeight();
        puWeight = ev->getPUWeight();
        lumiWeight=1000*xsec*genNEventsWeight;
        NUP=ev->getLHEnOutPartons();
        fillStitchingWeight( ev->getSampleType() );
        fillScalefactors();
        fillLeptonFakeRateWeights();

    } else
    {
        effweight = DEFWEIGHT;
        xsec = DEFWEIGHT;
        genWeight = DEFWEIGHT;        
        genNEventsWeight = DEFWEIGHT;
        puWeight = DEFWEIGHT;
        weight = DEFWEIGHT;
        NUP=-1;
        lumiWeight=DEFWEIGHT;
    }

    // Make sure this weight is up to date. 
    // Trigger sf need to be applied according to what triggers are used
    weight = puWeight*stitchedWeight*genWeight*eleTauFakeRateWeight*muTauFakeRateWeight*topPtReweightWeightRun1*zPtReweightWeight*idisoweight_1*idisoweight_2*sf_trk*sf_reco;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillStitchingWeight(HTTEvent::sampleTypeEnum sampleType)
{
    if(sampleType == HTTEvent::DY)
    {
        if(NUP == 0) stitchedWeight = 0.0591296350328;
        if(NUP == 1) stitchedWeight = 0.0101562362959;
        if(NUP == 2) stitchedWeight = 0.0215030982864;
        if(NUP == 3) stitchedWeight = 0.0135063197418;
        if(NUP == 4) stitchedWeight = 0.00922459522375;

    } else if(sampleType == HTTEvent::WJets)
    {
        if(NUP == 0) stitchedWeight = 0.790555230048;
        if(NUP == 1) stitchedWeight = 0.150361226486;
        if(NUP == 2) stitchedWeight = 0.307137166339;
        if(NUP == 3) stitchedWeight = 0.0555843884964;
        if(NUP == 4) stitchedWeight = 0.052271728229;

    } else stitchedWeight=lumiWeight;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillLeg1Branches()
{
    UChar_t bitmask;

    pt_1=leg1P4.Pt();
    phi_1=leg1P4.Phi();
    eta_1=leg1P4.Eta();
    eta_SC_1=eta_1+leg1.getProperty(PropertyEnum::deltaEtaSC);
    m_1=leg1P4.M();
    q_1=leg1.getCharge();
    d0_1=leg1.getProperty(PropertyEnum::dxy);
    dZ_1=leg1.getProperty(PropertyEnum::dz);


    if      (pdg1==15) iso_1=leg1.getProperty( HTTEvent::usePropertyFor.at("tauIsolation") );
    else if (pdg1==13) iso_1=leg1.getProperty( HTTEvent::usePropertyFor.at("muonIsolation") );
    else if (pdg1==11) iso_1=leg1.getProperty( HTTEvent::usePropertyFor.at("electronIsolation") );;

    // bitmask=leg1.getProperty(PropertyEnum::idAntiEle);

    againstElectronMVA6_1 = leg1.getProperty(PropertyEnum::idAntiEle);
    againstElectronVLooseMVA6_1 =(againstElectronMVA6_1 & 0x1 )>0;
    againstElectronLooseMVA6_1  =(againstElectronMVA6_1 & 0x2 )>0;
    againstElectronMediumMVA6_1 =(againstElectronMVA6_1 & 0x4 )>0;
    againstElectronTightMVA6_1  =(againstElectronMVA6_1 & 0x8 )>0;
    againstElectronVTightMVA6_1 =(againstElectronMVA6_1 & 0x10)>0;

    againstMuon3_1=leg1.getProperty(PropertyEnum::idAntiMu);
    againstMuonLoose3_1=(againstMuon3_1 & 0x1)>0;
    againstMuonTight3_1=(againstMuon3_1 & 0x2)>0;

    byCombinedIsolationDeltaBetaCorrRaw3Hits_1=leg1.getProperty(PropertyEnum::rawIso);
    byLooseCombinedIsolationDeltaBetaCorr3Hits_1=DEF; //not in nanoAOD
    byMediumCombinedIsolationDeltaBetaCorr3Hits_1=DEF; //not in nanoAOD
    byTightCombinedIsolationDeltaBetaCorr3Hits_1=DEF; //not in nanoAOD
    // byIsolationMVA3newDMwoLTraw_1=DEF;
    // byIsolationMVA3newDMwLTraw_1=DEF;
    byIsolationMVA3oldDMwLTraw_1 =leg1.getProperty(HTTEvent::usePropertyFor.at("tauID")); //same as above!?

    byIsolationMVArun2017v2DBoldDMwLTraw2017_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauIsolation")); // Raw value
    byIsolationMVArun2017v2DBoldDMwLT2017_1=leg1.getProperty(HTTEvent::usePropertyFor.at("tauID"));    // Bitmask

    byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1 = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x1)>0;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x2)>0;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_1   = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x4)>0;;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x8)>0;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_1   = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x10)>0;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x20)>0;
    byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1  = (byIsolationMVArun2017v2DBoldDMwLT2017_1 & 0x40)>0;

    // byVLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byMediumIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byTightIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byVTightIsolationMVArun2v1DBnewDMwLT_1=DEF;

    // NewMVAIDVLoose_1=DEF; //?
    // NewMVAIDLoose_1=DEF;
    // NewMVAIDMedium_1=DEF;
    // NewMVAIDTight_1=DEF;
    // NewMVAIDVTight_1=DEF;
    // NewMVAIDVVTight_1=DEF;
    // idMVANewDM_1=DEF; //?

    chargedIsoPtSum_1=leg1.getProperty(PropertyEnum::chargedIso);
    neutralIsoPtSum_1=leg1.getProperty(PropertyEnum::neutralIso);
    puCorrPtSum_1=leg1.getProperty(PropertyEnum::puCorr);
    decayModeFindingOldDMs_1=leg1.getProperty(PropertyEnum::idDecayMode);
    decayMode_1=leg1.getProperty(PropertyEnum::decayMode);

    if (pdg1==13) id_m_loose_1=1; //already filtered at NanoAOD production
    id_m_medium_1=leg1.getProperty(PropertyEnum::mediumId);
    id_m_tight_1=leg1.getProperty(PropertyEnum::tightId);
    id_m_tightnovtx_1=DEF;
    bitmask=leg1.getProperty(PropertyEnum::highPtId);
    id_m_highpt_1=(bitmask & 0x2)>0;

    int intmask=leg1.getProperty(PropertyEnum::cutBased);
    id_e_cut_veto_1=intmask>=1;
    id_e_cut_loose_1=intmask>=2;
    id_e_cut_medium_1=intmask>=3;
    id_e_cut_tight_1=intmask>=4;
    id_e_mva_nt_loose_1=DEF;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillLeg2Branches()
{
    UChar_t bitmask;

    pt_2=leg2P4.Pt();
    phi_2=leg2P4.Phi();
    eta_2=leg2P4.Eta();
    m_2=leg2P4.M();
    q_2=leg2.getCharge();
    d0_2=leg2.getProperty(PropertyEnum::dxy);
    dZ_2=leg2.getProperty(PropertyEnum::dz);


    if      (pdg2==15) iso_2=leg2.getProperty( HTTEvent::usePropertyFor.at("tauIsolation") );
    else if (pdg2==13) iso_2=leg2.getProperty( HTTEvent::usePropertyFor.at("muonIsolation") );
    else if (pdg2==11) iso_2=leg2.getProperty( HTTEvent::usePropertyFor.at("electronIsolation") );;


    againstElectronMVA6_2=leg2.getProperty(PropertyEnum::idAntiEle);
    againstElectronVLooseMVA6_2 = (againstElectronMVA6_2 & 0x1)>0;
    againstElectronLooseMVA6_2  = (againstElectronMVA6_2 & 0x2)>0;
    againstElectronMediumMVA6_2 = (againstElectronMVA6_2 & 0x4)>0;
    againstElectronTightMVA6_2  = (againstElectronMVA6_2 & 0x8)>0;
    againstElectronVTightMVA6_2 = (againstElectronMVA6_2 & 0x10)>0;

    againstMuon3_2=leg2.getProperty(PropertyEnum::idAntiMu);
    againstMuonLoose3_2=(againstMuon3_2 & 0x1)>0;
    againstMuonTight3_2=(againstMuon3_2 & 0x2)>0;

    byCombinedIsolationDeltaBetaCorrRaw3Hits_2=leg2.getProperty(PropertyEnum::rawIso);
    byLooseCombinedIsolationDeltaBetaCorr3Hits_2=DEF; //not in nanoAOD
    byMediumCombinedIsolationDeltaBetaCorr3Hits_2=DEF; //not in nanoAOD
    byTightCombinedIsolationDeltaBetaCorr3Hits_2=DEF; //not in nanoAOD
    // byIsolationMVA3newDMwoLTraw_2=DEF;
    // byIsolationMVA3newDMwLTraw_2=DEF;
    byIsolationMVA3oldDMwLTraw_2 =leg2.getProperty(PropertyEnum::rawMVAoldDM2017v2); //same as above!?

    byIsolationMVArun2017v2DBoldDMwLTraw2017_2=leg2.getProperty( HTTEvent::usePropertyFor.at("tauIsolation") ); // Raw value
    byIsolationMVArun2017v2DBoldDMwLT2017_2=leg2.getProperty( HTTEvent::usePropertyFor.at("tauID") );     // Bitmask

    byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2 = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x1)>0;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x2)>0;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_2   = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x4)>0;;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x8)>0;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_2   = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x10)>0;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x20)>0;
    byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2  = (byIsolationMVArun2017v2DBoldDMwLT2017_2 & 0x40)>0;

    // byVLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byMediumIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byVTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // NewMVAIDVLoose_2=DEF;
    // NewMVAIDLoose_2=DEF;
    // NewMVAIDMedium_2=DEF;
    // NewMVAIDTight_2=DEF;
    // NewMVAIDVTight_2=DEF;
    // NewMVAIDVVTight_2=DEF;
    // idMVANewDM_2=DEF;

    chargedIsoPtSum_2=leg2.getProperty(PropertyEnum::chargedIso);
    neutralIsoPtSum_2=leg2.getProperty(PropertyEnum::neutralIso);
    puCorrPtSum_2=leg2.getProperty(PropertyEnum::puCorr);
    decayModeFindingOldDMs_2=leg2.getProperty(PropertyEnum::idDecayMode);
    decayMode_2=leg2.getProperty(PropertyEnum::decayMode);

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillJetBranches(HTTJetCollection *jets)
{
    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {

        jets->setCurrentUncertShift( jecShifts[shift].second.first, jecShifts[shift].second.second);
        njets[shift]        = jets->getNJets(30);
        njetspt20[shift]    = jets->getNJets(20);
        njetingap[shift]    = jets->getNJetInGap(30);
        njetingap20[shift]  = jets->getNJetInGap(20);        


        if ( njetspt20[shift] >=1 )
        {
            jpt_1[shift]   = jets->getJet(0).Pt();
            jeta_1[shift]  = jets->getJet(0).Eta();
            jphi_1[shift]  = jets->getJet(0).Phi();
            if( strcmp(jecShifts[shift].first.c_str(), "") == 0 )
            {
                jm_1=   jets->getJet(0).M();
                jrawf_1=jets->getJet(0).getProperty(PropertyEnum::rawFactor);
                jmva_1= jets->getJet(0).getProperty(PropertyEnum::btagCMVA);
                jcsv_1= jets->getJet(0).getProperty(PropertyEnum::btagCSVV2);
            }
        }
        if ( njetspt20[shift] >=2 )
        {

            jpt_2[shift]   =  jets->getJet(1).Pt();
            jeta_2[shift]  = jets->getJet(1).Eta();
            jphi_2[shift]  = jets->getJet(1).Phi();
            
            mjj[shift]      = jets->getDiJetP4().M();
            dijetpt[shift]  = jets->getDiJetP4().Pt();
            dijetphi[shift] = jets->getDiJetP4().Phi();
            jdeta[shift]    = std::abs( jets->getJet(0).Eta()-jets->getJet(1).Eta() );
            jdphi[shift]    = jets->getJet(0).P4().DeltaPhi( jets->getJet(1).P4() );

            if( strcmp(jecShifts[shift].first.c_str(), "") == 0 )
            {
                jm_2=   jets->getJet(1).M();
                jrawf_2=jets->getJet(1).getProperty(PropertyEnum::rawFactor);
                jmva_2= jets->getJet(1).getProperty(PropertyEnum::btagCMVA);
                jcsv_2= jets->getJet(1).getProperty(PropertyEnum::btagCSVV2);
                jeta1eta2=jeta_1[shift]*jeta_2[shift];
                lep_etacentrality=TMath::Exp( -4/pow(jeta_1[shift]-jeta_2[shift],2) * pow( (eta_1-( jeta_1[shift]+jeta_2[shift] )*0.5), 2 ) );
            }
        }
    }
    jets->setCurrentUncertShift("",true);
    for(unsigned int shift = 0; shift<btagShifts.size(); ++shift )
    {
        jets->btagPromoteDemote( btagShifts[shift].second.first, btagShifts[shift].second.second );
        nbtag[shift]= jets->getNBtag();
        if (nbtag[shift] >= 1)
        {
            bpt_1[shift]=   jets->getBtagJet(0).Pt();
            beta_1[shift]=  jets->getBtagJet(0).Eta();
            bphi_1[shift]=  jets->getBtagJet(0).Phi();
            brawf_1[shift]= jets->getBtagJet(0).getProperty(PropertyEnum::rawFactor);
            bmva_1[shift]=  jets->getBtagJet(0).getProperty(PropertyEnum::btagCMVA);
            bcsv_1[shift]=  jets->getBtagJet(0).getProperty(PropertyEnum::btagCSVV2);
        }

        if (nbtag[shift] >= 2)
        {
            bpt_2[shift]=   jets->getBtagJet(1).Pt();
            beta_2[shift]=  jets->getBtagJet(1).Eta();
            bphi_2[shift]=  jets->getBtagJet(1).Phi();
            brawf_2[shift]= jets->getBtagJet(1).getProperty(PropertyEnum::rawFactor);
            bmva_2[shift]=  jets->getBtagJet(1).getProperty(PropertyEnum::btagCMVA);
            bcsv_2[shift]=  jets->getBtagJet(1).getProperty(PropertyEnum::btagCSVV2);
        }


        if (pdg1==15) gen_match_jetId_1=getGenMatch_jetId(leg1P4,jets);
        if (pdg2==15) gen_match_jetId_2=getGenMatch_jetId(leg2P4,jets);
    } 

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillPairBranches(HTTPair *pair, HTTJetCollection *jets)
{


    metcov00=pair->getMETMatrix().at(0);
    metcov01=pair->getMETMatrix().at(1);
    metcov10=pair->getMETMatrix().at(2);
    metcov11=pair->getMETMatrix().at(3);

    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {
        jets->setCurrentUncertShift( jecShifts[shift].second.first, jecShifts[shift].second.second );
        pair->setCurrentMETShift( jecShifts[shift].first );

        met[shift]   =pair->getMET().Mod();
        met_ex[shift]=pair->getMET().X();
        met_ey[shift]=pair->getMET().Y();
        metphi[shift]= TVector2::Phi_mpi_pi( pair->getMET().Phi() );

        mt_1[shift]=pair->getMTLeg1();
        mt_2[shift]=pair->getMTLeg2();

        m_sv[shift]=pair->getP4(HTTPair::SVFit).M();
        pt_sv[shift]=pair->getP4(HTTPair::SVFit).Pt();

        pt_tt[shift]=pair->getP4( HTTPair::Simple ).Pt();
        mt_tot[shift]=pair->getMTTOT();
        pt_sum[shift]=pt_1+pt_2+met[shift];

        htxs_reco_ggf[shift] = (int)getStage1Category(HiggsProdMode::GGF, pair->getP4( HTTPair::Simple ), jets);
        htxs_reco_vbf[shift] = (int)getStage1Category(HiggsProdMode::VBF, pair->getP4( HTTPair::Simple ), jets);
        //////////////////////////////////////////////////////////////////
        TLorentzVector ttjj; ttjj.SetPtEtaPhiM(-10,-10,-10,-10);
        if ( njetspt20[0] >=2 )
        {
            ttjj = pair->getP4( HTTPair::Simple ) + jets->getJet(0).P4() + jets->getJet(1).P4();
        }
        pt_ttjj[shift] = ttjj.Pt();
        m_ttjj[shift] = ttjj.M();
        //////////////////////////////////////////////////////////////////
    }
    jets->setCurrentUncertShift( "", true );
    pair->setCurrentMETShift( "" );
    //////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////
    pt_vis=pair->getPTVis();
    m_vis=pair->getMVis();
    dphi=leg1P4.DeltaPhi(leg2P4);
    dr_leptau=leg1P4.DeltaR(leg2P4);    
    //////////////////////////////////////////////////////////////////
    double zetaX = TMath::Cos(phi_1) + TMath::Cos(phi_2);
    double zetaY = TMath::Sin(phi_1) + TMath::Sin(phi_2);
    double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
    if ( zetaR > 0. ) {
      zetaX /= zetaR;
      zetaY /= zetaR;
    }
    pzetavis =  (leg1P4.Px() + leg2P4.Px())*zetaX + (leg1P4.Py() + leg2P4.Py())*zetaY;
    pzetamiss = met_ex[0]*zetaX + met_ey[0]*zetaY;
    dzeta = pzetamiss + (pzetavis - 1.85 * pzetavis);
    //////////////////////////////////////////////////////////////////
    double x1 = pt_1 / ( pt_1 + met[0] );
    double x2 = pt_2 / ( pt_2 + met[0] );
    if( TMath::Cos( dphi ) > -0.95
        && ( x1 > 0 && x1 < 1)
        && ( x2 > 0 && x2 < 1) )
    {
        m_coll =  m_vis / sqrt( x1 * x2 ) ;
    }
    else m_coll = DEF;
    //////////////////////////////////////////////////////////////////
    vector<TLorentzVector> objs;
    objs.push_back(leg1P4);
    objs.push_back(leg2P4);
    if ( njetspt20[0]>0 ) objs.push_back(jets->getJet(0).P4());
    if ( njetspt20[0]>1 ) objs.push_back(jets->getJet(1).P4());
    sphericity=calcSphericity(objs);  
    //////////////////////////////////////////////////////////////////
    TLorentzVector v0=leg1P4*( 1/sqrt(  pow(leg1P4.Px(),2)+pow(leg1P4.Py(),2)  ) ); //lep, normalized in transverse plane
    TLorentzVector v1=leg2P4*( 1/sqrt(  pow(leg2P4.Px(),2)+pow(leg2P4.Py(),2)  ) ); //tau, normalized in transverse plane
    float omega=v1.DeltaPhi(v0);
    float theta=-v0.Phi();
    float x=(     met_ex[0] * TMath::Sin(omega-theta)  - met_ey[0]*TMath::Cos(omega-theta)   ) / TMath::Sin(omega); //x coord in lep-tau system
    float y=(     met_ex[0] * TMath::Sin(theta)        + met_ey[0]*TMath::Cos(theta)         ) / TMath::Sin(omega); //y coord in lep-tau system
    met_centrality=( x+y ) / sqrt(x*x + y*y);
    //////////////////////////////////////////////////////////////////
 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillAdditionalLeptons( std::vector<HTTParticle> leptons, HTTPair *pair )
{
    unsigned int indexLeg1 = (unsigned int)pair->getIndexLeg1();
    unsigned int indexLeg2 = (unsigned int)pair->getIndexLeg2();

    UChar_t bitmask;

    for(unsigned int i = 0; i<leptons.size(); i++)
    {
        if( i == indexLeg1 || i == indexLeg2 || !leptons[i].isAdditionalLepton() ) continue;

        int lepPDGId = leptons[i].getPDGid();
        if( std::abs(lepPDGId) == 15 &&  channel == HTTAnalysis::TauTau ) continue;

        addlepton_p4.push_back( leptons[i].getP4() );
        addlepton_pt.push_back( leptons[i].getP4().Pt() );
        addlepton_eta.push_back( leptons[i].getP4().Eta() );
        addlepton_phi.push_back( leptons[i].getP4().Phi() );
        addlepton_m.push_back( leptons[i].getP4().M() );        
        addlepton_pdgId.push_back( leptons[i].getPDGid()*(-1) );
        addlepton_mc_match.push_back( leptons[i].getProperty(PropertyEnum::mc_match) );

        addlepton_tauAntiEle.push_back( (int)leptons[i].getProperty(PropertyEnum::idAntiEle) );
        addlepton_tauAntiMu.push_back( (int)leptons[i].getProperty(PropertyEnum::idAntiMu)  );

        addlepton_d0.push_back( leptons[i].getProperty(PropertyEnum::dxy) );
        addlepton_dZ.push_back( leptons[i].getProperty(PropertyEnum::dz) );
        addlepton_mt.push_back( leptons[i].getMT( pair->getMET() ) );

        if( std::abs(lepPDGId) == 13 )
        {
            nadditionalMu++;
            addmuon_pt.push_back( leptons[i].getP4().Pt() );
            addmuon_eta.push_back( leptons[i].getP4().Eta() );
            addmuon_phi.push_back( leptons[i].getP4().Phi() );
            addmuon_m.push_back( leptons[i].getP4().M() );
            addmuon_q.push_back(leptons[i].getProperty(PropertyEnum::charge));
            addmuon_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("muonIsolation") ) );
            addmuon_gen_match.push_back(leptons[i].getProperty(PropertyEnum::mc_match) );

            addlepton_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("muonIsolation") ) );
            addlepton_mvis.push_back( ( leg2P4 + leptons[i].getP4() ).M() );

            addlepton_tauID.push_back( 0. );
            addlepton_tauDM.push_back( -1. );
            addlepton_tauCombIso.push_back( 0. );

        }else if( std::abs(lepPDGId) == 11 )
        {
            nadditionalEle++;
            addele_pt.push_back( leptons[i].getP4().Pt() );
            addele_eta.push_back( leptons[i].getP4().Eta() );
            addele_phi.push_back( leptons[i].getP4().Phi() );
            addele_m.push_back( leptons[i].getP4().M() );
            addele_q.push_back(leptons[i].getProperty(PropertyEnum::charge));
            addele_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("electronIsolation") ) );
            addele_gen_match.push_back(leptons[i].getProperty(PropertyEnum::mc_match) );

            addlepton_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("electronIsolation") ) );
            addlepton_mvis.push_back( ( leg2P4 + leptons[i].getP4() ).M() );

            addlepton_tauID.push_back( 0. );
            addlepton_tauDM.push_back( -1. );
            addlepton_tauCombIso.push_back( 0. );


        }else if( std::abs(lepPDGId) == 15 )
        {
            nadditionalTau++;
            addtau_pt.push_back( leptons[i].getP4().Pt() );
            addtau_eta.push_back( leptons[i].getP4().Eta() );
            addtau_phi.push_back( leptons[i].getP4().Phi() );
            addtau_m.push_back( leptons[i].getP4().M() );
            addtau_q.push_back(leptons[i].getProperty(PropertyEnum::charge));
            addtau_byIsolationMVArun2v1DBoldDMwLTraw.push_back(leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolation") ));
            addtau_gen_match.push_back( leptons[i].getProperty(PropertyEnum::mc_match) );

            bitmask = leptons[i].getProperty(PropertyEnum::idAntiEle);
            bool antiEle = ( bitmask & 0x8) > 0;
            bitmask = leptons[i].getProperty(PropertyEnum::idAntiMu);
            bool antiMu  = (bitmask & 0x1 ) > 0;
            addtau_passesTauLepVetos.push_back( antiEle & antiMu );

            addtau_decayMode.push_back( leptons[i].getProperty(PropertyEnum::decayMode) );
            addtau_d0.push_back( leptons[i].getProperty(PropertyEnum::dxy) );
            addtau_dZ.push_back( leptons[i].getProperty(PropertyEnum::dz) );

            addtau_mt.push_back( leptons[i].getMT( pair->getMET() ) );
            addtau_mvis.push_back( ( leg1P4 + leptons[i].getP4() ).M() );

            addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits.push_back( leptons[i].getProperty(PropertyEnum::rawIso) );

            bitmask=leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauID"));
            addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x1)>0 );
            addtau_byVLooseIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x2)>0 );
            addtau_byLooseIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x4)>0 );
            addtau_byMediumIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x8)>0 );
            addtau_byTightIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x10)>0 );
            addtau_byVTightIsolationMVArun2v1DBoldDMwLT.push_back( (bitmask & 0x20)>0 );

            addlepton_iso.push_back( leptons[i].getProperty( HTTEvent::usePropertyFor.at("tauIsolation") ) );
            addlepton_mvis.push_back( ( leg1P4 + leptons[i].getP4() ).M() );

            addlepton_tauID.push_back(  (int)leptons[i].getProperty(HTTEvent::usePropertyFor.at("tauID")) );
            addlepton_tauDM.push_back(  (int)leptons[i].getProperty(PropertyEnum::decayMode) );
            addlepton_tauCombIso.push_back( leptons[i].getProperty(PropertyEnum::rawIso) );


        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillLeptonFakeRateWeights()
{
    eleTauFakeRateWeight = 1.0;
    muTauFakeRateWeight = 1.0;
    antilep_tauscaling = 1.0;

    if(channel == HTTAnalysis::MuTau)
    {

        if((gen_match_2 == 1 || gen_match_2 == 3) && againstElectronVLooseMVA6_2 > 0.5 )
        {
            if( std::abs(eta_2) < 1.460 )       eleTauFakeRateWeight *= 1.09;
            else if ( std::abs(eta_2) > 1.558 ) eleTauFakeRateWeight *= 1.19;
        }
        if((gen_match_2 == 2 || gen_match_2 == 4) && againstMuonTight3_2 > 0.5 )
        {
            if( std::abs(eta_2) < 0.4 )       muTauFakeRateWeight *= 1.17;
            else if ( std::abs(eta_2) < 0.8 ) muTauFakeRateWeight *= 1.29;
            else if ( std::abs(eta_2) < 1.2 ) muTauFakeRateWeight *= 1.14;
            else if ( std::abs(eta_2) < 1.7 ) muTauFakeRateWeight *= 0.93;
            else if ( std::abs(eta_2) < 2.3 ) muTauFakeRateWeight *= 1.61;
        }        
    }
    if(channel == HTTAnalysis::EleTau)
    {
        if((gen_match_2 == 1 || gen_match_2 == 3) && againstElectronTightMVA6_2 > 0.5 )
        {
            if( std::abs(eta_2) < 1.460 )       eleTauFakeRateWeight *= 1.80;
            else if ( std::abs(eta_2) > 1.558 ) eleTauFakeRateWeight *= 1.53;
        }
        if((gen_match_2 == 2 || gen_match_2 == 4) && againstMuonLoose3_2> 0.5 )
        {
            if( std::abs(eta_2) < 0.4 )       muTauFakeRateWeight *= 1.06;
            else if ( std::abs(eta_2) < 0.8 ) muTauFakeRateWeight *= 1.02;
            else if ( std::abs(eta_2) < 1.2 ) muTauFakeRateWeight *= 1.10;
            else if ( std::abs(eta_2) < 1.7 ) muTauFakeRateWeight *= 1.03;
            else if ( std::abs(eta_2) < 2.3 ) muTauFakeRateWeight *= 1.94;
        }        
    }
    if(channel == HTTAnalysis::TauTau)
    {

        if((gen_match_1 == 1 || gen_match_1 == 3) && againstElectronVLooseMVA6_1 > 0.5 )
        {
            if( std::abs(eta_1) < 1.460 )       eleTauFakeRateWeight *= 1.09;
            else if ( std::abs(eta_1) > 1.558 ) eleTauFakeRateWeight *= 1.19;
        }
        if((gen_match_1 == 2 || gen_match_1 == 4) && againstMuonLoose3_1> 0.5 )
        {
            if( std::abs(eta_1) < 0.4 )       muTauFakeRateWeight *= 1.06;
            else if ( std::abs(eta_1) < 0.8 ) muTauFakeRateWeight *= 1.02;
            else if ( std::abs(eta_1) < 1.2 ) muTauFakeRateWeight *= 1.10;
            else if ( std::abs(eta_1) < 1.7 ) muTauFakeRateWeight *= 1.03;
            else if ( std::abs(eta_1) < 2.3 ) muTauFakeRateWeight *= 1.94;
        } 

        if((gen_match_2 == 1 || gen_match_2 == 3) && againstElectronVLooseMVA6_2 > 0.5 )
        {
            if( std::abs(eta_2) < 1.460 )       eleTauFakeRateWeight *= 1.09;
            else if ( std::abs(eta_2) > 1.558 ) eleTauFakeRateWeight *= 1.19;
        }
        if((gen_match_2 == 2 || gen_match_2 == 4) && againstMuonLoose3_2> 0.5 )
        {
            if( std::abs(eta_2) < 0.4 )       muTauFakeRateWeight *= 1.06;
            else if ( std::abs(eta_2) < 0.8 ) muTauFakeRateWeight *= 1.02;
            else if ( std::abs(eta_2) < 1.2 ) muTauFakeRateWeight *= 1.10;
            else if ( std::abs(eta_2) < 1.7 ) muTauFakeRateWeight *= 1.03;
            else if ( std::abs(eta_2) < 2.3 ) muTauFakeRateWeight *= 1.94;
        }        
    }
    antilep_tauscaling = eleTauFakeRateWeight * muTauFakeRateWeight;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::fillScalefactors()
{
    // from https://github.com/CMS-HTT/CorrectionsWorkspace/tree/2017_17NovReRecoData_Fall17MC
    singleTriggerSFLeg1 = DEFWEIGHT;
    singleTriggerSFLeg2 = DEFWEIGHT;

    double s1_data = DEFWEIGHT;
    double s1_mc   = DEFWEIGHT;
    double x1_data = DEFWEIGHT;
    double x1_mc   = DEFWEIGHT;
    double x2_data = DEFWEIGHT;
    double x2_mc   = DEFWEIGHT;

    isoWeight_1 = DEFWEIGHT;
    isoWeight_2 = DEFWEIGHT;
    idWeight_1 = DEFWEIGHT;
    idWeight_2 = DEFWEIGHT;

    idisoweight_1 = DEFWEIGHT;
    idisoweight_2 = DEFWEIGHT;

    // from https://github.com/truggles/TauTriggerSFs2017/tree/tauTriggers2017_MCv2_PreReMiniaod
    xTriggerSFLeg1 = DEFWEIGHT;
    xTriggerSFLeg2 = DEFWEIGHT;

    effweight = DEFWEIGHT;

    sf_trk =DEFWEIGHT;
    sf_reco=DEFWEIGHT;
    sf_SingleOrCrossTrigger = DEFWEIGHT;
    sf_SingleXorCrossTrigger = DEFWEIGHT;
    sf_SingleTrigger = DEFWEIGHT;
    sf_DoubleTauTight = DEFWEIGHT;
    sf_DoubleTauVTight = DEFWEIGHT;

    if( channel == HTTAnalysis::MuTau )
    {
        w->var("m_pt")->setVal(  pt_1  );
        w->var("m_eta")->setVal( eta_1 );

        singleTriggerSFLeg1 = w->function("m_trg_SingleMu_Mu24ORMu27_desy_ratio")->getVal();
        s1_data = w->function("m_trg_SingleMu_Mu24ORMu27_desy_data")->getVal();
        s1_mc   = w->function("m_trg_SingleMu_Mu24ORMu27_desy_mc")->getVal();

        if( std::abs(eta_1) < 2.1 )
        {
            xTriggerSFLeg1 = w->function("m_trg_MuTau_Mu20Leg_desy_ratio")->getVal();
            x1_data = w->function("m_trg_MuTau_Mu20Leg_desy_data")->getVal();
            x1_mc   = w->function("m_trg_MuTau_Mu20Leg_desy_mc")->getVal();
        }
        if( std::abs(eta_2) < 2.1 )
        {
            xTriggerSFLeg2 = tauTrigSFTight->getMuTauScaleFactor( pt_2 ,  eta_2 ,  phi_2 );
            x2_data        = tauTrigSFTight->getMuTauEfficiencyData( pt_2 ,  eta_2 ,  phi_2 );
            x2_mc          = tauTrigSFTight->getMuTauEfficiencyMC( pt_2 ,  eta_2 ,  phi_2 );
        }

        idWeight_1  = w->function("m_id_ratio")->getVal();
        isoWeight_1 = w->function("m_iso_ratio")->getVal();
        idisoweight_1 =  idWeight_1 * isoWeight_1 ;

        sf_trk = w->function("m_trk_ratio")->getVal();
        sf_SingleOrCrossTrigger = (s1_data*(1 - x2_data ) + x1_data*x2_data ) / (s1_mc*(1 - x2_mc ) + x1_mc*x2_mc );
        sf_SingleXorCrossTrigger = singleTriggerSFLeg1*xTriggerSFLeg2;
        sf_SingleTrigger = singleTriggerSFLeg1;
    }

    if( channel == HTTAnalysis::EleTau )
    {
        w->var("e_pt")->setVal(  pt_1  );
        w->var("e_eta")->setVal( eta_1 );

        singleTriggerSFLeg1 = w->function("e_trg_27_32_35_ratio")->getVal();
        s1_data = w->function("e_trg_27_32_35_data")->getVal();
        s1_mc   = w->function("e_trg_27_32_35_mc")->getVal();        

        if( std::abs(eta_1) < 2.1 )
        {
            xTriggerSFLeg1 = w->function("e_trg24_ratio")->getVal();
            x1_data = w->function("e_trg24_data")->getVal();
            x1_mc   = w->function("e_trg24_mc")->getVal();
        }       
        if( std::abs(eta_2) < 2.1 )
        {
            xTriggerSFLeg2 = tauTrigSFTight->getETauScaleFactor( pt_2 ,  eta_2 ,  phi_2 );
            x2_data        = tauTrigSFTight->getETauEfficiencyData( pt_2 ,  eta_2 ,  phi_2 );
            x2_mc          = tauTrigSFTight->getETauEfficiencyMC( pt_2 ,  eta_2 ,  phi_2 );
        }    

        idWeight_1  = w->function("e_id_ratio")->getVal();
        isoWeight_1 = w->function("e_iso_ratio")->getVal();
        idisoweight_1 =  idWeight_1 * isoWeight_1 ;

        sf_reco = w->function("e_reco_ratio")->getVal();
        sf_SingleOrCrossTrigger = (s1_data*(1 - x2_data ) + x1_data*x2_data ) / (s1_mc*(1 - x2_mc ) + x1_mc*x2_mc );
        sf_SingleXorCrossTrigger = singleTriggerSFLeg1*xTriggerSFLeg1*xTriggerSFLeg2;
        sf_SingleTrigger = singleTriggerSFLeg1;
    }

    if( channel == HTTAnalysis::TauTau )
    {
        xTriggerSFLeg1 = tauTrigSFTight->getDiTauScaleFactor(  pt_1 ,  eta_1 ,  phi_1  );
        xTriggerSFLeg2 = tauTrigSFTight->getDiTauScaleFactor(  pt_2 ,  eta_2 ,  phi_2  );

        sf_DoubleTauTight = xTriggerSFLeg1*xTriggerSFLeg2;

        xTriggerSFLeg1 = tauTrigSFVTight->getDiTauScaleFactor(  pt_1 ,  eta_1 ,  phi_1  );
        xTriggerSFLeg2 = tauTrigSFVTight->getDiTauScaleFactor(  pt_2 ,  eta_2 ,  phi_2  );

        sf_DoubleTauVTight = xTriggerSFLeg1*xTriggerSFLeg2;
    }

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EventWriter::calcSphericity(std::vector<TLorentzVector> p){

  TMatrixD S(3,3);

  std::vector<std::vector<double> > pvec;

  float denom=0;
  for (unsigned k=0; k<p.size(); k++){
    double dtmp[3]={p.at(k).Px(),p.at(k).Py(),p.at(k).Pz()};
    std::vector<double> vtmp(dtmp, dtmp+3);
    pvec.push_back(vtmp);
    denom+=dtmp[0]*dtmp[0]+dtmp[1]*dtmp[1]+dtmp[2]*dtmp[2];
  }

  for (int i=0; i<3; i++){
    for (int j=0; j<3; j++){
      float num=0;
      for (unsigned k=0; k<pvec.size(); k++){
        num+=pvec.at(k)[i]*pvec.at(k)[j];
      }
      S(i,j)=num/denom;
    }
  }
  return calcSphericityFromMatrix(S);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EventWriter::calcSphericityFromMatrix(TMatrixD M) {

  //  TMatrixD M(3,3);
  //  M.SetMatrixArray(A);

  TMatrixDEigen V(M);
  TMatrixD Eval = V.GetEigenValues();

  double e1 = TMatrixDRow(Eval,0)(0);
  double e2 = TMatrixDRow(Eval,1)(1);
  double e3 = TMatrixDRow(Eval,2)(2);

  std::vector<double> evalvector;
  evalvector.push_back(e1);
  evalvector.push_back(e2);
  evalvector.push_back(e3);

  //sort eigenvalues to get the lowest two for use in sphericity (lowest to highest order)
  sort (evalvector.begin(), evalvector.end());

  //error-checking
  //this number should equal zero as the off-diagonal elements should all be zero in the eigenvalue matrix
  //returns error value of -1
  double check = TMatrixDRow(Eval,0)(1) + TMatrixDRow(Eval,0)(2) + TMatrixDRow(Eval,1)(0) + TMatrixDRow(Eval,1)(2) + TMatrixDRow(Eval,2)(0) + TMatrixDRow(Eval,2)(1);
  if (check != 0.0) {double err = -1; return err;}

  //for formula, see: http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node213.html
  
  double value = evalvector.at(0)+evalvector.at(1);
  double spher = 1.5*value;

  return spher;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int EventWriter::getGenMatch_jetId(TLorentzVector selObj, HTTJetCollection *jets){
  float minDR=1;
  int whichjet=0;

  for (unsigned i=0; i<(unsigned int)jets->getNJets(20); i++){
    TLorentzVector p4=jets->getJet(i).P4();
    if(p4.Pt() > 20 && fabs(p4.Eta() ) < 4.7 ){
      float tmpDR = calcDR( selObj.Eta(), selObj.Phi(), p4.Eta(), p4.Phi() );
      if( tmpDR < minDR ){
    minDR = tmpDR;
    whichjet=i;
      }
    }
  }

  if( minDR < 0.5 ) return jets->getJet(whichjet).getProperty(PropertyEnum::partonFlavour);
  return -99;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double EventWriter::calcDR(double eta1, double phi1, double eta2, double phi2){
  double deta = eta1-eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  return TMath::Sqrt( deta*deta+dphi*dphi );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::setDefault(){

    lumiWeight=DEFWEIGHT;
    run_syncro=DEF;
    lumi_syncro=DEF;
    evt_syncro=0; //unsigned
    entry=DEF;
    fileEntry=DEF;

    npv=DEF;
    npvGood=DEF;
    npu=DEF;
    rho=DEF;
    gen_match_1=DEF;
    gen_match_2=DEF;
    genPt_1=DEF;
    genPt_2=DEF;
    gen_match_jetId_1=DEF;
    gen_match_jetId_2=DEF;
    NUP=DEF;
    //evtWeight=DEF;

    singleTriggerSFLeg1 = DEFWEIGHT;
    singleTriggerSFLeg2 = DEFWEIGHT;
    xTriggerSFLeg1 = DEFWEIGHT;
    xTriggerSFLeg2 = DEFWEIGHT;

    isoWeight_1 = DEFWEIGHT;
    isoWeight_2 = DEFWEIGHT;
    idWeight_1 = DEFWEIGHT;
    idWeight_2 = DEFWEIGHT;

    weight=DEFWEIGHT;
    puWeight=DEFWEIGHT;
    genWeight=DEFWEIGHT;
    trigweight_1=DEFWEIGHT;
    anti_trigweight_1=DEFWEIGHT;
    trigweight_2=DEFWEIGHT;
    idisoweight_1=DEFWEIGHT;
    anti_idisoweight_1=DEFWEIGHT;
    idisoweight_2=DEFWEIGHT;
    sf_trk=DEFWEIGHT;
    sf_reco=DEFWEIGHT;
    effweight=DEFWEIGHT;
    sf_SingleOrCrossTrigger=DEFWEIGHT;
    stitchedWeight=DEFWEIGHT;
    topPtReweightWeightRun2=DEFWEIGHT;
    topPtReweightWeightRun1=DEFWEIGHT;
    zPtReweightWeight=DEFWEIGHT;
    zpt_weight_nom=DEFWEIGHT;
    zpt_weight_esup=DEFWEIGHT;
    zpt_weight_esdown=DEFWEIGHT;
    zpt_weight_ttup=DEFWEIGHT;
    zpt_weight_ttdown=DEFWEIGHT;
    zpt_weight_statpt0up=DEFWEIGHT;
    zpt_weight_statpt0down=DEFWEIGHT;
    zpt_weight_statpt40up=DEFWEIGHT;
    zpt_weight_statpt40down=DEFWEIGHT;
    zpt_weight_statpt80up=DEFWEIGHT;
    zpt_weight_statpt80down=DEFWEIGHT;

    
    eleTauFakeRateWeight=DEFWEIGHT;
    muTauFakeRateWeight=DEFWEIGHT;
    antilep_tauscaling=DEFWEIGHT;
//////////////////////////////////////////////////////////////////
    gen_Mll=DEF;
    gen_ll_px=DEF;
    gen_ll_py=DEF;
    gen_ll_pz=DEF;
    gen_vis_Mll=DEF;
    gen_ll_vis_px=DEF;
    gen_ll_vis_py=DEF;
    gen_ll_vis_pz=DEF;
    gen_top_pt_1=DEF;
    gen_top_pt_2=DEF;
    genJets=DEF;
    genJet_match_1=DEF;
    genJet_match_2=DEF;
    //////////////////////////////////////////////////////////////////  
    trg_singletau_leading=false;
    trg_singletau_trailing=false;
    trg_singlemuon_27=false;
    trg_singlemuon_24=false;
    trg_crossmuon_mu20tau27=false;
    trg_singleelectron_35=false;
    trg_singleelectron_32=false;
    trg_singleelectron_27=false;
    trg_crossele_ele24tau30=false;
    trg_doubletau_40_tightiso=false;
    trg_doubletau_40_mediso_tightid=false;
    trg_doubletau_35_tightiso_tightid=false;
    trg_muonelectron=DEF; //fires HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
    Flag_HBHENoiseFilter=DEF;
    Flag_HBHENoiseIsoFilter=DEF;
    Flag_EcalDeadCellTriggerPrimitiveFilter=DEF;
    Flag_goodVertices=DEF;
    Flag_eeBadScFilter=DEF;
    Flag_globalTightHalo2016Filter=DEF;
    failBadGlobalMuonTagger=DEF;
    failCloneGlobalMuonTagger=DEF;

    //////////////////////////////////////////////////////////////////  
    pt_1=DEF;
    phi_1=DEF;
    eta_1=DEF;
    eta_SC_1=DEF;
    m_1=DEF;
    q_1=DEF;
    d0_1=DEF;
    dZ_1=DEF;
    iso_1=DEF;
    againstElectronLooseMVA6_1=DEF;
    againstElectronMediumMVA6_1=DEF;
    againstElectronTightMVA6_1=DEF;
    againstElectronVLooseMVA6_1=DEF;
    againstElectronVTightMVA6_1=DEF;
    againstMuonLoose3_1=DEF;
    againstMuonTight3_1=DEF;
    byCombinedIsolationDeltaBetaCorrRaw3Hits_1=DEF;
    byLooseCombinedIsolationDeltaBetaCorr3Hits_1=DEF;
    byMediumCombinedIsolationDeltaBetaCorr3Hits_1=DEF;
    byTightCombinedIsolationDeltaBetaCorr3Hits_1=DEF;
    byIsolationMVA3newDMwoLTraw_1=DEF;
    byIsolationMVA3oldDMwoLTraw_1=DEF;
    byIsolationMVA3newDMwLTraw_1=DEF;
    byIsolationMVA3oldDMwLTraw_1=DEF;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_1=DEF;
    // byVLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byMediumIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byTightIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // byVTightIsolationMVArun2v1DBnewDMwLT_1=DEF;
    // NewMVAIDVLoose_1=DEF;
    // NewMVAIDLoose_1=DEF;
    // NewMVAIDMedium_1=DEF;
    // NewMVAIDTight_1=DEF;
    // NewMVAIDVTight_1=DEF;
    // NewMVAIDVVTight_1=DEF;
    // idMVANewDM_1=DEF;
    chargedIsoPtSum_1=DEF;
    neutralIsoPtSum_1=DEF;
    puCorrPtSum_1=DEF;
    decayModeFindingOldDMs_1=DEF;
    decayMode_1=DEF;
    id_e_mva_nt_loose_1=DEF;
    id_m_loose_1=DEF;
    id_m_medium_1=DEF;
    id_m_tight_1=DEF;
    id_m_tightnovtx_1=DEF;
    id_m_highpt_1=DEF;
    id_e_cut_veto_1=DEF;
    id_e_cut_loose_1=DEF;
    id_e_cut_medium_1=DEF;
    id_e_cut_tight_1=DEF;
    
    //////////////////////////////////////////////////////////////////
    pt_2=DEF;
    phi_2=DEF;
    eta_2=DEF;
    m_2=DEF;
    q_2=DEF;
    d0_2=DEF;
    dZ_2=DEF;
    iso_2=DEF;
    againstElectronLooseMVA6_2=DEF;
    againstElectronMediumMVA6_2=DEF;
    againstElectronTightMVA6_2=DEF;
    againstElectronVLooseMVA6_2=DEF;
    againstElectronVTightMVA6_2=DEF;
    againstMuonLoose3_2=DEF;
    againstMuonTight3_2=DEF;
    byCombinedIsolationDeltaBetaCorrRaw3Hits_2=DEF;
    byLooseCombinedIsolationDeltaBetaCorr3Hits_2=DEF;
    byMediumCombinedIsolationDeltaBetaCorr3Hits_2=DEF;
    byTightCombinedIsolationDeltaBetaCorr3Hits_2=DEF;
    byIsolationMVA3newDMwoLTraw_2=DEF;
    byIsolationMVA3oldDMwoLTraw_2=DEF;
    byIsolationMVA3newDMwLTraw_2=DEF;
    byIsolationMVA3oldDMwLTraw_2=DEF;
    byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byLooseIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byMediumIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byTightIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    byVTightIsolationMVArun2017v2DBoldDMwLT2017_2=DEF;
    // byVLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byMediumIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // byVTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
    // NewMVAIDVLoose_2=DEF;
    // NewMVAIDLoose_2=DEF;
    // NewMVAIDMedium_2=DEF;
    // NewMVAIDTight_2=DEF;
    // NewMVAIDVTight_2=DEF;
    // NewMVAIDVVTight_2=DEF;
    // idMVANewDM_2=DEF;
    chargedIsoPtSum_2=DEF;
    neutralIsoPtSum_2=DEF;
    puCorrPtSum_2=DEF;
    decayModeFindingOldDMs_2=DEF;
    decayMode_2=DEF;
    //////////////////////////////////////////////////////////////////
    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {
        met[shift]=DEF;
        metphi[shift]=DEF;
        met_ex[shift]=DEF;
        met_ey[shift]=DEF;
        m_sv[shift]=DEF;
        pt_sv[shift]=DEF;
        pt_tt[shift]=DEF;
        pt_ttjj[shift]=DEF;
        m_ttjj[shift]=DEF;
        mt_1[shift]=DEF;
        mt_2[shift]=DEF;
        mt_tot[shift]=DEF;
        pt_sum[shift]=DEF;

        njets[shift]=DEF;
        njetspt20[shift]=DEF;
        njetingap[shift]=0;
        njetingap20[shift]=0;
        mjj[shift]=DEF;
        jdeta[shift]=DEF;
        dijetpt[shift]=DEF;
        dijetphi[shift]=DEF;
        jdphi[shift]=DEF;
        jpt_1[shift]=DEF;
        jpt_2[shift]=DEF;
        jeta_1[shift]=DEF;
        jeta_2[shift]=DEF;
        jphi_1[shift]=DEF;
        jphi_2[shift]=DEF;
        htxs_reco_ggf[shift]=0;
        htxs_reco_vbf[shift]=0;
    }

    jm_1=DEF;
    jm_2=DEF;
    jrawf_1=DEF;
    jrawf_2=DEF;
    jmva_1=DEF;
    jmva_2=DEF;
    jcsv_1=DEF;
    jcsv_2=DEF;

    for(int i =0 ; i < 5; ++i){
        nbtag[i]=DEF;
        bpt_1[i]=DEF;
        beta_1[i]=DEF;
        bphi_1[i]=DEF;
        brawf_1[i]=DEF;
        bmva_1[i]=DEF;
        bcsv_1[i]=DEF;
        bpt_2[i]=DEF;
        beta_2[i]=DEF;
        bphi_2[i]=DEF;
        brawf_2[i]=DEF;
        bmva_2[i]=DEF;
        bcsv_2[i]=DEF;
    }
    //////////////////////////////////////////////////////////////////
    uncorrmet=DEF;
    corrmet=DEF;
    corrmet_ex=DEF;
    corrmet_ey=DEF;
    corrmetphi=DEF;
    mvamet=DEF;
    mvamet_ex=DEF;
    mvamet_ey=DEF;
    mvametphi=DEF;
    corrmvamet_ex=DEF;
    corrmvamet_ey=DEF;
    corrmvamet=DEF;
    corrmvametphi=DEF;
    mvacov00=DEF;
    mvacov01=DEF;
    mvacov10=DEF;
    mvacov11=DEF;
    metcov00=DEF;
    metcov01=DEF;
    metcov10=DEF;
    metcov11=DEF;
    //////////////////////////////////////////////////////////////////
    passesIsoCuts=DEF;
    passesLepIsoCuts=DEF;
    passesTauLepVetos=DEF;
    passesThirdLepVeto=DEF;
    diMuonVeto=DEF;
    diElectronVeto=DEF;
    //////////////////////////////////////////////////////////////////
    dilepton_veto=DEF;
    extramuon_veto=DEF;
    extraelec_veto=DEF;
    //////////////////////////////////////////////////////////////////
    pzetavis=DEF;
    pzetamiss=DEF;
    dzeta=DEF;
    //////////////////////////////////////////////////////////////////
    pt_vis=DEF;
    m_vis=DEF;
    m_coll=DEF;
    dphi=DEF;
    //////////////////////////////////////////////////////////////////
    pfmt_1=DEF;
    pfmt_2=DEF;
    pfpt_sum=DEF;
    dr_leptau=DEF;
    jeta1eta2=DEF;
    met_centrality=DEF;
    mvamet_centrality=DEF;
    lep_etacentrality=DEF;
    sphericity=DEF;
    //////////////////////////////////////////////////////////////////

    addlepton_p4.clear();
    addlepton_pt.clear();
    addlepton_eta.clear();
    addlepton_phi.clear();
    addlepton_m.clear();    
    addlepton_iso.clear();
    addlepton_pdgId.clear();
    addlepton_mc_match.clear();
    addlepton_d0.clear();
    addlepton_dZ.clear();
    addlepton_mt.clear();
    addlepton_mvis.clear();
    addlepton_tauDM.clear();
    addlepton_tauCombIso.clear();
    addlepton_tauID.clear();
    addlepton_tauAntiEle.clear();
    addlepton_tauAntiMu.clear();

    //////////////////////////////////////////////////////////////////
    nadditionalMu = 0;
    addmuon_pt.clear();
    addmuon_eta.clear();
    addmuon_phi.clear();
    addmuon_m.clear();
    addmuon_q.clear();
    addmuon_iso.clear();
    addmuon_gen_match.clear();
    //////////////////////////////////////////////////////////////////
    nadditionalEle = 0;
    addele_pt.clear();
    addele_eta.clear();
    addele_phi.clear();
    addele_m.clear();
    addele_q.clear();
    addele_iso.clear();
    addele_gen_match.clear();
    //////////////////////////////////////////////////////////////////
    nadditionalTau = 0;
    addtau_pt.clear();
    addtau_eta.clear();
    addtau_phi.clear();
    addtau_m.clear();
    addtau_q.clear();
    addtau_byIsolationMVArun2v1DBoldDMwLTraw.clear();
    addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits.clear();
    addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
    addtau_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();
    addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
    addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byLooseIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byMediumIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byTightIsolationMVArun2v1DBoldDMwLT.clear();
    addtau_byVTightIsolationMVArun2v1DBoldDMwLT.clear();
    // addtau_byVLooseIsolationMVArun2v1DBnewDMwLT.clear();
    // addtau_byLooseIsolationMVArun2v1DBnewDMwLT.clear();
    // addtau_byMediumIsolationMVArun2v1DBnewDMwLT.clear();
    // addtau_byTightIsolationMVArun2v1DBnewDMwLT.clear();
    // addtau_byVTightIsolationMVArun2v1DBnewDMwLT.clear();
    // addtau_NewMVAIDVLoose.clear();
    // addtau_NewMVAIDLoose.clear();
    // addtau_NewMVAIDMedium.clear();
    // addtau_NewMVAIDTight.clear();
    // addtau_NewMVAIDVTight.clear();
    // addtau_NewMVAIDVVTight.clear();
    addtau_passesTauLepVetos.clear();
    addtau_decayMode.clear();
    addtau_d0.clear();
    addtau_dZ.clear();
    addtau_gen_match.clear();
    addtau_mt.clear();
    addtau_mvis.clear();
    //////////////////////////////////////////////////////////////////

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void EventWriter::initTree(TTree *t, vector< pair< string, pair<string,bool> > > jecShifts_,  bool isMC_, bool isSync_){

    isMC=isMC_;
    isSync=isSync_;

    jecShifts = jecShifts_;
    btagShifts.clear();
    btagShifts.push_back( make_pair( "",           make_pair("central","central") ) );
    btagShifts.push_back( make_pair( "MistagUp",   make_pair("up","central") ) );
    btagShifts.push_back( make_pair( "MistagDown", make_pair("down","central") ) );
    btagShifts.push_back( make_pair( "BtagUp",     make_pair("central","up") ) );
    btagShifts.push_back( make_pair( "BtagDown",   make_pair("central","down") ) );


    tauTrigSFTight = new TauTriggerSFs2017("utils/TauTriggerSFs2017/data/tauTriggerEfficiencies2017.root","tight");
    tauTrigSFVTight = new TauTriggerSFs2017("utils/TauTriggerSFs2017/data/tauTriggerEfficiencies2017.root","vtight");

    TFile wsp("utils/CorrectionWorkspaces/htt_scalefactors_v17_1.root");
    w = (RooWorkspace*)wsp.Get("w");



    if(!isSync){
        t->Branch("pfpt_sum", &pfpt_sum);
        t->Branch("dr_leptau", &dr_leptau);

        t->Branch("jeta1eta2", &jeta1eta2);
        t->Branch("met_centrality", &met_centrality);
        t->Branch("lep_etacentrality", &lep_etacentrality);
        t->Branch("sphericity", &sphericity);

        t->Branch("addlepton_p4", &addlepton_p4);
        t->Branch("addlepton_pt", &addlepton_pt);
        t->Branch("addlepton_eta", &addlepton_eta);
        t->Branch("addlepton_phi", &addlepton_phi);
        t->Branch("addlepton_m", &addlepton_m);
        t->Branch("addlepton_iso", &addlepton_iso);
        t->Branch("addlepton_pdgId", &addlepton_pdgId);
        t->Branch("addlepton_mc_match", &addlepton_mc_match);
        t->Branch("addlepton_d0", &addlepton_d0);
        t->Branch("addlepton_dZ", &addlepton_dZ);
        t->Branch("addlepton_mt", &addlepton_mt);
        t->Branch("addlepton_mvis", &addlepton_mvis);
        t->Branch("addlepton_tauCombIso", &addlepton_tauCombIso);
        t->Branch("addlepton_tauID", &addlepton_tauID);
        t->Branch("addlepton_tauDM", &addlepton_tauDM);
        t->Branch("addlepton_tauAntiEle", &addlepton_tauAntiEle);
        t->Branch("addlepton_tauAntiMu", &addlepton_tauAntiMu);

        t->Branch("nadditionalMu", &nadditionalMu);
        t->Branch("addmuon_pt", &addmuon_pt);
        t->Branch("addmuon_eta", &addmuon_eta);
        t->Branch("addmuon_phi", &addmuon_phi);
        t->Branch("addmuon_m", &addmuon_m);
        t->Branch("addmuon_q", &addmuon_q);
        t->Branch("addmuon_iso", &addmuon_iso);
        t->Branch("addmuon_gen_match", &addmuon_gen_match);

        t->Branch("nadditionalEle", &nadditionalEle);
        t->Branch("addele_pt", &addele_pt);
        t->Branch("addele_eta", &addele_eta);
        t->Branch("addele_phi", &addele_phi);
        t->Branch("addele_m", &addele_m);
        t->Branch("addele_q", &addele_q);
        t->Branch("addele_iso", &addele_iso);
        t->Branch("addele_gen_match", &addele_gen_match);

        t->Branch("nadditionalTau", &nadditionalTau);
        t->Branch("addtau_pt", &addtau_pt);
        t->Branch("addtau_eta", &addtau_eta);
        t->Branch("addtau_phi", &addtau_phi);
        t->Branch("addtau_m", &addtau_m);
        t->Branch("addtau_q", &addtau_q);
        t->Branch("addtau_byIsolationMVArun2v1DBoldDMwLTraw", &addtau_byIsolationMVArun2v1DBoldDMwLTraw);
        t->Branch("addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
        t->Branch("addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
        t->Branch("addtau_byTightCombinedIsolationDeltaBetaCorr3Hits", &addtau_byTightCombinedIsolationDeltaBetaCorr3Hits);
        t->Branch("addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
        t->Branch("addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byVVLooseIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byVLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byVLooseIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byLooseIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byMediumIsolationMVArun2v1DBoldDMwLT", &addtau_byMediumIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byTightIsolationMVArun2v1DBoldDMwLT", &addtau_byTightIsolationMVArun2v1DBoldDMwLT);
        t->Branch("addtau_byVTightIsolationMVArun2v1DBoldDMwLT", &addtau_byVTightIsolationMVArun2v1DBoldDMwLT);
        // t->Branch("addtau_byVLooseIsolationMVArun2v1DBnewDMwLT", &addtau_byVLooseIsolationMVArun2v1DBnewDMwLT);
        // t->Branch("addtau_byLooseIsolationMVArun2v1DBnewDMwLT", &addtau_byLooseIsolationMVArun2v1DBnewDMwLT);
        // t->Branch("addtau_byMediumIsolationMVArun2v1DBnewDMwLT", &addtau_byMediumIsolationMVArun2v1DBnewDMwLT);
        // t->Branch("addtau_byTightIsolationMVArun2v1DBnewDMwLT", &addtau_byTightIsolationMVArun2v1DBnewDMwLT);
        // t->Branch("addtau_byVTightIsolationMVArun2v1DBnewDMwLT", &addtau_byVTightIsolationMVArun2v1DBnewDMwLT);

        // t->Branch("addtau_NewMVAIDVLoose", &addtau_NewMVAIDVLoose);
        // t->Branch("addtau_NewMVAIDLoose", &addtau_NewMVAIDLoose);
        // t->Branch("addtau_NewMVAIDMedium", &addtau_NewMVAIDMedium);
        // t->Branch("addtau_NewMVAIDTight", &addtau_NewMVAIDTight);
        // t->Branch("addtau_NewMVAIDVTight", &addtau_NewMVAIDVTight);
        // t->Branch("addtau_NewMVAIDVVTight", &addtau_NewMVAIDVVTight);
      
        t->Branch("addtau_passesTauLepVetos", &addtau_passesTauLepVetos);
        t->Branch("addtau_decayMode", &addtau_decayMode);
        t->Branch("addtau_d0", &addtau_d0);
        t->Branch("addtau_dZ", &addtau_dZ);
        t->Branch("addtau_gen_match", &addtau_gen_match);
        t->Branch("addtau_mt", &addtau_mt);
        t->Branch("addtau_mvis", &addtau_mvis);
    }

    t->Branch("gen_Mll", &gen_Mll);
    t->Branch("genpX", &gen_ll_px);
    t->Branch("genpY", &gen_ll_py);
    t->Branch("genpZ", &gen_ll_pz);
    t->Branch("gen_top_pt_1", &gen_top_pt_1);
    t->Branch("gen_top_pt_2", &gen_top_pt_2);
    t->Branch("gen_vis_Mll", &gen_vis_Mll);
    t->Branch("vispX", &gen_ll_vis_px);
    t->Branch("vispY", &gen_ll_vis_py);
    t->Branch("vispZ", &gen_ll_vis_pz);
    t->Branch("npv", &npv);
    t->Branch("npu", &npu);
    t->Branch("rho", &rho);
    t->Branch("NUP", &NUP);
    
    t->Branch("flagMETFilter", &flagMETFilter);
    t->Branch("Flag_goodVertices", &Flag_goodVertices);
    t->Branch("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
    t->Branch("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
    t->Branch("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
    t->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
    t->Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
    t->Branch("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter);
    t->Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter);
    t->Branch("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
    t->Branch("Flag_METFilters", &Flag_METFilters);

    if(isMC){
        t->Branch("Flag_badMuons", &failBadGlobalMuonTagger);
        t->Branch("Flag_duplicateMuons", &failCloneGlobalMuonTagger);
    }
    else{
        t->Branch("Flag_badMuons", &Flag_badMuons);
        t->Branch("Flag_duplicateMuons", &Flag_duplicateMuons);
    }

    t->Branch("gen_match_1", &gen_match_1);
    t->Branch("gen_match_2", &gen_match_2);
    t->Branch("gen_match_jetId_1", &gen_match_jetId_1);
    t->Branch("gen_match_jetId_2", &gen_match_jetId_2);
    t->Branch("genJets", &genJets);
    t->Branch("genPt_1", &genPt_1);
    t->Branch("genPt_2", &genPt_2);
    t->Branch("genJet_match_1", &genJet_match_1);
    t->Branch("genJet_match_2", &genJet_match_2);
    t->Branch("pdg_1", &pdg1);
    t->Branch("pdg_2", &pdg2);

    t->Branch("pt_1", &pt_1);
    t->Branch("phi_1", &phi_1);
    t->Branch("eta_1", &eta_1);
    t->Branch("eta_SC_1", &eta_SC_1);
    t->Branch("m_1", &m_1);
    t->Branch("q_1", &q_1);
    t->Branch("d0_1", &d0_1);
    t->Branch("dZ_1", &dZ_1);
    t->Branch("pfmt_1", &pfmt_1);
    t->Branch("iso_1", &iso_1);
    t->Branch("againstElectronMVA6_1", &againstElectronMVA6_1);
    t->Branch("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1);
    t->Branch("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1);
    t->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1);
    t->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1);
    t->Branch("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1);
    t->Branch("againstMuon3_1", &againstMuon3_1);    
    t->Branch("againstMuonLoose3_1", &againstMuonLoose3_1);
    t->Branch("againstMuonTight3_1", &againstMuonTight3_1);
    t->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
    t->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_1", &byLooseCombinedIsolationDeltaBetaCorr3Hits_1);
    t->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_1", &byMediumCombinedIsolationDeltaBetaCorr3Hits_1);
    t->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_1", &byTightCombinedIsolationDeltaBetaCorr3Hits_1);
    t->Branch("byIsolationMVA3newDMwoLTraw_1", &byIsolationMVA3newDMwoLTraw_1);
    t->Branch("byIsolationMVA3oldDMwoLTraw_1", &byIsolationMVA3oldDMwoLTraw_1);
    t->Branch("byIsolationMVA3newDMwLTraw_1", &byIsolationMVA3newDMwLTraw_1);
    t->Branch("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017_1", &byIsolationMVArun2017v2DBoldDMwLTraw2017_1);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLT2017_1", &byIsolationMVArun2017v2DBoldDMwLT2017_1);   
    t->Branch("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017_1", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017_1", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byTightIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
    t->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_1);
    // t->Branch("byVLooseIsolationMVArun2v1DBnewDMwLT_1", &byVLooseIsolationMVArun2v1DBnewDMwLT_1);
    // t->Branch("byLooseIsolationMVArun2v1DBnewDMwLT_1", &byLooseIsolationMVArun2v1DBnewDMwLT_1);
    // t->Branch("byMediumIsolationMVArun2v1DBnewDMwLT_1", &byMediumIsolationMVArun2v1DBnewDMwLT_1);
    // t->Branch("byTightIsolationMVArun2v1DBnewDMwLT_1", &byTightIsolationMVArun2v1DBnewDMwLT_1);
    // t->Branch("byVTightIsolationMVArun2v1DBnewDMwLT_1", &byVTightIsolationMVArun2v1DBnewDMwLT_1);

    // t->Branch("byRerunMVAIdVLoose_1", &NewMVAIDVLoose_1);
    // t->Branch("byRerunMVAIdLoose_1", &NewMVAIDLoose_1);
    // t->Branch("byRerunMVAIdMedium_1", &NewMVAIDMedium_1);
    // t->Branch("byRerunMVAIdTight_1", &NewMVAIDTight_1);
    // t->Branch("byRerunMVAIdVTight_1", &NewMVAIDVTight_1);
    // t->Branch("byRerunMVAIdVVTight_1", &NewMVAIDVVTight_1);
    // t->Branch("idMVANewDM_1", &idMVANewDM_1);

    t->Branch("chargedIsoPtSum_1", &chargedIsoPtSum_1);
    t->Branch("neutralIsoPtSum_1", &neutralIsoPtSum_1);
    t->Branch("puCorrPtSum_1", &puCorrPtSum_1);
    t->Branch("decayModeFindingOldDMs_1", &decayModeFindingOldDMs_1);
    t->Branch("decayMode_1", &decayMode_1);
    t->Branch("id_e_mva_nt_loose_1", &id_e_mva_nt_loose_1);

    t->Branch("id_m_loose_1", &id_m_loose_1);
    t->Branch("id_m_medium_1", &id_m_medium_1);
    t->Branch("id_m_tight_1", &id_m_tight_1);
    t->Branch("id_m_tightnovtx_1", &id_m_tightnovtx_1);
    t->Branch("id_m_highpt_1", &id_m_highpt_1);
    t->Branch("id_e_cut_veto_1", &id_e_cut_veto_1);
    t->Branch("id_e_cut_loose_1", &id_e_cut_loose_1);
    t->Branch("id_e_cut_medium_1", &id_e_cut_medium_1);
    t->Branch("id_e_cut_tight_1", &id_e_cut_tight_1);

    
    t->Branch("pt_2", &pt_2);
    t->Branch("phi_2", &phi_2);
    t->Branch("eta_2", &eta_2); 
    t->Branch("m_2", &m_2);
    t->Branch("q_2", &q_2);
    t->Branch("d0_2", &d0_2);
    t->Branch("dZ_2", &dZ_2);

    t->Branch("pfmt_2", &pfmt_2);
    t->Branch("iso_2", &iso_2);
    t->Branch("againstElectronMVA6_2", &againstElectronMVA6_2);
    t->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2);
    t->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2);
    t->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2);
    t->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2);
    t->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2);
    t->Branch("againstMuon3_2", &againstMuon3_2);    
    t->Branch("againstMuonLoose3_2", &againstMuonLoose3_2);
    t->Branch("againstMuonTight3_2", &againstMuonTight3_2);
    t->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
    t->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2);
    t->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
    t->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &byTightCombinedIsolationDeltaBetaCorr3Hits_2);
    t->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2);
    t->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2);
    t->Branch("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2);
    t->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLTraw2017_2", &byIsolationMVArun2017v2DBoldDMwLTraw2017_2);
    t->Branch("byIsolationMVArun2017v2DBoldDMwLT2017_2", &byIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byVVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byVLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byLooseIsolationMVArun2017v2DBoldDMwLT2017_2", &byLooseIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byMediumIsolationMVArun2017v2DBoldDMwLT2017_2", &byMediumIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byTightIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVTightIsolationMVArun2017v2DBoldDMwLT2017_2);
    t->Branch("byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2", &byVVTightIsolationMVArun2017v2DBoldDMwLT2017_2);
    // t->Branch("byVLooseIsolationMVArun2v1DBnewDMwLT_2", &byVLooseIsolationMVArun2v1DBnewDMwLT_2);
    // t->Branch("byLooseIsolationMVArun2v1DBnewDMwLT_2", &byLooseIsolationMVArun2v1DBnewDMwLT_2);
    // t->Branch("byMediumIsolationMVArun2v1DBnewDMwLT_2", &byMediumIsolationMVArun2v1DBnewDMwLT_2);
    // t->Branch("byTightIsolationMVArun2v1DBnewDMwLT_2", &byTightIsolationMVArun2v1DBnewDMwLT_2);
    // t->Branch("byVTightIsolationMVArun2v1DBnewDMwLT_2", &byVTightIsolationMVArun2v1DBnewDMwLT_2);

    // t->Branch("byRerunMVAIdVLoose_2", &NewMVAIDVLoose_2);
    // t->Branch("byRerunMVAIdLoose_2", &NewMVAIDLoose_2);
    // t->Branch("byRerunMVAIdMedium_2", &NewMVAIDMedium_2);
    // t->Branch("byRerunMVAIdTight_2", &NewMVAIDTight_2);
    // t->Branch("byRerunMVAIdVTight_2", &NewMVAIDVTight_2);
    // t->Branch("byRerunMVAIdVVTight_2", &NewMVAIDVVTight_2);
    // t->Branch("idMVANewDM_2", &idMVANewDM_2);

    t->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2);
    t->Branch("neutralIsoPtSum_2", &neutralIsoPtSum_2);
    t->Branch("puCorrPtSum_2", &puCorrPtSum_2);
    t->Branch("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2);
    t->Branch("decayMode_2", &decayMode_2);

    t->Branch("pzetavis", &pzetavis);
    t->Branch("pzetamiss", &pzetamiss);
    t->Branch("dzeta", &dzeta);
    t->Branch("m_vis", &m_vis);
    t->Branch("m_coll", &m_coll);
    t->Branch("pt_vis", &pt_vis);
    t->Branch("dphi", &dphi);

    t->Branch("passesIsoCuts", &passesIsoCuts);
    t->Branch("passesLepIsoCuts", &passesLepIsoCuts);
    t->Branch("passesTauLepVetos", &passesTauLepVetos);
    t->Branch("passesThirdLepVeto", &passesThirdLepVeto);
    t->Branch("passesDiMuonVeto", &passesDiMuonVeto);
    t->Branch("passesDiElectronVeto", &passesDiElectronVeto);
    t->Branch("diMuonVeto", &diMuonVeto);
    t->Branch("diElectronVeto", &diElectronVeto);    
    t->Branch("dilepton_veto", &dilepton_veto);
    t->Branch("extraelec_veto", &extraelec_veto);
    t->Branch("extramuon_veto", &extramuon_veto);

    t->Branch("uncorrmet", &uncorrmet );
    t->Branch("corrmet", &corrmet);
    t->Branch("corrmetphi", &corrmetphi);
    t->Branch("corrmet_ex", &corrmet_ex);
    t->Branch("corrmet_ey", &corrmet_ey);
    t->Branch("mvamet", &mvamet);
    t->Branch("mvametphi", &mvametphi);
    t->Branch("mvamet_ex", &mvamet_ex);
    t->Branch("mvamet_ey", &mvamet_ey);
    t->Branch("corrmvamet", &corrmvamet);
    t->Branch("corrmvametphi", &corrmvametphi);
    t->Branch("corrmvamet_ex", &corrmvamet_ex);
    t->Branch("corrmvamet_ey", &corrmvamet_ey);
    t->Branch("mvacov00", &mvacov00);
    t->Branch("mvacov01", &mvacov01);
    t->Branch("mvacov10", &mvacov10);
    t->Branch("mvacov11", &mvacov11);
    t->Branch("metcov00", &metcov00);
    t->Branch("metcov01", &metcov01);
    t->Branch("metcov10", &metcov10);
    t->Branch("metcov11", &metcov11);


    for(unsigned int shift = 0; shift<jecShifts.size(); ++shift )
    {

        t->Branch( ("met"+jecShifts[shift].first).c_str(),      &met[shift]);
        t->Branch( ("metphi"+jecShifts[shift].first).c_str(),   &metphi[shift]);
        t->Branch( ("met_ex"+jecShifts[shift].first).c_str(),   &met_ex[shift]);
        t->Branch( ("met_ey"+jecShifts[shift].first).c_str(),   &met_ey[shift]);
        t->Branch( ("m_sv"+jecShifts[shift].first).c_str(),     &m_sv[shift]);
        t->Branch( ("pt_sv"+jecShifts[shift].first).c_str(),    &pt_sv[shift]);
        t->Branch( ("pt_tt"+jecShifts[shift].first).c_str(),    &pt_tt[shift]);
        t->Branch( ("pt_ttjj"+jecShifts[shift].first).c_str(),  &pt_ttjj[shift]);
        t->Branch( ("m_ttjj"+jecShifts[shift].first).c_str(),   &m_ttjj[shift]);
        t->Branch( ("pt_sum"+jecShifts[shift].first).c_str(),   &pt_sum[shift]);
        t->Branch( ("mt_1"+jecShifts[shift].first).c_str(),     &mt_1[shift]);
        t->Branch( ("mt_2"+jecShifts[shift].first).c_str(),     &mt_2[shift]);
        t->Branch( ("mt_tot"+jecShifts[shift].first).c_str(),   &mt_tot[shift]);

        t->Branch( ("njets"+jecShifts[shift].first).c_str() ,      &njets[shift]);

        t->Branch( ("njetspt20"+jecShifts[shift].first).c_str(),   &njetspt20[shift]);
        t->Branch( ("njetingap"+jecShifts[shift].first).c_str(),   &njetingap[shift]);
        t->Branch( ("njetingap20"+jecShifts[shift].first).c_str(), &njetingap20[shift]);
        t->Branch( ("dijetpt"+jecShifts[shift].first).c_str(),     &dijetpt[shift]);
        t->Branch( ("dijetphi"+jecShifts[shift].first).c_str(),    &dijetphi[shift]);
        t->Branch( ("jdphi"+jecShifts[shift].first).c_str(),       &jdphi[shift]);
        t->Branch( ("jdeta"+jecShifts[shift].first).c_str(),       &jdeta[shift]);
        t->Branch( ("mjj"+jecShifts[shift].first).c_str(),         &mjj[shift]);

        t->Branch( ("jpt_1"+jecShifts[shift].first).c_str(),       &jpt_1[shift]);
        t->Branch( ("jpt_2"+jecShifts[shift].first).c_str(),       &jpt_2[shift]);
        t->Branch( ("jeta_1"+jecShifts[shift].first).c_str(),      &jeta_1[shift]);
        t->Branch( ("jeta_2"+jecShifts[shift].first).c_str(),      &jeta_2[shift]);
        t->Branch( ("jphi_1"+jecShifts[shift].first).c_str(),      &jphi_1[shift]);
        t->Branch( ("jphi_2"+jecShifts[shift].first).c_str(),      &jphi_2[shift]);

        t->Branch( ("htxs_reco_ggf"+jecShifts[shift].first).c_str(), &htxs_reco_ggf[shift]);
        t->Branch( ("htxs_reco_vbf"+jecShifts[shift].first).c_str(), &htxs_reco_vbf[shift]);
    }
    
    t->Branch("jm_1", &jm_1);
    t->Branch("jm_2", &jm_2);
    t->Branch("jrawf_1", &jrawf_1);
    t->Branch("jrawf_2", &jrawf_2);    
    t->Branch("jmva_1", &jmva_1);
    t->Branch("jmva_2",&jmva_2);    
    t->Branch("jcsv_1", &jcsv_1);
    t->Branch("jcsv_2",&jcsv_2);    

    for(unsigned int shift=0; shift < btagShifts.size(); ++shift )
    {
        t->Branch( ("nbtag"+btagShifts[shift].first).c_str(),   &nbtag[shift]);
        t->Branch( ("bpt_1"+btagShifts[shift].first).c_str(),   &bpt_1[shift]);
        t->Branch( ("bpt_2"+btagShifts[shift].first).c_str(),   &bpt_2[shift]);
        t->Branch( ("beta_1"+btagShifts[shift].first).c_str(),  &beta_1[shift]);
        t->Branch( ("beta_2"+btagShifts[shift].first).c_str(),  &beta_2[shift]);
        t->Branch( ("bphi_1"+btagShifts[shift].first).c_str(),  &bphi_1[shift]);
        t->Branch( ("bphi_2"+btagShifts[shift].first).c_str(),  &bphi_2[shift]);
        t->Branch( ("brawf_1"+btagShifts[shift].first).c_str(), &brawf_1[shift]);
        t->Branch( ("brawf_2"+btagShifts[shift].first).c_str(), &brawf_2[shift]);
        t->Branch( ("bmva_1"+btagShifts[shift].first).c_str(),  &bmva_1[shift]);
        t->Branch( ("bmva_2"+btagShifts[shift].first).c_str(),  &bmva_2[shift]);
        t->Branch( ("bcsv_1"+btagShifts[shift].first).c_str(),  &bcsv_1[shift]);
        t->Branch( ("bcsv_2"+btagShifts[shift].first).c_str(),  &bcsv_2[shift]);
    }

    t->Branch("weight", &weight);
    t->Branch("eventWeight", &weight);
    t->Branch("lumiWeight", &lumiWeight);
    t->Branch("puweight", &puWeight);
    t->Branch("genweight", &genWeight);
    t->Branch("xsec", &xsec);
    t->Branch("genNEventsWeight", &genNEventsWeight);

    t->Branch("singleTriggerSFLeg1",&singleTriggerSFLeg1);
    t->Branch("singleTriggerSFLeg2",&singleTriggerSFLeg2);
    t->Branch("xTriggerSFLeg1",&xTriggerSFLeg1);
    t->Branch("xTriggerSFLeg2",&xTriggerSFLeg2);

    t->Branch("isoWeight_1",&isoWeight_1);
    t->Branch("isoWeight_2",&isoWeight_2);
    t->Branch("idWeight_1",&idWeight_1);
    t->Branch("idWeight_2",&idWeight_2);

    t->Branch("trigweight_1", &trigweight_1 );
    t->Branch("trigweight_2", &trigweight_2 );
    t->Branch("anti_trigweight_1", &anti_trigweight_1);
    t->Branch("idisoweight_1", &idisoweight_1);
    t->Branch("anti_idisoweight_1", &anti_idisoweight_1);
    t->Branch("idisoweight_2", &idisoweight_2);
    t->Branch("sf_trk", &sf_trk);
    t->Branch("sf_reco", &sf_reco);
    t->Branch("effweight", &effweight);

    t->Branch("sf_SingleOrCrossTrigger", &sf_SingleOrCrossTrigger);
    t->Branch("sf_SingleXorCrossTrigger", &sf_SingleXorCrossTrigger);
    t->Branch("sf_SingleTrigger", &sf_SingleTrigger);
    t->Branch("sf_DoubleTauTight", &sf_DoubleTauTight);
    t->Branch("sf_DoubleTauVTight", &sf_DoubleTauVTight);

    t->Branch("stitchedWeight", &stitchedWeight);
    t->Branch("topPtReweightWeightRun2", &topPtReweightWeightRun2);
    t->Branch("topPtReweightWeightRun1", &topPtReweightWeightRun1);
    t->Branch("zPtReweightWeight", &zPtReweightWeight);
    t->Branch("eleTauFakeRateWeight", &eleTauFakeRateWeight);
    t->Branch("muTauFakeRateWeight", &muTauFakeRateWeight);
    t->Branch("antilep_tauscaling", &antilep_tauscaling);

    t->Branch("zpt_weight_nom",&zpt_weight_nom);
    t->Branch("zpt_weight_esup",&zpt_weight_esup);
    t->Branch("zpt_weight_esdown",&zpt_weight_esdown);
    t->Branch("zpt_weight_ttup",&zpt_weight_ttup);
    t->Branch("zpt_weight_ttdown",&zpt_weight_ttdown);
    t->Branch("zpt_weight_statpt0up",&zpt_weight_statpt0up);
    t->Branch("zpt_weight_statpt0down",&zpt_weight_statpt0down);
    t->Branch("zpt_weight_statpt40up",&zpt_weight_statpt40up);
    t->Branch("zpt_weight_statpt40down",&zpt_weight_statpt40down);
    t->Branch("zpt_weight_statpt80up",&zpt_weight_statpt80up);
    t->Branch("zpt_weight_statpt80down",&zpt_weight_statpt80down);

    t->Branch("trg_singletau_leading", &trg_singletau_leading);
    t->Branch("trg_singletau_trailing", &trg_singletau_trailing);
    t->Branch("trg_singlemuon_27", &trg_singlemuon_27);
    t->Branch("trg_singlemuon_24", &trg_singlemuon_24);
    t->Branch("trg_crossmuon_mu20tau27", &trg_crossmuon_mu20tau27);
    t->Branch("trg_singleelectron_35", &trg_singleelectron_35);
    t->Branch("trg_singleelectron_32", &trg_singleelectron_32);
    t->Branch("trg_singleelectron_27", &trg_singleelectron_27);
    t->Branch("trg_crossele_ele24tau30", &trg_crossele_ele24tau30);
    t->Branch("trg_doubletau_40_tightiso", &trg_doubletau_40_tightiso);
    t->Branch("trg_doubletau_40_mediso_tightid", &trg_doubletau_40_mediso_tightid);
    t->Branch("trg_doubletau_35_tightiso_tightid", &trg_doubletau_35_tightiso_tightid);

    t->Branch("fileEntry", &fileEntry);
    t->Branch("entry", &entry);
    t->Branch("run", &run_syncro);
    t->Branch("lumi", &lumi_syncro);
    t->Branch("evt", &evt_syncro);


}
