#include "syncDATA.h"

const float DEF = -10.;

  const int gen_el_map[24]={ 6, 1,6,6,6,6, 6,6,6,6,6, 6,6,6,6,3, 6,6,6,6,6, 6,6,6 }; 
  const int gen_mu_map[24]={ 6, 2,6,6,6,6, 6,6,6,6,6, 6,6,6,6,4, 6,6,6,6,6, 6,6,6 }; 

void syncDATA::fill(HTTEvent *ev, std::vector<HTTParticle> jets, HTTPair *pair){

  lumiWeight=DEF;
  run_syncro=ev->getRunId();
  lumi_syncro=ev->getLSId();
  evt_syncro=ev->getEventId();
  //  entry=DEF; //filled in main loop
  //  fileEntry=DEF; //filled in main loop
  matchXTrig_obj=DEF;

  pdg1=std::abs(pair->getLeg1().getProperty(PropertyEnum::pdgId));
  pdg2=std::abs(pair->getLeg2().getProperty(PropertyEnum::pdgId));
  unsigned genFlav1=pair->getLeg1().getProperty(PropertyEnum::genPartFlav);
  unsigned genFlav2=pair->getLeg2().getProperty(PropertyEnum::genPartFlav);

  npv=ev->getNPV();
  npvGood=DEF;
  npu=ev->getNPU();
  rho=ev->getRho();

  gen_match_1=pair->getLeg1().getProperty(PropertyEnum::mc_match);
  gen_match_2=pair->getLeg2().getProperty(PropertyEnum::mc_match);

  /*
  //  if (pdg1==15)      gen_match_1=genFlav1;
  if (pdg1==15)      gen_match_1=pair->getLeg1().getProperty(PropertyEnum::mc_match);
  else if (pdg1==13) gen_match_1=gen_mu_map[genFlav1];
  else if (pdg1==11) gen_match_1=gen_el_map[genFlav1];

  //  if (pdg2==15)      gen_match_2=genFlav2;
  if (pdg2==15)      gen_match_2=pair->getLeg2().getProperty(PropertyEnum::mc_match);
  else if (pdg2==13) gen_match_2=gen_mu_map[genFlav2];
  else if (pdg2==11) gen_match_2=gen_el_map[genFlav2];
  */

  genPt_1=DEF;
  genPt_2=DEF;

  if (isMC){
    trk_sf = 1.;
    trigweight_1 = 1;
    idisoweight_1 = 1.;
    anti_idisoweight_1 = 1.;
    trigweight_2 = 1;
    idisoweight_2 = 1.;
    effweight = 1.;
    puWeight = 1.;
    weight = 1.;
    NUP=ev->getLHEnOutPartons();
  } else{
    trk_sf = 1.;
    trigweight_1 = 1;
    idisoweight_1 = 1.;
    anti_idisoweight_1 = 1.;
    trigweight_2 = 1;
    idisoweight_2 = 1.;
    effweight = 1.;
    puWeight = 1.;
    weight = 1.;
    NUP=-1;
  }

  stitchedWeight=1.;
  topWeight=ev->getPtReWeight();
  topWeight_run1=ev->getPtReWeightR1();
  ZWeight=ev->getPtReWeightSUSY();
  zpt_weight_nom=DEF;
  zpt_weight_esup=DEF;
  zpt_weight_esdown=DEF;
  zpt_weight_ttup=DEF;
  zpt_weight_ttdown=DEF;
  zpt_weight_statpt0up=DEF;
  zpt_weight_statpt0down=DEF;
  zpt_weight_statpt40up=DEF;
  zpt_weight_statpt40down=DEF;
  zpt_weight_statpt80up=DEF;
  zpt_weight_statpt80down=DEF;

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

  gen_top_pt_1=DEF;
  gen_top_pt_2=DEF;
  genJets=DEF;
  matchedJetPt03_1=DEF;
  matchedJetPt05_1=DEF;
  matchedJetPt03_2=DEF;
  matchedJetPt05_2=DEF;
  //////////////////////////////////////////////////////////////////  

  //this is quite slow, calling the function for each trigger item...
  if ( pdg1==13 && pdg2==15 ){ //mu-tau
    trg_singlemuon=  
      pair->getLeg1().hasTriggerMatch(TriggerEnum::HLT_IsoMu24);
    trg_mutaucross=
      ( pair->getLeg1().hasTriggerMatch(TriggerEnum::HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1) && pair->getLeg2().hasTriggerMatch(TriggerEnum::HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1) );
  } else if ( pdg1==11 && pdg2==15 ){ //e-tau
    trg_singleelectron=
      pair->getLeg1().hasTriggerMatch(TriggerEnum::HLT_Ele32_WPTight_Gsf);
  } else if ( pdg1==15 && pdg2==15 ){ //tau-tau
    trg_singletau=
      false;
    trg_doubletau= ( pair->getLeg1().hasTriggerMatch(TriggerEnum::HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg) && pair->getLeg2().hasTriggerMatch(TriggerEnum::HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg) )
                   || ( pair->getLeg1().hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg) && pair->getLeg2().hasTriggerMatch(TriggerEnum::HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg) );
  }

  trg_muonelectron=DEF; //fires HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL

  passBadMuonFilter=ev->getFilter(FilterEnum::Flag_muonBadTrackFilter); //?
  passBadChargedHadronFilter=ev->getFilter(FilterEnum::Flag_chargedHadronTrackResolutionFilter); //?

  Flag_HBHENoiseFilter=ev->getFilter(FilterEnum::Flag_HBHENoiseFilter);
  Flag_HBHENoiseIsoFilter=ev->getFilter(FilterEnum::Flag_HBHENoiseIsoFilter);
  Flag_EcalDeadCellTriggerPrimitiveFilter=ev->getFilter(FilterEnum::Flag_EcalDeadCellTriggerPrimitiveFilter);
  Flag_goodVertices=ev->getFilter(FilterEnum::Flag_goodVertices);
  Flag_eeBadScFilter=ev->getFilter(FilterEnum::Flag_eeBadScFilter);
  Flag_globalTightHalo2016Filter=ev->getFilter(FilterEnum::Flag_globalTightHalo2016Filter);

  Flag_badMuons=ev->getFilter(FilterEnum::Flag_muonBadTrackFilter); //??

  failBadGlobalMuonTagger=DEF;
  failCloneGlobalMuonTagger=DEF;
  Flag_duplicateMuons=DEF;
  Flag_noBadMuons=DEF;

  passesMetMuonFilter=passBadMuonFilter && passBadChargedHadronFilter && Flag_HBHENoiseFilter && Flag_HBHENoiseIsoFilter && Flag_EcalDeadCellTriggerPrimitiveFilter && Flag_goodVertices && Flag_eeBadScFilter && Flag_globalTightHalo2016Filter;

  //////////////////////////////////////////////////////////////////  
  HTTParticle leg1=pair->getLeg1();
  TLorentzVector leg1P4=leg1.getP4();
  pt_1=leg1P4.Pt();
  phi_1=leg1P4.Phi();
  eta_1=leg1P4.Eta();
  eta_SC_1=eta_1+leg1.getProperty(PropertyEnum::deltaEtaSC);
  m_1=leg1P4.M();
  q_1=leg1.getCharge();
  d0_1=leg1.getProperty(PropertyEnum::dxy);
  dZ_1=leg1.getProperty(PropertyEnum::dz);
  mt_1=pair->getMTLeg1();

  if (pdg1==15)      iso_1=leg1.getProperty(PropertyEnum::rawMVAoldDM2017v2);
  else if (pdg1==13) iso_1=leg1.getProperty(PropertyEnum::pfRelIso04_all);
  else if (pdg1==11) iso_1=leg1.getProperty(PropertyEnum::pfRelIso03_all);;

  UChar_t bitmask=leg1.getProperty(PropertyEnum::idAntiEle);
  againstElectronVLooseMVA6_1 =(bitmask & 0x1 )>0;
  againstElectronLooseMVA6_1  =(bitmask & 0x2 )>0;
  againstElectronMediumMVA6_1 =(bitmask & 0x4 )>0;
  againstElectronTightMVA6_1  =(bitmask & 0x8 )>0;
  againstElectronVTightMVA6_1 =(bitmask & 0x10)>0;

  bitmask=leg1.getProperty(PropertyEnum::idAntiMu);
  againstMuonLoose3_1=(bitmask & 0x1)>0;
  againstMuonTight3_1=(bitmask & 0x2)>0;

  byCombinedIsolationDeltaBetaCorrRaw3Hits_1=leg1.getProperty(PropertyEnum::rawIso);
  byLooseCombinedIsolationDeltaBetaCorr3Hits_1=DEF; //not in nanoAOD
  byMediumCombinedIsolationDeltaBetaCorr3Hits_1=DEF; //not in nanoAOD
  byTightCombinedIsolationDeltaBetaCorr3Hits_1=DEF; //not in nanoAOD
  byIsolationMVA3newDMwoLTraw_1=DEF;
  byIsolationMVA3oldDMwoLTraw_1=leg1.getProperty(PropertyEnum::rawMVAoldDM2017v2);
  byIsolationMVA3newDMwLTraw_1=DEF;
  byIsolationMVA3oldDMwLTraw_1 =leg1.getProperty(PropertyEnum::rawMVAoldDM2017v2); //same as above!?

  bitmask=leg1.getProperty(PropertyEnum::idMVAoldDM2017v2);
  byVLooseIsolationMVArun2v1DBoldDMwLT_1=(bitmask & 0x1)>0;
  byLooseIsolationMVArun2v1DBoldDMwLT_1=(bitmask & 0x2)>0;;
  byMediumIsolationMVArun2v1DBoldDMwLT_1=(bitmask & 0x4)>0;
  byTightIsolationMVArun2v1DBoldDMwLT_1=(bitmask & 0x8)>0;
  byVTightIsolationMVArun2v1DBoldDMwLT_1=(bitmask & 0x10)>0;

  

  byVLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byMediumIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byTightIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byVTightIsolationMVArun2v1DBnewDMwLT_1=DEF;

  NewMVAIDVLoose_1=DEF; //?
  NewMVAIDLoose_1=DEF;
  NewMVAIDMedium_1=DEF;
  NewMVAIDTight_1=DEF;
  NewMVAIDVTight_1=DEF;
  NewMVAIDVVTight_1=DEF;

  idMVANewDM_1=DEF; //?

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

  if (pdg1==15) gen_match_jetId_1=getGenMatch_jetId(leg1P4,jets);
  
  //////////////////////////////////////////////////////////////////
  HTTParticle leg2=pair->getLeg2();
  TLorentzVector leg2P4=leg2.getP4();
  pt_2=leg2P4.Pt();
  phi_2=leg2P4.Phi();
  eta_2=leg2P4.Eta();
  m_2=leg2P4.M();
  q_2=leg2.getCharge();
  d0_2=leg2.getProperty(PropertyEnum::dxy);
  dZ_2=leg2.getProperty(PropertyEnum::dz);
  mt_2=pair->getMTLeg2();

  if (pdg2==15)      iso_2=leg2.getProperty(PropertyEnum::rawMVAoldDM2017v2);
  else if (pdg2==13) iso_2=leg2.getProperty(PropertyEnum::pfRelIso04_all);
  else if (pdg2==11) iso_2=leg2.getProperty(PropertyEnum::pfRelIso03_all);;

  bitmask=leg2.getProperty(PropertyEnum::idAntiEle);
  againstElectronVLooseMVA6_2=(bitmask & 0x1)>0;
  againstElectronLooseMVA6_2= (bitmask & 0x2)>0;
  againstElectronMediumMVA6_2=(bitmask & 0x4)>0;
  againstElectronTightMVA6_2= (bitmask & 0x8)>0;
  againstElectronVTightMVA6_2=(bitmask & 0x10)>0;

  bitmask=leg2.getProperty(PropertyEnum::idAntiMu);
  againstMuonLoose3_2=(bitmask & 0x1)>0;
  againstMuonTight3_2=(bitmask & 0x2)>0;

  byCombinedIsolationDeltaBetaCorrRaw3Hits_2=leg2.getProperty(PropertyEnum::rawIso);
  byLooseCombinedIsolationDeltaBetaCorr3Hits_2=DEF; //not in nanoAOD
  byMediumCombinedIsolationDeltaBetaCorr3Hits_2=DEF; //not in nanoAOD
  byTightCombinedIsolationDeltaBetaCorr3Hits_2=DEF; //not in nanoAOD
  byIsolationMVA3newDMwoLTraw_2=DEF;
  byIsolationMVA3oldDMwoLTraw_2=leg2.getProperty(PropertyEnum::rawMVAoldDM2017v2);
  byIsolationMVA3newDMwLTraw_2=DEF;
  byIsolationMVA3oldDMwLTraw_2 =leg2.getProperty(PropertyEnum::rawMVAoldDM2017v2); //same as above!?

  bitmask=leg2.getProperty(PropertyEnum::idMVAoldDM2017v2);
  byVLooseIsolationMVArun2v1DBoldDMwLT_2=(bitmask & 0x1)>0;
  byLooseIsolationMVArun2v1DBoldDMwLT_2=(bitmask & 0x2)>0;;
  byMediumIsolationMVArun2v1DBoldDMwLT_2=(bitmask & 0x4)>0;
  byTightIsolationMVArun2v1DBoldDMwLT_2=(bitmask & 0x8)>0;
  byVTightIsolationMVArun2v1DBoldDMwLT_2=(bitmask & 0x10)>0;

  byVLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byMediumIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byVTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
  NewMVAIDVLoose_2=DEF;
  NewMVAIDLoose_2=DEF;
  NewMVAIDMedium_2=DEF;
  NewMVAIDTight_2=DEF;
  NewMVAIDVTight_2=DEF;
  NewMVAIDVVTight_2=DEF;
  idMVANewDM_2=DEF;

  chargedIsoPtSum_2=leg2.getProperty(PropertyEnum::chargedIso);
  neutralIsoPtSum_2=leg2.getProperty(PropertyEnum::neutralIso);
  puCorrPtSum_2=leg2.getProperty(PropertyEnum::puCorr);
  decayModeFindingOldDMs_2=leg2.getProperty(PropertyEnum::idDecayMode);
  decayMode_2=leg2.getProperty(PropertyEnum::decayMode);

  if (pdg2==15) gen_match_jetId_2=getGenMatch_jetId(leg2P4,jets);
  //////////////////////////////////////////////////////////////////
  nbtag=0;
  njets=0; 
  int ind_b1=-1;
  int ind_b2=-1;
  for (unsigned ij=0; ij<jets.size(); ij++){ 
    if (jets.at(ij).getP4().Pt()>30) njets++; 
    if ( std::abs(jets.at(ij).getP4().Eta())<2.4 && jets.at(ij).getProperty(PropertyEnum::btagCSVV2)>0.8484 ){
      nbtag++; 
      if ( ind_b1>=0 && ind_b2<0 ) ind_b2=ij;
      if ( ind_b1<0 )              ind_b1=ij;
    }
    if (evt_syncro==1279980){ std::cout << ij << " " << ind_b1 << " " << jets.at(ij).getProperty(PropertyEnum::btagCSVV2) << " " << jets.at(ij).getP4().Pt() << " " << jets.at(ij).getP4().Eta()  << " " <<std::endl; }
  }
  njetsUp=njets;
  njetsDown=njets;
  njetspt20=jets.size();
  TLorentzVector j1;
  TLorentzVector j2;

  if ( jets.size()>=1 ){
    jpt_1=jets.at(0).getP4().Pt();
    jptUp_1=jpt_1;
    jptDown_1=jpt_1;
    jeta_1=jets.at(0).getP4().Eta();
    jphi_1=jets.at(0).getP4().Phi();
    jm_1=jets.at(0).getP4().M();
    jrawf_1=jets.at(0).getProperty(PropertyEnum::rawFactor);
    jmva_1=jets.at(0).getProperty(PropertyEnum::btagCMVA);
    jcsv_1=jets.at(0).getProperty(PropertyEnum::btagCSVV2);
    //    gen_match_jetId_1=jets.at(0).getProperty(PropertyEnum::partonFlavour);
    genJet_match_1=0;
    j1=jets.at(0).getP4();
  }
  if ( jets.size()>=2 ){
    jpt_2=jets.at(1).getP4().Pt();
    jptUp_2=jpt_2;
    jptDown_2=jpt_2;
    jeta_2=jets.at(1).getP4().Eta();
    jphi_2=jets.at(1).getP4().Phi();
    jm_2=jets.at(1).getP4().M();
    jrawf_2=jets.at(1).getProperty(PropertyEnum::rawFactor);
    jmva_2=jets.at(1).getProperty(PropertyEnum::btagCMVA);
    jcsv_2=jets.at(1).getProperty(PropertyEnum::btagCSVV2);
    //    gen_match_jetId_2=jets.at(1).getProperty(PropertyEnum::partonFlavour);
    genJet_match_2=0;
    jeta1eta2=jeta_1*jeta_2;
    lep_etacentrality=TMath::Exp( -4/pow(jeta_1-jeta_2,2) * pow( (eta_1-( jeta_1+jeta_2 )*0.5), 2 ) );

    j2=jets.at(1).getP4();
    TLorentzVector jj=j1+j2;

    mjj=jj.M();
    dijetpt=jj.Pt();
    dijetphi=jj.Phi();
    
    jdeta=std::abs( j1.Eta()-j2.Eta() );
    jdphi=j1.DeltaPhi(j2);

    njetingap=0;
    njetingap20=0;

    for (unsigned ij=2; ij<jets.size(); ij++){
      float aj_eta=jets.at(ij).getP4().Eta();
      //      if ( ( aj_eta<j1.Eta() && aj_eta>j2.Eta() ) || ( aj_eta>j1.Eta() && aj_eta<j2.Eta() ) ) ){

      if ( ( aj_eta<j1.Eta() && aj_eta>j2.Eta() ) || ( aj_eta>j1.Eta() && aj_eta<j2.Eta() ) ){
	  njetingap20++;
	  if ( jets.at(ij).getP4().Pt()>30 ) njetingap++;
      }
    }
  }
  mjjUp=mjj;
  mjjDown=mjj;
  jdetaUp=jdeta;
  jdetaDown=jdeta;

  if (ind_b1>=0){
    bpt_1=jets.at(ind_b1).getP4().Pt();
    beta_1=jets.at(ind_b1).getP4().Eta();
    bphi_1=jets.at(ind_b1).getP4().Phi();
    brawf_1=jets.at(ind_b1).getProperty(PropertyEnum::rawFactor);
    bmva_1=jets.at(ind_b1).getProperty(PropertyEnum::btagCMVA);
    bcsv_1=jets.at(ind_b1).getProperty(PropertyEnum::btagCSVV2);
    if (evt_syncro==1279980){ std::cout << ind_b1 << " " << ind_b1 << " " << jets.at(ind_b1).getProperty(PropertyEnum::btagCSVV2) << " " << jets.at(ind_b1).getP4().Pt() << " " << jets.at(ind_b1).getP4().Eta()  << " " << bcsv_1 << " XX " <<std::endl; }
  }

  if (ind_b2>=0){
    bpt_2=jets.at(ind_b2).getP4().Pt();
    beta_2=jets.at(ind_b2).getP4().Eta();
    bphi_2=jets.at(ind_b2).getP4().Phi();
    brawf_2=jets.at(ind_b2).getProperty(PropertyEnum::rawFactor);
    bmva_2=jets.at(ind_b2).getProperty(PropertyEnum::btagCMVA);
    bcsv_2=jets.at(ind_b2).getProperty(PropertyEnum::btagCSVV2);
  }

  //////////////////////////////////////////////////////////////////
  met   =pair->getMET().Mod();
  met_ex=pair->getMET().X();
  met_ey=pair->getMET().Y();
  metphi=pair->getMET().Phi();
  if (metphi>TMath::Pi()) metphi-=2*TMath::Pi();

  metcov00=pair->getMETMatrix().at(0);
  metcov01=pair->getMETMatrix().at(1);
  metcov10=pair->getMETMatrix().at(2);
  metcov11=pair->getMETMatrix().at(3);

  uncorrmet=ev->getMET_uncorr().Mod();

  corrmet=met;
  corrmet_ex=met_ex;
  corrmet_ey=met_ey;
  corrmetphi=metphi;

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

  m_sv=pair->getP4().M();
  pt_sv=pair->getP4().Pt();
  //////////////////////////////////////////////////////////////////
  eleTauFakeRateWeight=DEF;
  muTauFakeRateWeight=DEF;
  antilep_tauscaling=DEF;
  //////////////////////////////////////////////////////////////////
  if (pdg2==15){
    if (pdg1==15){passesIsoCuts=byVLooseIsolationMVArun2v1DBoldDMwLT_1 && byVLooseIsolationMVArun2v1DBoldDMwLT_2;} //tautau
    else{passesIsoCuts=byVLooseIsolationMVArun2v1DBoldDMwLT_2;} //etau/mutau
  }
  if (pdg1==11 || pdg1==13) passesLepIsoCuts=(iso_1<0.1);

  if( pdg1==13 && pdg2==15 ) passesTauLepVetos = againstElectronLooseMVA6_2 && againstMuonTight3_2;
  if( pdg1==11 && pdg2==15 ) passesTauLepVetos = againstElectronTightMVA6_2 && againstMuonLoose3_2;
  if( pdg1==15 && pdg2==15 ) passesTauLepVetos = againstElectronVLooseMVA6_1 && againstMuonLoose3_1 && againstElectronVLooseMVA6_2 && againstMuonLoose3_2;

  passesThirdLepVeto=!( ev->checkSelectionBit(SelectionBitsEnum::extraMuonVeto) && ev->checkSelectionBit(SelectionBitsEnum::extraElectronVeto) );
  passesDiMuonVeto=!( ev->checkSelectionBit(SelectionBitsEnum::diMuonVeto) );
  passesDiElectronVeto=!( ev->checkSelectionBit(SelectionBitsEnum::diElectronVeto) );
  //////////////////////////////////////////////////////////////////
  dilepton_veto=!(passesDiMuonVeto && passesDiElectronVeto);
  extramuon_veto=ev->checkSelectionBit(SelectionBitsEnum::extraMuonVeto);
  extraelec_veto=ev->checkSelectionBit(SelectionBitsEnum::extraElectronVeto);
  //////////////////////////////////////////////////////////////////
  double zetaX = TMath::Cos(phi_1) + TMath::Cos(phi_2);
  double zetaY = TMath::Sin(phi_1) + TMath::Sin(phi_2);
  double zetaR = TMath::Sqrt(zetaX*zetaX + zetaY*zetaY);
  if ( zetaR > 0. ) {
    zetaX /= zetaR;
    zetaY /= zetaR;
  }
  pzetavis =  (leg1P4.Px() + leg2P4.Px())*zetaX + (leg1P4.Py() + leg2P4.Py())*zetaY;
  pzetamiss = met_ex*zetaX + met_ey*zetaY;
  dzeta = this->pzetamiss + (pzetavis - 1.85 * pzetavis);
  //////////////////////////////////////////////////////////////////
  pt_tt=pair->getPTTOT();
  pt_vis=pair->getPTVis();
  mt_tot=pair->getMTTOT();
  m_vis=pair->getMVis();
  pfpt_tt=pt_tt;
  dphi=leg1P4.DeltaPhi(leg2P4);

  double x1 = pt_1 / ( pt_1 + met );
  double x2 = pt_2 / ( pt_2 + met );
  if( TMath::Cos( dphi ) > -0.95
      && ( x1 > 0 && x1 < 1)
      && ( x2 > 0 && x2 < 1)
      ){
    m_coll =  m_vis / sqrt( x1 * x2 ) ;
  }
  else m_coll = -999;

  //////////////////////////////////////////////////////////////////
  pfmt_1=mt_1;
  pfmt_2=mt_2;
  pt_sum=pt_1+pt_2+met;
  pfpt_sum=pt_sum;
  dr_leptau=leg1P4.DeltaR(leg2P4);
  mvamet_centrality=DEF;

  TLorentzVector v0=leg1P4*( 1/sqrt(  pow(leg1P4.Px(),2)+pow(leg1P4.Py(),2)  ) ); //lep, normalized in transverse plane
  TLorentzVector v1=leg2P4*( 1/sqrt(  pow(leg2P4.Px(),2)+pow(leg2P4.Py(),2)  ) ); //tau, normalized in transverse plane
  float omega=v1.DeltaPhi(v0);
  float theta=-v0.Phi();
  float x=(     met_ex * TMath::Sin(omega-theta)  - met_ey*TMath::Cos(omega-theta)   ) / TMath::Sin(omega); //x coord in lep-tau system
  float y=(     met_ex * TMath::Sin(theta)        + met_ey*TMath::Cos(theta)         ) / TMath::Sin(omega); //y coord in lep-tau system
  met_centrality=( x+y ) / sqrt(x*x + y*y);

  vector<TLorentzVector> objs;
  objs.push_back(leg1P4);
  objs.push_back(leg2P4);
  if ( njetspt20>0 ) objs.push_back(j1);
  if ( njetspt20>1 ) objs.push_back(j2);
  sphericity=calcSphericity(objs);

}

double syncDATA::calcSphericity(std::vector<TLorentzVector> p){

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

double syncDATA::calcSphericityFromMatrix(TMatrixD M) {

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

int syncDATA::getGenMatch_jetId(TLorentzVector selObj, std::vector<HTTParticle> jets){
  float minDR=1;
  int whichjet=0;

  for (unsigned i=0; i<jets.size(); i++){
    TLorentzVector p4=jets.at(i).getP4();
    if(p4.Pt() > 20 && fabs(p4.Eta() ) < 4.7 ){
      float tmpDR = calcDR( selObj.Eta(), selObj.Phi(), p4.Eta(), p4.Phi() );
      if( tmpDR < minDR ){
	minDR = tmpDR;
	whichjet=i;
      }
    }
  }

  if( minDR < 0.5 ) return jets.at(whichjet).getProperty(PropertyEnum::partonFlavour);
  return -99;
}

double syncDATA::calcDR(double eta1, double phi1, double eta2, double phi2){
  double deta = eta1-eta2;
  double dphi = TVector2::Phi_mpi_pi(phi1-phi2);
  return TMath::Sqrt( deta*deta+dphi*dphi );
}

void syncDATA::setDefault(){

  lumiWeight=DEF;
  run_syncro=DEF;
  lumi_syncro=DEF;
  evt_syncro=0; //unsigned
  entry=DEF;
  fileEntry=DEF;
  matchXTrig_obj=DEF;

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
  weight=DEF;
  puWeight=DEF;
  genWeight=DEF;
  trigweight_1=DEF;
  anti_trigweight_1=DEF;
  trigweight_2=DEF;
  idisoweight_1=DEF;
  anti_idisoweight_1=DEF;
  idisoweight_2=DEF;
  trk_sf=DEF;
  effweight=DEF;
  stitchedWeight=DEF;
  topWeight=DEF;
  topWeight_run1=DEF;
  ZWeight=DEF;
  zpt_weight_nom=DEF;
  zpt_weight_esup=DEF;
  zpt_weight_esdown=DEF;
  zpt_weight_ttup=DEF;
  zpt_weight_ttdown=DEF;
  zpt_weight_statpt0up=DEF;
  zpt_weight_statpt0down=DEF;
  zpt_weight_statpt40up=DEF;
  zpt_weight_statpt40down=DEF;
  zpt_weight_statpt80up=DEF;
  zpt_weight_statpt80down=DEF;
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
  matchedJetPt03_1=DEF;
  matchedJetPt05_1=DEF;
  matchedJetPt03_2=DEF;
  matchedJetPt05_2=DEF;
  //////////////////////////////////////////////////////////////////  
  trg_singlemuon=DEF; //fires OR of HLT_IsoMu22, HLT_IsoTkMu22, HLT_IsoMu22eta2p1, HLT_IsoTkMu22_eta2p1
  trg_mutaucross=DEF;
  trg_singleelectron=DEF; //fires HLT_Ele25_eta2p1_WPTight_Gsf
  trg_singletau=DEF; //fires HLT_VLooseIsoPFTau120_Trk50_eta2p1
  trg_doubletau=DEF; //fires HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg or HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg
  trg_muonelectron=DEF; //fires HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
  passBadMuonFilter=DEF;
  passBadChargedHadronFilter=DEF;
  Flag_HBHENoiseFilter=DEF;
  Flag_HBHENoiseIsoFilter=DEF;
  Flag_EcalDeadCellTriggerPrimitiveFilter=DEF;
  Flag_goodVertices=DEF;
  Flag_eeBadScFilter=DEF;
  Flag_globalTightHalo2016Filter=DEF;
  failBadGlobalMuonTagger=DEF;
  failCloneGlobalMuonTagger=DEF;
  Flag_duplicateMuons=DEF;
  Flag_badMuons=DEF;
  Flag_noBadMuons=DEF;
  passesMetMuonFilter=DEF;
  //////////////////////////////////////////////////////////////////  
  pt_1=DEF;
  phi_1=DEF;
  eta_1=DEF;
  eta_SC_1=DEF;
  m_1=DEF;
  q_1=DEF;
  d0_1=DEF;
  dZ_1=DEF;
  mt_1=DEF;
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
  byVLooseIsolationMVArun2v1DBoldDMwLT_1=DEF;
  byLooseIsolationMVArun2v1DBoldDMwLT_1=DEF;
  byMediumIsolationMVArun2v1DBoldDMwLT_1=DEF;
  byTightIsolationMVArun2v1DBoldDMwLT_1=DEF;
  byVTightIsolationMVArun2v1DBoldDMwLT_1=DEF;
  byVLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byLooseIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byMediumIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byTightIsolationMVArun2v1DBnewDMwLT_1=DEF;
  byVTightIsolationMVArun2v1DBnewDMwLT_1=DEF;
  NewMVAIDVLoose_1=DEF;
  NewMVAIDLoose_1=DEF;
  NewMVAIDMedium_1=DEF;
  NewMVAIDTight_1=DEF;
  NewMVAIDVTight_1=DEF;
  NewMVAIDVVTight_1=DEF;
  idMVANewDM_1=DEF;
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
  mt_2=DEF;
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
  byVLooseIsolationMVArun2v1DBoldDMwLT_2=DEF;
  byLooseIsolationMVArun2v1DBoldDMwLT_2=DEF;
  byMediumIsolationMVArun2v1DBoldDMwLT_2=DEF;
  byTightIsolationMVArun2v1DBoldDMwLT_2=DEF;
  byVTightIsolationMVArun2v1DBoldDMwLT_2=DEF;
  byVLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byLooseIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byMediumIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
  byVTightIsolationMVArun2v1DBnewDMwLT_2=DEF;
  NewMVAIDVLoose_2=DEF;
  NewMVAIDLoose_2=DEF;
  NewMVAIDMedium_2=DEF;
  NewMVAIDTight_2=DEF;
  NewMVAIDVTight_2=DEF;
  NewMVAIDVVTight_2=DEF;
  idMVANewDM_2=DEF;
  chargedIsoPtSum_2=DEF;
  neutralIsoPtSum_2=DEF;
  puCorrPtSum_2=DEF;
  decayModeFindingOldDMs_2=DEF;
  decayMode_2=DEF;
  //////////////////////////////////////////////////////////////////
  nbtag=DEF;
  njets=DEF;
  njetsUp=DEF;
  njetsDown=DEF;
  njetspt20=DEF;
  mjj=DEF;
  mjjUp=DEF;
  mjjDown=DEF;
  jdeta=DEF;
  jdetaUp=DEF;
  jdetaDown=DEF;
  njetingap=0;
  njetingap20=0;
  dijetpt=DEF;
  dijetphi=DEF;
  jdphi=DEF;
  jpt_1=DEF;
  jptUp_1=DEF;
  jptDown_1=DEF;
  jeta_1=DEF;
  jphi_1=DEF;
  jm_1=DEF;
  jrawf_1=DEF;
  jmva_1=DEF;
  jcsv_1=DEF;
  jpt_2=DEF;
  jptUp_2=DEF;
  jptDown_2=DEF;
  jeta_2=DEF;
  jphi_2=DEF;
  jm_2=DEF;
  jrawf_2=DEF;
  jmva_2=DEF;
  jcsv_2=DEF;
  bpt_1=DEF;
  beta_1=DEF;
  bphi_1=DEF;
  brawf_1=DEF;
  bmva_1=DEF;
  bcsv_1=DEF;
  bpt_2=DEF;
  beta_2=DEF;
  bphi_2=DEF;
  brawf_2=DEF;
  bmva_2=DEF;
  bcsv_2=DEF;
  //////////////////////////////////////////////////////////////////
  met=DEF;
  uncorrmet=DEF;
  met_ex=DEF;
  met_ey=DEF;
  metphi=DEF;
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
  m_sv=DEF;
  pt_sv=DEF;
  //////////////////////////////////////////////////////////////////
  eleTauFakeRateWeight=DEF;
  muTauFakeRateWeight=DEF;
  antilep_tauscaling=DEF;
  //////////////////////////////////////////////////////////////////
  passesIsoCuts=DEF;
  passesLepIsoCuts=DEF;
  passesTauLepVetos=DEF;
  passesThirdLepVeto=DEF;
  passesDiMuonVeto=DEF;
  passesDiElectronVeto=DEF;
  //////////////////////////////////////////////////////////////////
  dilepton_veto=DEF;
  extramuon_veto=DEF;
  extraelec_veto=DEF;
  //////////////////////////////////////////////////////////////////
  pzetavis=DEF;
  pzetamiss=DEF;
  dzeta=DEF;
  //////////////////////////////////////////////////////////////////
  pt_tt=DEF;
  pt_vis=DEF;
  mt_tot=DEF;
  pfpt_tt=DEF;
  m_vis=DEF;
  m_coll=DEF;
  dphi=DEF;
  //////////////////////////////////////////////////////////////////
  pfmt_1=DEF;
  pfmt_2=DEF;
  pfpt_sum=DEF;
  pt_sum=DEF;
  dr_leptau=DEF;
  jeta1eta2=DEF;
  met_centrality=DEF;
  mvamet_centrality=DEF;
  lep_etacentrality=DEF;
  sphericity=DEF;
}

void syncDATA::initTree(TTree *t, bool isMC_, bool isSync_){

  isMC=isMC_;
  isSync=isSync_;

  t->Branch("fileEntry", &fileEntry);
  t->Branch("entry", &entry);
  t->Branch("run", &run_syncro);
  t->Branch("lumi", &lumi_syncro);
  t->Branch("evt", &evt_syncro);
  t->Branch("weight", &weight);
  t->Branch("eventWeight", &weight);
  t->Branch("lumiWeight", &lumiWeight);
  t->Branch("puweight", &puWeight);
  t->Branch("genweight", &genWeight);
  t->Branch("trigweight_1", &trigweight_1);
  t->Branch("anti_trigweight_1", &anti_trigweight_1);
  t->Branch("trigweight_2", &trigweight_2);
  t->Branch("idisoweight_1", &idisoweight_1);
  t->Branch("anti_idisoweight_1", &anti_idisoweight_1);
  t->Branch("idisoweight_2", &idisoweight_2);
  t->Branch("trk_sf", &trk_sf);
  t->Branch("effweight", &effweight);
  t->Branch("stitchedWeight", &stitchedWeight);
  t->Branch("topWeight", &topWeight);
  t->Branch("topWeight_run1", &topWeight_run1);
  t->Branch("zPtReweightWeight", &ZWeight);

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

  t->Branch("trg_singlemuon", &trg_singlemuon);
  t->Branch("trg_mutaucross", &trg_mutaucross);
  t->Branch("trg_singleelectron", &trg_singleelectron);
  t->Branch("trg_singletau", &trg_singletau);
  t->Branch("trg_doubletau", &trg_doubletau);
  t->Branch("trg_muonelectron", &trg_muonelectron);
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
  
  t->Branch("passBadMuonFilter", &passBadMuonFilter);
  t->Branch("passBadChargedHadronFilter", &passBadChargedHadronFilter);
  t->Branch("flagHBHENoiseFilter", &Flag_HBHENoiseFilter);
  t->Branch("flagHBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
  t->Branch("flagEcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
  t->Branch("flagGoodVertices", &Flag_goodVertices);
  t->Branch("flagEeBadScFilter", &Flag_eeBadScFilter);
  t->Branch("flagGlobalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);

  if(isMC){
    t->Branch("Flag_badMuons", &failBadGlobalMuonTagger);
    t->Branch("Flag_duplicateMuons", &failCloneGlobalMuonTagger);
  }
  else{
    t->Branch("Flag_badMuons", &Flag_badMuons);
    t->Branch("Flag_duplicateMuons", &Flag_duplicateMuons);
  }

  t->Branch("passesFilter", &passesMetMuonFilter);

  t->Branch("matchedJetPt03_1", &matchedJetPt03_1);
  t->Branch("matchedJetPt05_1", &matchedJetPt05_1);
  t->Branch("matchedJetPt03_2", &matchedJetPt03_2);
  t->Branch("matchedJetPt05_2", &matchedJetPt05_2);

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
  t->Branch("mt_1", &mt_1);
  t->Branch("pfmt_1", &pfmt_1);
  t->Branch("iso_1", &iso_1);
  t->Branch("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1);
  t->Branch("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1);
  t->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1);
  t->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1);
  t->Branch("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1);
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
  t->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_1", &byVLooseIsolationMVArun2v1DBoldDMwLT_1);
  t->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_1", &byLooseIsolationMVArun2v1DBoldDMwLT_1);
  t->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_1", &byMediumIsolationMVArun2v1DBoldDMwLT_1);
  t->Branch("byTightIsolationMVArun2v1DBoldDMwLT_1", &byTightIsolationMVArun2v1DBoldDMwLT_1);
  t->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_1", &byVTightIsolationMVArun2v1DBoldDMwLT_1);
  t->Branch("byVLooseIsolationMVArun2v1DBnewDMwLT_1", &byVLooseIsolationMVArun2v1DBnewDMwLT_1);
  t->Branch("byLooseIsolationMVArun2v1DBnewDMwLT_1", &byLooseIsolationMVArun2v1DBnewDMwLT_1);
  t->Branch("byMediumIsolationMVArun2v1DBnewDMwLT_1", &byMediumIsolationMVArun2v1DBnewDMwLT_1);
  t->Branch("byTightIsolationMVArun2v1DBnewDMwLT_1", &byTightIsolationMVArun2v1DBnewDMwLT_1);
  t->Branch("byVTightIsolationMVArun2v1DBnewDMwLT_1", &byVTightIsolationMVArun2v1DBnewDMwLT_1);

  t->Branch("byRerunMVAIdVLoose_1", &NewMVAIDVLoose_1);
  t->Branch("byRerunMVAIdLoose_1", &NewMVAIDLoose_1);
  t->Branch("byRerunMVAIdMedium_1", &NewMVAIDMedium_1);
  t->Branch("byRerunMVAIdTight_1", &NewMVAIDTight_1);
  t->Branch("byRerunMVAIdVTight_1", &NewMVAIDVTight_1);
  t->Branch("byRerunMVAIdVVTight_1", &NewMVAIDVVTight_1);


  t->Branch("idMVANewDM_1", &idMVANewDM_1);
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

  t->Branch("antilep_tauscaling", &antilep_tauscaling);
  
  t->Branch("pt_2", &pt_2);
  t->Branch("phi_2", &phi_2);
  t->Branch("eta_2", &eta_2); 
  t->Branch("m_2", &m_2);
  t->Branch("q_2", &q_2);
  t->Branch("d0_2", &d0_2);
  t->Branch("dZ_2", &dZ_2);
  t->Branch("mt_2", &mt_2);
  t->Branch("pfmt_2", &pfmt_2);
  t->Branch("iso_2", &iso_2);
  t->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2);
  t->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2);
  t->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2);
  t->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2);
  t->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2);
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
  t->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &byVLooseIsolationMVArun2v1DBoldDMwLT_2);
  t->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_2", &byLooseIsolationMVArun2v1DBoldDMwLT_2);
  t->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2", &byMediumIsolationMVArun2v1DBoldDMwLT_2);
  t->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2", &byTightIsolationMVArun2v1DBoldDMwLT_2);
  t->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_2", &byVTightIsolationMVArun2v1DBoldDMwLT_2);
  t->Branch("byVLooseIsolationMVArun2v1DBnewDMwLT_2", &byVLooseIsolationMVArun2v1DBnewDMwLT_2);
  t->Branch("byLooseIsolationMVArun2v1DBnewDMwLT_2", &byLooseIsolationMVArun2v1DBnewDMwLT_2);
  t->Branch("byMediumIsolationMVArun2v1DBnewDMwLT_2", &byMediumIsolationMVArun2v1DBnewDMwLT_2);
  t->Branch("byTightIsolationMVArun2v1DBnewDMwLT_2", &byTightIsolationMVArun2v1DBnewDMwLT_2);
  t->Branch("byVTightIsolationMVArun2v1DBnewDMwLT_2", &byVTightIsolationMVArun2v1DBnewDMwLT_2);

  t->Branch("byRerunMVAIdVLoose_2", &NewMVAIDVLoose_2);
  t->Branch("byRerunMVAIdLoose_2", &NewMVAIDLoose_2);
  t->Branch("byRerunMVAIdMedium_2", &NewMVAIDMedium_2);
  t->Branch("byRerunMVAIdTight_2", &NewMVAIDTight_2);
  t->Branch("byRerunMVAIdVTight_2", &NewMVAIDVTight_2);
  t->Branch("byRerunMVAIdVVTight_2", &NewMVAIDVVTight_2);

  t->Branch("idMVANewDM_2", &idMVANewDM_2);
  t->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2);
  t->Branch("neutralIsoPtSum_2", &neutralIsoPtSum_2);
  t->Branch("puCorrPtSum_2", &puCorrPtSum_2);
  t->Branch("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2);
  t->Branch("decayMode_2", &decayMode_2);

  t->Branch("pzetavis", &pzetavis);
  t->Branch("pzetamiss", &pzetamiss);
  t->Branch("dzeta", &dzeta);
  
  t->Branch("pt_tt", &pt_tt);
  t->Branch("pt_vis", &pt_vis);
  t->Branch("dphi", &dphi);
  t->Branch("mt_tot", &mt_tot);
  t->Branch("pfpt_tt", &pfpt_tt);
  t->Branch("m_vis", &m_vis);
  t->Branch("m_coll", &m_coll);

  t->Branch("eleTauFakeRateWeight", &eleTauFakeRateWeight);
  t->Branch("muTauFakeRateWeight", &muTauFakeRateWeight);

  t->Branch("passesIsoCuts", &passesIsoCuts);
  t->Branch("passesLepIsoCuts", &passesLepIsoCuts);
  t->Branch("passesTauLepVetos", &passesTauLepVetos);
  t->Branch("passesThirdLepVeto", &passesThirdLepVeto);
  t->Branch("passesDiMuonVeto", &passesDiMuonVeto);
  t->Branch("passesDiElectronVeto", &passesDiElectronVeto);

  t->Branch("matchXTrig_obj", &matchXTrig_obj);
  t->Branch("dilepton_veto", &dilepton_veto);
  t->Branch("extraelec_veto", &extraelec_veto);
  t->Branch("extramuon_veto", &extramuon_veto);
  t->Branch("uncorrmet", &uncorrmet );
  t->Branch("met", &met);
  t->Branch("metphi", &metphi);
  t->Branch("met_ex", &met_ex);
  t->Branch("met_ey", &met_ey);
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

  t->Branch("m_sv", &m_sv);
  t->Branch("pt_sv", &pt_sv);

  t->Branch("mjj", &mjj);
  t->Branch("mjjUp", &mjjUp);
  t->Branch("mjjDown", &mjjDown);
  t->Branch("jdeta", &jdeta);
  t->Branch("jdetaUp", &jdetaUp);
  t->Branch("jdetaDown", &jdetaDown);
  t->Branch("njetingap", &njetingap);
  t->Branch("njetingap20", &njetingap20);
  t->Branch("dijetpt", &dijetpt);
  t->Branch("dijetphi", &dijetphi);
  t->Branch("jdphi", &jdphi);
  t->Branch("nbtag", &nbtag);
  t->Branch("njets", &njets);
  t->Branch("njetsUp", &njetsUp);
  t->Branch("njetsDown", &njetsDown);
  t->Branch("njetspt20", &njetspt20);
  t->Branch("jpt_1", &jpt_1);
  t->Branch("jptUp_1", &jptUp_1);
  t->Branch("jptDown_1", &jptDown_1);
  t->Branch("jeta_1", &jeta_1);
  t->Branch("jphi_1", &jphi_1);
  t->Branch("jm_1", &jm_1);
  t->Branch("jrawf_1", &jrawf_1);
  t->Branch("jmva_1", &jmva_1);
  t->Branch("jcsv_1", &jcsv_1);
  t->Branch("jpt_2", &jpt_2);
  t->Branch("jptUp_2", &jptUp_2);
  t->Branch("jptDown_2", &jptDown_2);
  t->Branch("jeta_2", &jeta_2);
  t->Branch("jphi_2", &jphi_2);
  t->Branch("jm_2", &jm_2);
  t->Branch("jrawf_2", &jrawf_2);
  t->Branch("jmva_2",&jmva_2);
  t->Branch("jcsv_2",&bcsv_2);
  t->Branch("bpt_1", &bpt_1);
  t->Branch("beta_1", &beta_1);
  t->Branch("bphi_1", &bphi_1);
  t->Branch("brawf_1",&brawf_1);
  t->Branch("bmva_1",&bmva_1);
  t->Branch("bcsv_1", &bcsv_1);
  t->Branch("bpt_2", &bpt_2);
  t->Branch("beta_2", &beta_2);
  t->Branch("bphi_2", &bphi_2);
  t->Branch("brawf_2",&brawf_2);
  t->Branch("bmva_2",&bmva_2);
  t->Branch("bcsv_2", &bcsv_2);

  if(!isSync){
    t->Branch("pfpt_sum", &pfpt_sum);
    t->Branch("pt_sum", &pt_sum);
    t->Branch("dr_leptau", &dr_leptau);

    t->Branch("jeta1eta2", &jeta1eta2);
    t->Branch("met_centrality", &met_centrality);
    t->Branch("lep_etacentrality", &lep_etacentrality);
    t->Branch("sphericity", &sphericity);

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
    t->Branch("addtau_byIsolationMVArun2v1DBnewDMwLTraw", &addtau_byIsolationMVArun2v1DBnewDMwLTraw);
    t->Branch("addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits", &addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits);
    t->Branch("addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
    t->Branch("addtau_byTightCombinedIsolationDeltaBetaCorr3Hits", &addtau_byTightCombinedIsolationDeltaBetaCorr3Hits);
    t->Branch("addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    t->Branch("addtau_byVLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byVLooseIsolationMVArun2v1DBoldDMwLT);
    t->Branch("addtau_byLooseIsolationMVArun2v1DBoldDMwLT", &addtau_byLooseIsolationMVArun2v1DBoldDMwLT);
    t->Branch("addtau_byMediumIsolationMVArun2v1DBoldDMwLT", &addtau_byMediumIsolationMVArun2v1DBoldDMwLT);
    t->Branch("addtau_byTightIsolationMVArun2v1DBoldDMwLT", &addtau_byTightIsolationMVArun2v1DBoldDMwLT);
    t->Branch("addtau_byVTightIsolationMVArun2v1DBoldDMwLT", &addtau_byVTightIsolationMVArun2v1DBoldDMwLT);
    t->Branch("addtau_byVLooseIsolationMVArun2v1DBnewDMwLT", &addtau_byVLooseIsolationMVArun2v1DBnewDMwLT);
    t->Branch("addtau_byLooseIsolationMVArun2v1DBnewDMwLT", &addtau_byLooseIsolationMVArun2v1DBnewDMwLT);
    t->Branch("addtau_byMediumIsolationMVArun2v1DBnewDMwLT", &addtau_byMediumIsolationMVArun2v1DBnewDMwLT);
    t->Branch("addtau_byTightIsolationMVArun2v1DBnewDMwLT", &addtau_byTightIsolationMVArun2v1DBnewDMwLT);
    t->Branch("addtau_byVTightIsolationMVArun2v1DBnewDMwLT", &addtau_byVTightIsolationMVArun2v1DBnewDMwLT);

    t->Branch("addtau_NewMVAIDVLoose", &addtau_NewMVAIDVLoose);
    t->Branch("addtau_NewMVAIDLoose", &addtau_NewMVAIDLoose);
    t->Branch("addtau_NewMVAIDMedium", &addtau_NewMVAIDMedium);
    t->Branch("addtau_NewMVAIDTight", &addtau_NewMVAIDTight);
    t->Branch("addtau_NewMVAIDVTight", &addtau_NewMVAIDVTight);
    t->Branch("addtau_NewMVAIDVVTight", &addtau_NewMVAIDVVTight);
  
    t->Branch("addtau_passesTauLepVetos", &addtau_passesTauLepVetos);
    t->Branch("addtau_decayMode", &addtau_decayMode);
    t->Branch("addtau_d0", &addtau_d0);
    t->Branch("addtau_dZ", &addtau_dZ);
    t->Branch("addtau_gen_match", &addtau_gen_match);
    t->Branch("addtau_mt", &addtau_mt);
    t->Branch("addtau_mvis", &addtau_mvis);
  }

}
