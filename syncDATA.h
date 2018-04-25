#include "HTTEvent.h"

#ifndef __syncDATA__
#define __syncDATA__

class syncDATA
{
 public:
  //ClassDef(syncDATA,0);

  Float_t lumiWeight;
  Int_t run_syncro;
  Float_t lumi_syncro;
  ULong64_t evt_syncro;
  int entry;
  Int_t fileEntry;
  bool matchXTrig_obj;

  int npv;
  int npvGood;
  float npu;
  float rho;
  int gen_match_1;
  int gen_match_2;
  float genPt_1;
  float genPt_2;
  int gen_match_jetId_1;
  int gen_match_jetId_2;
  int NUP;
  //float evtWeight;
  float weight;
  float puWeight;
  float genWeight;
  float trigweight_1;
  float anti_trigweight_1;
  float trigweight_2;
  float idisoweight_1;
  float anti_idisoweight_1;
  float idisoweight_2;
  float trk_sf;
  float effweight;
  float stitchedWeight;
  float topWeight;
  float topWeight_run1;
  float ZWeight;
  float zpt_weight_nom;
  float zpt_weight_esup;
  float zpt_weight_esdown;
  float zpt_weight_ttup;
  float zpt_weight_ttdown;
  float zpt_weight_statpt0up;
  float zpt_weight_statpt0down;
  float zpt_weight_statpt40up;
  float zpt_weight_statpt40down;
  float zpt_weight_statpt80up;
  float zpt_weight_statpt80down;
  float gen_Mll;
  float gen_ll_px;
  float gen_ll_py;
  float gen_ll_pz;
  float gen_vis_Mll;
  float gen_ll_vis_px;
  float gen_ll_vis_py;
  float gen_ll_vis_pz;
  float gen_top_pt_1;
  float gen_top_pt_2;
  int genJets;
  int genJet_match_1;
  int genJet_match_2;
  float matchedJetPt03_1;
  float matchedJetPt05_1;
  float matchedJetPt03_2;
  float matchedJetPt05_2;
  //////////////////////////////////////////////////////////////////  
  int trg_singlemuon; //fires OR of HLT_IsoMu22, HLT_IsoTkMu22, HLT_IsoMu22eta2p1, HLT_IsoTkMu22_eta2p1
  int trg_mutaucross;
  int trg_singleelectron; //fires HLT_Ele25_eta2p1_WPTight_Gsf
  int trg_singletau; //fires HLT_VLooseIsoPFTau120_Trk50_eta2p1
  int trg_doubletau; //fires HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg or HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg
  int trg_muonelectron; //fires HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
  int passBadMuonFilter;
  int passBadChargedHadronFilter;
  int Flag_HBHENoiseFilter;
  int Flag_HBHENoiseIsoFilter;
  int Flag_EcalDeadCellTriggerPrimitiveFilter;
  int Flag_goodVertices;
  int Flag_eeBadScFilter;
  int Flag_globalTightHalo2016Filter;
  int failBadGlobalMuonTagger;
  int failCloneGlobalMuonTagger;
  int Flag_duplicateMuons;
  int Flag_badMuons;
  int Flag_noBadMuons;
  int passesMetMuonFilter;
  //////////////////////////////////////////////////////////////////  
  float pt_1;
  float phi_1;
  float eta_1;
  float eta_SC_1;
  float m_1;
  int q_1;
  float d0_1;
  float dZ_1;
  float mt_1;
  float iso_1;
  int againstElectronLooseMVA6_1;
  int againstElectronMediumMVA6_1;
  int againstElectronTightMVA6_1;
  int againstElectronVLooseMVA6_1;
  int againstElectronVTightMVA6_1;
  int againstMuonLoose3_1;
  int againstMuonTight3_1;
  float byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
  int byLooseCombinedIsolationDeltaBetaCorr3Hits_1;
  int byMediumCombinedIsolationDeltaBetaCorr3Hits_1;
  int byTightCombinedIsolationDeltaBetaCorr3Hits_1;
  int byIsolationMVA3newDMwoLTraw_1;
  int byIsolationMVA3oldDMwoLTraw_1;
  float byIsolationMVA3newDMwLTraw_1;
  float byIsolationMVA3oldDMwLTraw_1;
  int byVLooseIsolationMVArun2v1DBoldDMwLT_1;
  int byLooseIsolationMVArun2v1DBoldDMwLT_1;
  int byMediumIsolationMVArun2v1DBoldDMwLT_1;
  int byTightIsolationMVArun2v1DBoldDMwLT_1;
  int byVTightIsolationMVArun2v1DBoldDMwLT_1;
  int byVLooseIsolationMVArun2v1DBnewDMwLT_1;
  int byLooseIsolationMVArun2v1DBnewDMwLT_1;
  int byMediumIsolationMVArun2v1DBnewDMwLT_1;
  int byTightIsolationMVArun2v1DBnewDMwLT_1;
  int byVTightIsolationMVArun2v1DBnewDMwLT_1;
  int NewMVAIDVLoose_1;
  int NewMVAIDLoose_1;
  int NewMVAIDMedium_1;
  int NewMVAIDTight_1;
  int NewMVAIDVTight_1;
  int NewMVAIDVVTight_1;
  float idMVANewDM_1;
  float chargedIsoPtSum_1;
  float neutralIsoPtSum_1;
  float puCorrPtSum_1;
  int decayModeFindingOldDMs_1;
  int decayMode_1;
  float id_e_mva_nt_loose_1;
  float id_m_loose_1;
  float id_m_medium_1;
  float id_m_tight_1;
  float id_m_tightnovtx_1;
  float id_m_highpt_1;
  float id_e_cut_veto_1;
  float id_e_cut_loose_1;
  float id_e_cut_medium_1;
  float id_e_cut_tight_1;
  
  //////////////////////////////////////////////////////////////////
  float pt_2;
  float phi_2;
  float eta_2;
  float m_2;
  int q_2;
  float d0_2;
  float dZ_2;
  float mt_2;
  float iso_2;
  int againstElectronLooseMVA6_2;
  int againstElectronMediumMVA6_2;
  int againstElectronTightMVA6_2;
  int againstElectronVLooseMVA6_2;
  int againstElectronVTightMVA6_2;
  int againstMuonLoose3_2;
  int againstMuonTight3_2;
  float byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
  int byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
  int byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
  int byTightCombinedIsolationDeltaBetaCorr3Hits_2;
  int byIsolationMVA3newDMwoLTraw_2;
  int byIsolationMVA3oldDMwoLTraw_2;
  float byIsolationMVA3newDMwLTraw_2;
  float byIsolationMVA3oldDMwLTraw_2;
  int byVLooseIsolationMVArun2v1DBoldDMwLT_2;
  int byLooseIsolationMVArun2v1DBoldDMwLT_2;
  int byMediumIsolationMVArun2v1DBoldDMwLT_2;
  int byTightIsolationMVArun2v1DBoldDMwLT_2;
  int byVTightIsolationMVArun2v1DBoldDMwLT_2;
  int byVLooseIsolationMVArun2v1DBnewDMwLT_2;
  int byLooseIsolationMVArun2v1DBnewDMwLT_2;
  int byMediumIsolationMVArun2v1DBnewDMwLT_2;
  int byTightIsolationMVArun2v1DBnewDMwLT_2;
  int byVTightIsolationMVArun2v1DBnewDMwLT_2;
  int NewMVAIDVLoose_2;
  int NewMVAIDLoose_2;
  int NewMVAIDMedium_2;
  int NewMVAIDTight_2;
  int NewMVAIDVTight_2;
  int NewMVAIDVVTight_2;
  float idMVANewDM_2;
  float chargedIsoPtSum_2;
  float neutralIsoPtSum_2;
  float puCorrPtSum_2;
  int decayModeFindingOldDMs_2;
  int decayMode_2;
  //////////////////////////////////////////////////////////////////
  int nbtag;
  int njets;
  int njetsUp;
  int njetsDown;
  int njetspt20;
  float mjj;
  float mjjUp;
  float mjjDown;
  float jdeta;
  float jdetaUp;
  float jdetaDown;
  int njetingap;
  int njetingap20;
  float dijetpt;
  float dijetphi;
  float jdphi;
  float jpt_1;
  float jptUp_1;
  float jptDown_1;
  float jeta_1;
  float jphi_1;
  float jm_1;
  float jrawf_1;
  float jmva_1;
  float jcsv_1;
  float jpt_2;
  float jptUp_2;
  float jptDown_2;
  float jeta_2;
  float jphi_2;
  float jm_2;
  float jrawf_2;
  float jmva_2;
  float jcsv_2;
  float bpt_1;
  float beta_1;
  float bphi_1;
  float brawf_1;
  float bmva_1;
  float bcsv_1;
  float bpt_2;
  float beta_2;
  float bphi_2;
  float brawf_2;
  float bmva_2;
  float bcsv_2;
  //////////////////////////////////////////////////////////////////
  float met;
  float uncorrmet;
  float met_ex;
  float met_ey;
  float metphi;
  float corrmet;
  float corrmet_ex;
  float corrmet_ey;
  float corrmetphi;
  float mvamet;
  float mvamet_ex;
  float mvamet_ey;
  float mvametphi;
  float corrmvamet_ex;
  float corrmvamet_ey;
  float corrmvamet;
  float corrmvametphi;
  float mvacov00;
  float mvacov01;
  float mvacov10;
  float mvacov11;
  float metcov00;
  float metcov01;
  float metcov10;
  float metcov11;
  float m_sv;
  float pt_sv;
  //////////////////////////////////////////////////////////////////
  float eleTauFakeRateWeight;
  float muTauFakeRateWeight;
  float antilep_tauscaling;
  //////////////////////////////////////////////////////////////////
  bool passesIsoCuts;
  bool passesLepIsoCuts;
  bool passesTauLepVetos;
  bool passesThirdLepVeto;
  bool passesDiMuonVeto;
  bool passesDiElectronVeto;
  //////////////////////////////////////////////////////////////////
  bool dilepton_veto;
  bool extramuon_veto;
  bool extraelec_veto;
  //////////////////////////////////////////////////////////////////
  float pzetavis;
  float pzetamiss;
  float dzeta;
  //////////////////////////////////////////////////////////////////
  float pt_tt;
  float pt_vis;
  float mt_3;
  float mt_tot;
  float pfpt_tt;
  float m_vis;
  float m_coll;
  float dphi;
  //////////////////////////////////////////////////////////////////
  float pfmt_1;
  float pfmt_2;
  float pfpt_sum;
  float pt_sum;
  float dr_leptau;
  float jeta1eta2;
  float met_centrality;
  float mvamet_centrality;
  float lep_etacentrality;
  float sphericity;

  syncDATA(){}  
  ~syncDATA(){}  

  void setDefault();
  void fill(HTTEvent *ev);

  /*
  
  //////////////////////////////////////////////////////////////////
  int nadditionalMu;
  vector<double> addmuon_pt;
  vector<double> addmuon_eta;
  vector<double> addmuon_phi;
  vector<double> addmuon_m;
  vector<int> addmuon_q;
  vector<double> addmuon_iso;
  vector<int> addmuon_gen_match;
  //////////////////////////////////////////////////////////////////
  int nadditionalEle;
  vector<double> addele_pt;
  vector<double> addele_eta;
  vector<double> addele_phi;
  vector<double> addele_m;
  vector<int> addele_q;
  vector<double> addele_iso;
  vector<int> addele_gen_match;
  //////////////////////////////////////////////////////////////////
  int nadditionalTau;
  vector<double> addtau_pt;
  vector<double> addtau_eta;
  vector<double> addtau_phi;
  vector<double> addtau_m;
  vector<double> addtau_q;
  vector<double> addtau_byIsolationMVArun2v1DBnewDMwLTraw;
  vector<double> addtau_byCombinedIsolationDeltaBetaCorrRaw3Hits;
  vector<int> addtau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
  vector<int> addtau_byTightCombinedIsolationDeltaBetaCorr3Hits;
  vector<int> addtau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  vector<int> addtau_byVLooseIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byLooseIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byMediumIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byTightIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byVTightIsolationMVArun2v1DBoldDMwLT;
  vector<int> addtau_byVLooseIsolationMVArun2v1DBnewDMwLT;
  vector<int> addtau_byLooseIsolationMVArun2v1DBnewDMwLT;
  vector<int> addtau_byMediumIsolationMVArun2v1DBnewDMwLT;
  vector<int> addtau_byTightIsolationMVArun2v1DBnewDMwLT;
  vector<int> addtau_byVTightIsolationMVArun2v1DBnewDMwLT;
  vector<int> addtau_NewMVAIDVLoose;
  vector<int> addtau_NewMVAIDLoose;
  vector<int> addtau_NewMVAIDMedium;
  vector<int> addtau_NewMVAIDTight;
  vector<int> addtau_NewMVAIDVTight;
  vector<int> addtau_NewMVAIDVVTight;
  vector<int> addtau_passesTauLepVetos;
  vector<int> addtau_decayMode;
  vector<double> addtau_d0;
  vector<double> addtau_dZ;
  vector<int> addtau_gen_match;
  vector<double> addtau_mt;
  vector<double> addtau_mvis;
  //////////////////////////////////////////////////////////////////
  BTagCalibration calib;
  BTagCalibrationReader reader;
  int whichDilepton;
  int whichDilepton_addTaus;
  vector<float> top_pt;
  vector<RecObj> v_mu;
  vector<RecObj> v_trailMu;
  vector<RecObj> v_mu_add;
  vector<RecObj> v_el_add;
  vector<RecObj> v_mu_veto;
  vector<RecObj> v_tau;
  vector<RecObj> v_el;
  vector<RecObj> v_el_veto;
  vector<RecObj> v_jet;
  vector<RecObj> v_jetUp;
  vector<RecObj> v_jetDown;
  vector<RecObj> v_bjet;
  vector<RecObj> v_antibjet;
  // vector<RecObj> v_ele_pos;
  // vector<RecObj> v_ele_neg;
  // vector<RecObj> v_mu_pos;
  // vector<RecObj> v_mu_neg;
  RecObj s_leg1;
  RecObj s_leg2;
  RecObj s_jet;
  RecObj s_bjet;
  RecObj s_MVAmet;
  RecObj s_PFmet_NoRecoilCorr;
  RecObj s_PFmet;

  */

};
#endif
