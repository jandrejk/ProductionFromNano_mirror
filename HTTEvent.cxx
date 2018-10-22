#ifdef PROJECT_NAME
#include "m2n/HTTDataFormats/interface/HTTEvent.h"
#include "m2n/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#endif

std::map<string,PropertyEnum> HTTEvent::usePropertyFor = {};
HTTAnalysis::sysEffects HTTParticle::corrType = HTTAnalysis::NOMINAL;

////////////////////////////////////////////////
////////////////////////////////////////////////
const TString HTTEvent::tauIDStrings[ntauIds] = {
  "againstMuonLoose3",
  "againstMuonTight3",
  "againstElectronVLooseMVA6",
  "againstElectronLooseMVA6",
  "againstElectronMediumMVA6",
  "againstElectronTightMVA6",
  "againstElectronVTightMVA6",
  "byVLooseIsolationMVArun2v1DBoldDMwLT",
  "byLooseIsolationMVArun2v1DBoldDMwLT",
  "byMediumIsolationMVArun2v1DBoldDMwLT",
  "byTightIsolationMVArun2v1DBoldDMwLT",
  "byVTightIsolationMVArun2v1DBoldDMwLT",
  "byVVTightIsolationMVArun2v1DBoldDMwLT"
};
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTEvent::clear(){

  runId = 0;
  eventId = 0;
  lsId = 0;

  nPU = 0;
  nPV = 0;
  rho = 0;
  puWeight = 1.;

  mcWeight = 1.0;

  isSetTopP4 = false;
  topPtReWeight =1.0;
  topPtReWeightR1 =1.0;

  isSetGenBosonP4 = false;
  zPtReWeight =1.0;
  zPtReWeightSUSY =1.0;

  lheHt = 1.0;
  lheNOutPartons = 0;
  aMCatNLOweight = 1.0;

  genPV*=0;
  AODPV*=0;
  refittedPV*=0;

  met*=0;
  bosP4*=0;
  bosVisP4*=0;

  isRefit = false;

  nTracksInRefit = 0;

  metFilterDecision = 0;
  selectionWord.ResetAllBits();



}
void HTTEvent::setSampleType(string sampletype)
{

  if(sampletype == "SingleMuon")     sampleType = MuonData;
  if(sampletype == "SingleElectron") sampleType = EleData;
  if(sampletype == "Tau")            sampleType = TauData;
  if(sampletype == "dy")             sampleType = DY;
  if(sampletype == "dy_lowmass")     sampleType = DYLowM;
  if(sampletype == "w")              sampleType = WJets;
  if(sampletype == "tt")             sampleType = TTbar;
  if(sampletype == "st")             sampleType = ST;
  if(sampletype == "diboson")        sampleType = Diboson;
  if(sampletype == "ewk")            sampleType = EWK;   
  if(sampletype == "signal")         sampleType = h;
  sampletype = DUMMY;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HTTParticle::clear(){

  p4*=0;
  chargedP4*=0;
  neutralP4*=0;

  pca*=0;
  pcaRefitPV*=0;
  pcaGenPV*=0;

  properties.clear();
  deltaVector.SetMagPhi(0.,0.);

  lastSystEffect = HTTAnalysis::NOMINAL;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTParticle::setP4(const TLorentzVector &aP4 )
{
    p4 = aP4;
    currentP4 = getShiftedP4(corrType);
    lastSystEffect = corrType;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTParticle::getP4(HTTAnalysis::sysEffects defaultType) const
{
    HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : corrType;

    if(type == lastSystEffect) return currentP4;
    lastSystEffect = type;
    currentP4 = getShiftedP4(type); // Set currentP4 before returning
    return currentP4;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector HTTParticle::getShiftedP4(HTTAnalysis::sysEffects shift) const
{

  TLorentzVector shiftedP4 = p4;
  int pdg = std::abs(getPDGid());
  int dm = getProperty(PropertyEnum::decayMode);
  int mc_match = getProperty(PropertyEnum::mc_match);

  float tauES = HTTAnalysis::getEnergyScale(pdg, mc_match, dm, shift); // Defined in AnalysisEnums.h
  float tauES_mass = tauES;
  if (dm == 0) tauES_mass=0;

  shiftedP4.SetPtEtaPhiM(p4.Pt() * (1.0+tauES),
                       p4.Eta(),
                       p4.Phi(),
                       p4.M() * (1.0+tauES_mass) );

  deltaVector.SetMagPhi( (shiftedP4 - p4).Pt(), p4.Phi() );

  return shiftedP4;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HTTPair::clear()
{
    met.clear();
    met[""] = TVector2(0.,0.);
    metCache = TVector2(0.,0.);

    p4Vector.clear();
    p4Vector[""].SetPtEtaPhiM(-10,-10,-10,-10);
    metMatrix.clear();

    leg1.clear();
    leg2.clear();

    lastMETShift = "";
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector HTTPair::getP4(P4Type type)
{
    if(type == P4Type::SVFit)
    {
      if( p4Vector.find(lastMETShift) == p4Vector.end() ) return p4Vector.at("");
      return p4Vector.at(lastMETShift);

    } else if(type == P4Type::Simple)
    {
      TLorentzVector vmet; 
      vmet.SetPtEtaPhiM(getMET().Mod(),0,getMET().Phi(),0);

      return (vmet + leg1.getP4() + leg2.getP4() );
    } else
    {
      return ( leg1.getP4() + leg2.getP4() );
    }
}
////////////////////////////////////////////////
////////////////////////////////////////////////
HTTAnalysis::finalState HTTPair::getFinalState()
{
    int pdg1=std::abs(getLeg1().getProperty(PropertyEnum::pdgId));
    int pdg2=std::abs(getLeg2().getProperty(PropertyEnum::pdgId));

    if(pdg1 == 11 && pdg2 == 15) return HTTAnalysis::EleTau;
    if(pdg1 == 13 && pdg2 == 15) return HTTAnalysis::MuTau;
    if(pdg1 == 15 && pdg2 == 15) return HTTAnalysis::TauTau;
    return HTTAnalysis::NONE;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
bool HTTPair::isInLooseSR()
{
    bool mt = getMTLeg1() < 50;
    bool os = leg1.getCharge()*leg2.getCharge() < 0;
    bool tauiso = true;

    HTTAnalysis::finalState fs = getFinalState();
    if( fs == HTTAnalysis::EleTau || fs == HTTAnalysis::MuTau )
    {
      tauiso = ((unsigned int)leg2.getProperty( HTTEvent::usePropertyFor.at("tauID") ) & 0x10) == 0x10;
    } else if( fs == HTTAnalysis::TauTau )
    {
      tauiso = ((unsigned int)leg1.getProperty( HTTEvent::usePropertyFor.at("tauID") ) & 0x10) == 0x10;
      tauiso &= ((unsigned int)leg2.getProperty( HTTEvent::usePropertyFor.at("tauID") ) & 0x10) == 0x10;
    }

    return (mt && os && tauiso);
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTPair::setCurrentMETShift(string uncert)
{
  metCache = met.at(uncert);
  lastMETShift = uncert;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
float HTTPair::getMTTOT(HTTAnalysis::sysEffects defaultType) const
{
  HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;

  float mt1 = 2. * leg1.getP4(type).Pt() * getMET().Mod() * (1. - TMath::Cos(leg1.getP4(type).Phi()-getMET().Phi()));
  float mt2 = 2. * leg2.getP4(type).Pt() * getMET().Mod() * (1. - TMath::Cos(leg2.getP4(type).Phi()-getMET().Phi()));
  float mt3 = 2. * leg1.getP4(type).Pt() * leg2.getP4(type).Pt() * (1. - TMath::Cos(leg1.getP4(type).Phi()-leg2.getP4(type).Phi()));
  return TMath::Sqrt( mt1 + mt2 + mt3 );
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HTTJet::clear()
{
    p4.SetPtEtaPhiM(-10.,-10.,-10.,-10.);
    currentP4.SetPtEtaPhiM(-10.,-10.,-10.,-10.);

    jecUncertSourceValues.clear();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTJet::getP4(string uncert, bool up)
{
    currentP4 = p4;
    if(uncert == "") return currentP4;

    double shift = up ? 1. : -1.;   
    double scale = jecUncertSourceValues[uncert]; 

    currentP4.SetPtEtaPhiM( p4.Pt()*(1 + scale*shift),
                            p4.Eta(),
                            p4.Phi(),
                            p4.M()*(1 + scale*shift)
                          );
    return currentP4;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HTTJetCollection::initCollection(bool isMC, bool applyRecoil, bool isSync)
{
  jecShifts = { {"",{"",true}} };

  if( HTTParticle::corrType == HTTAnalysis::NOMINAL && isMC && !applyRecoil && !isSync)
  {
      for(auto uncert : JecAfterSplitting)
      {
          for(auto shift : { make_pair("Up",true), make_pair("Down",false) } )
          {
            jecShifts.push_back( make_pair(uncert.first + shift.first, make_pair(uncert.first, shift.second) ) );
          }
      }
  }
  clear();
}
void HTTJetCollection::clear()
{
    dijet.SetPtEtaPhiM(-10.,-10.,-10.,-10.);

    jetCollection.clear();
    jetCurrentCollection.clear();
    btagCurrentCollection.clear();
    antibtagCurrentCollection.clear();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::initForPromoteDemote()
{
  usePromoteDemote = true;
  calib = BTagCalibration("DeepCSV", "utils/BTagCalibration/data/DeepCSV_94XSF_V3_B_F.csv");
  reader = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central",{"up","down"});
  reader.load(calib,  BTagEntry::FLAV_B, "comb");
  reader.load(calib,  BTagEntry::FLAV_C, "comb");
  reader.load(calib,  BTagEntry::FLAV_UDSG, "incl");

  eff_file = new TFile("utils/BTagCalibration/data/tagging_efficiencies_march2018_btageff-all_samp-inc-DeepCSV_medium.root");
  hb_eff = dynamic_cast<TH2D*>(eff_file->Get("btag_eff_b") );
  hc_eff = dynamic_cast<TH2D*>(eff_file->Get("btag_eff_c") );
  hoth_eff = dynamic_cast<TH2D*>(eff_file->Get("btag_eff_oth") );

}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::btagPromoteDemote(){

  TRandom3 rand;

  double scaleFactor;
  double scaleRatio;
  double jetPt;
  double jetEta;
  double tagging_efficiency;
  double rn;
  
  int index = 0;
  for(auto & jet : btagCurrentCollection){
    jetPt=  jet.Pt();    
    jetEta= jet.Eta();
    rand.SetSeed((int)((jetEta+5)*100000));
    rn = rand.Rndm();


    if( jet.getProperty(PropertyEnum::hadronFlavour) ==5 )       scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, jetEta, jetPt);
    else if( jet.getProperty(PropertyEnum::hadronFlavour) ==4 )  scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, jetEta, jetPt);
    else                                                         scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetEta, jetPt);
    
    if(scaleFactor < 1){
      scaleRatio = fabs( 1-scaleFactor );
      if( rn < scaleRatio ) btagCurrentCollection.erase(btagCurrentCollection.begin() + index);
    }
    ++index;
  }
    
  for(auto & jet : antibtagCurrentCollection){
    jetPt=  jet.Pt();
    jetEta= jet.Eta();

    rand.SetSeed((int)((jetEta+5)*100000));
    rn = rand.Rndm();

    if( jet.getProperty(PropertyEnum::hadronFlavour)==5 ){
      
      scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, jetEta, jetPt);
      if(jetPt>hb_eff->GetXaxis()->GetBinLowEdge(hb_eff->GetNbinsX()+1))     tagging_efficiency = hb_eff->GetBinContent(hb_eff->GetNbinsX(),hb_eff->GetYaxis()->FindBin(fabs( jetEta ) )); 
      else                                                                   tagging_efficiency = hb_eff->GetBinContent(hb_eff->GetXaxis()->FindBin(jetPt),hb_eff->GetYaxis()->FindBin(fabs( jetEta )));
    }
    else if( jet.getProperty(PropertyEnum::hadronFlavour)==4 ){

      scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, jetEta, jetPt);

      if(jetPt>hc_eff->GetXaxis()->GetBinLowEdge(hc_eff->GetNbinsX()+1))     tagging_efficiency = hc_eff->GetBinContent(hc_eff->GetNbinsX(),hc_eff->GetYaxis()->FindBin(fabs( jetEta )));
      else                                                                   tagging_efficiency = hc_eff->GetBinContent(hc_eff->GetXaxis()->FindBin(jetPt),hc_eff->GetYaxis()->FindBin(fabs( jetEta )));
    }
    else{
 
      scaleFactor = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, jetEta, jetPt);

      if(jetPt>hoth_eff->GetXaxis()->GetBinLowEdge(hoth_eff->GetNbinsX()+1)) tagging_efficiency = hoth_eff->GetBinContent(hoth_eff->GetNbinsX(),hoth_eff->GetYaxis()->FindBin(fabs( jetEta )));
      else                                                                   tagging_efficiency = hoth_eff->GetBinContent(hoth_eff->GetXaxis()->FindBin(jetPt),hoth_eff->GetYaxis()->FindBin(fabs( jetEta )));  
    }
    if(scaleFactor > 1){
      scaleRatio = fabs( ( 1-scaleFactor )/( 1-(1/ tagging_efficiency ) ) );
      if( rn < scaleRatio ) btagCurrentCollection.push_back( jet );
    }  
  }
  std::sort(btagCurrentCollection.begin(), btagCurrentCollection.end(), sortJets);
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::fillCurrentCollections(string uncert, bool up)
{
    jetCurrentCollection.clear();
    btagCurrentCollection.clear();
    antibtagCurrentCollection.clear();

    for(auto jet : jetCollection)
    {   
        jet.setUncertShift(uncert, up);
        if(jet.Pt() > 20)
          jetCurrentCollection.push_back(jet);

        if(jet.Pt() > 20 && abs(jet.Eta()) < 2.4 )
        {
          if( jet.getProperty(PropertyEnum::btagDeepB)>0.4941 )
          {
            btagCurrentCollection.push_back( jet );
          }else
          {
            antibtagCurrentCollection.push_back( jet );
          }
        }
    }
    if(usePromoteDemote) btagPromoteDemote();
    
    if(jetCurrentCollection.size() > 1)
        setDijetP4();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TVector2  HTTJetCollection::getTotalJetShift(string uncert, bool up)
{
  TLorentzVector before; before.SetPtEtaPhiE(0,0,0,0);
  TLorentzVector after; after.SetPtEtaPhiE(0,0,0,0);

  for(auto jet : jetCollection){
        before += jet.P4();
        jet.setUncertShift(uncert, up);
        after += jet.P4();
  }

  return TVector2( (after - before).Px(), (after - before).Py() )  ;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::setDijetP4()
{
    dijet = jetCurrentCollection.at(0).P4()+jetCurrentCollection.at(1).P4();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
int HTTJetCollection::getNJets(double pt)
{
    int njets = 0;
    for(auto &jet : jetCurrentCollection)
        if(jet.Pt() > pt )njets++;
    
    return njets;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
int HTTJetCollection::getNJetInGap(double pt)
{
    int njetingap = 0;
    if(jetCurrentCollection.size() < 2) return 0;

    double j1_eta =  jetCurrentCollection.at(0).Eta();
    double j2_eta =  jetCurrentCollection.at(1).Eta();

    for (unsigned ij=2; ij<jetCurrentCollection.size(); ij++)
    {
        TLorentzVector addJet = jetCurrentCollection.at(ij).P4();

        if ( addJet.Pt() < pt ) continue;

        float aj_eta=addJet.Eta();

        if( (aj_eta<j1_eta && aj_eta>j2_eta)
             || (aj_eta>j1_eta && aj_eta<j2_eta)
        )  njetingap++;

    }
    return njetingap;

}
////////////////////////////////////////////////
////////////////////////////////////////////////
bool HTTJetCollection::sortJets(HTTJet j1, HTTJet j2)
{
  if( j1.Pt() == j2.Pt() ) return true;
  return j1.Pt() > j2.Pt();
}

