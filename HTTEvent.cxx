#ifdef PROJECT_NAME
#include "m2n/HTTDataFormats/interface/HTTEvent.h"
#include "m2n/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
// #include "EnergyScales.h"
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

  lastSystEffect = HTTAnalysis::NOMINAL;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTParticle::getP4(HTTAnalysis::sysEffects defaultType) const
{
    HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : corrType;

    if(type==HTTAnalysis::DUMMY_SYS)
    {
        lastSystEffect = type;
        return p4;

    }else
    {
        lastSystEffect = type;

        p4Cache = p4;
        int pdg = std::abs(getPDGid());
        int dm = getProperty(PropertyEnum::decayMode);
        int mc_match = getProperty(PropertyEnum::mc_match);


        float tauES = HTTAnalysis::getEnergyScale(pdg, mc_match, dm, type); // Defined in AnalysisEnums.h
        float tauES_mass = tauES;
        if (dm == 0) tauES_mass=0;

        p4Cache.SetPtEtaPhiM(p4.Pt() * (1.0+tauES),
                             p4.Eta(),
                             p4.Phi(),
                             p4.M() * (1.0+tauES_mass) );

        return p4Cache;
    }

    p4Cache = p4;

    return p4;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTParticle::getShiftedP4(float scale, bool preserveMass) const
{

    if(!preserveMass){
        p4Cache = scale*p4;
        return p4Cache;
    }

    double pt = p4.Perp();
    double energy =  p4.E();
    pt*=scale;
    double shiftedMomentum = pt/sin(p4.Theta());
    energy = sqrt(p4.M2() + pow(shiftedMomentum,2));
    p4Cache = p4;
    p4Cache.SetRho(shiftedMomentum);
    p4Cache.SetE(energy);
    return p4Cache;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HTTPair::clear()
{

    if(!p4Vector.size()) p4Vector.resize(HTTAnalysis::DUMMY_SYS);
    if(!leg1p4Vector.size()) leg1p4Vector.resize(HTTAnalysis::DUMMY_SYS);
    if(!leg2p4Vector.size()) leg2p4Vector.resize(HTTAnalysis::DUMMY_SYS);
    if(!svMetVector.size()) svMetVector.resize(HTTAnalysis::DUMMY_SYS);

    for(auto &it:p4Vector) it*=0;
    for(auto &it:leg1p4Vector) it*=0;
    for(auto &it:leg2p4Vector) it*=0;
    for(auto &it:svMetVector) it*=0;

    metMatrix.clear();

    leg1.clear();
    leg2.clear();

    lastSystEffect = HTTAnalysis::NOMINAL;
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
const TLorentzVector & HTTPair::getP4(HTTAnalysis::sysEffects defaultType) const
{
    HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;
    if(p4Vector.size()>(unsigned int)defaultType) return p4Vector[(unsigned int)defaultType];
    return p4Vector[(unsigned int)HTTAnalysis::NOMINAL];
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTPair::getLeg1P4(HTTAnalysis::sysEffects defaultType) const
{
    HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;
    if(leg1p4Vector.size()>(unsigned int)defaultType) return leg1p4Vector[(unsigned int)defaultType];
    return leg1p4Vector[(unsigned int)HTTAnalysis::NOMINAL];
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTPair::getLeg2P4(HTTAnalysis::sysEffects defaultType) const
{
    HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;
    if(leg2p4Vector.size()>(unsigned int)defaultType) return leg2p4Vector[(unsigned int)defaultType];
    return leg2p4Vector[(unsigned int)HTTAnalysis::NOMINAL];
}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TVector2 & HTTPair::getSystScaleMET(HTTAnalysis::sysEffects defaultType) const
{

    HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;

    if(type==HTTAnalysis::NOMINAL
      || (unsigned int)type>(unsigned int)HTTAnalysis::DUMMY_SYS)
    {
      
        lastSystEffect = type;
        if( (std::abs(leg1.getPDGid())==15) 
            || (std::abs(leg2.getPDGid())==15) )
        {

            double metX = met.X();
            metX+=leg1.getP4(HTTAnalysis::DUMMY_SYS).X(); //uncor
            metX+=leg2.getP4(HTTAnalysis::DUMMY_SYS).X(); //uncor
            metX-=leg1.getP4(HTTAnalysis::NOMINAL).X();
            metX-=leg2.getP4(HTTAnalysis::NOMINAL).X();

            double metY = met.Y();
            metY+=leg1.getP4(HTTAnalysis::DUMMY_SYS).Y();
            metY+=leg2.getP4(HTTAnalysis::DUMMY_SYS).Y();
            metY-=leg1.getP4(HTTAnalysis::NOMINAL).Y();
            metY-=leg2.getP4(HTTAnalysis::NOMINAL).Y();

            metCache.SetX(metX);
            metCache.SetY(metY);
            return metCache;
        }
        return met;
    } else if(lastSystEffect==type) return metCache;

    double metX = met.X();
    metX+=leg1.getP4(HTTAnalysis::DUMMY_SYS).X();
    metX+=leg2.getP4(HTTAnalysis::DUMMY_SYS).X();
    metX-=leg1.getP4(type).X();
    metX-=leg2.getP4(type).X();

    double metY = met.Y();
    metY+=leg1.getP4(HTTAnalysis::DUMMY_SYS).Y();
    metY+=leg2.getP4(HTTAnalysis::DUMMY_SYS).Y();
    metY-=leg1.getP4(type).Y();
    metY-=leg2.getP4(type).Y();

    metCache.SetX(metX);
    metCache.SetY(metY);

    return metCache;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
float HTTPair::getMTTOT(HTTAnalysis::sysEffects defaultType) const
{
  HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;

  float mt1 = 2. * leg1.getP4(type).Pt() * getSystScaleMET(type).Mod() * (1. - TMath::Cos(leg1.getP4(type).Phi()-getSystScaleMET(type).Phi()));
  float mt2 = 2. * leg2.getP4(type).Pt() * getSystScaleMET(type).Mod() * (1. - TMath::Cos(leg2.getP4(type).Phi()-getSystScaleMET(type).Phi()));
  float mt3 = 2. * leg1.getP4(type).Pt() * leg2.getP4(type).Pt() * (1. - TMath::Cos(leg1.getP4(type).Phi()-leg2.getP4(type).Phi()));
  return TMath::Sqrt( mt1 + mt2 + mt3 );
}

float HTTPair::getPT_TT(HTTAnalysis::sysEffects defaultType) const
{
  HTTAnalysis::sysEffects type =  defaultType != HTTAnalysis::NOMINAL ? defaultType : HTTParticle::corrType;

  TLorentzVector vmet; 
  vmet.SetPtEtaPhiM(getSystScaleMET(type).Mod(),0,getSystScaleMET(type).Phi(),0);

  return (vmet + leg1.getP4(type) + leg2.getP4(type) ).Pt();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HTTJet::clear()
{

    p4.SetPtEtaPhiM(-10.,-10.,-10.,-10.);
    currentP4.SetPtEtaPhiM(-10.,-10.,-10.,-10.);

    jecUncertSourceValues.clear();
    jecUncertSourceValues.reserve( (unsigned int)JecUncertEnum::NONE );

    // Fallback when jec uncerts get asymmetric
    // jecUncertSourceValuesUp.clear();
    // jecUncertSourceValuesUp.reserve( (unsigned int)JecUncertEnum::NONE );

    // jecUncertSourceValuesDown.clear();
    // jecUncertSourceValuesDown.reserve( (unsigned int)JecUncertEnum::NONE );

}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTJet::getP4(JecUncertEnum uncert, bool up)
{
    currentP4 = p4;
    if(uncert == JecUncertEnum::NONE) return currentP4;

    double shift = up ? 1. : -1.;   
    double scale = jecUncertSourceValues[(unsigned int)uncert]; 

    // Fallback when jec uncerts get asymmetric
    // if(currentShift)
    // {
    //     scale = jecUncertSourceValuesUp[(unsigned int)currentUncert];
    //     shift = 1.;

    // }else
    // {
    //     scale = jecUncertSourceValuesDown[(unsigned int)currentUncert];
    //     shift = -1.;
    // }

    currentP4.SetPtEtaPhiM( p4.Pt()*(1 + scale*shift),
                            p4.Eta(),
                            p4.Phi(),
                            p4.M()*(1 + scale*shift)
                          );
    return currentP4;

}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJet::setJecUncertSourceValue(unsigned int uncert, double value, bool up)
{
    jecUncertSourceValues[uncert] = value;

    // Fallback when jec uncerts get asymmetric
    // if(up) jecUncertSourceValuesUp[uncert] = value;
    // else   jecUncertSourceValuesDown[uncert] = value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HTTJetCollection::clear()
{
    usePromoteDemote = false;

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

  eff_file = new TFile("utils/BTagCalibration/data/tagging_efficiencies_Moriond2017.root");
  hb_eff = dynamic_cast<TH2F*>(eff_file->Get("btag_eff_b") );
  hc_eff = dynamic_cast<TH2F*>(eff_file->Get("btag_eff_c") );
  hoth_eff = dynamic_cast<TH2F*>(eff_file->Get("btag_eff_oth") );

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
void HTTJetCollection::fillCurrentCollections(JecUncertEnum uncert, bool up)
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

