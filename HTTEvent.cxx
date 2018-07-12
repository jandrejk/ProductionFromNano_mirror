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

  topPtReWeight =1.0;
  topPtReWeightR1 =1.0;
  zPtReWeight =1.0;
  zPtReWeightSUSY =1.0;

  lheHt = 1.0;
  lheNOutPartons = 0;
  aMCatNLOweight = 1.0;

  sampleType = DUMMY;

#ifdef PROJECT_NAME
  decayModeMinus = WawGenInfoHelper::tauDecayModes::tauDecayOther;
  decayModePlus = WawGenInfoHelper::tauDecayModes::tauDecayOther;
#else
  decayModeMinus  = 99;
  decayModePlus  = 99;
#endif


#ifdef PROJECT_NAME
  decayModeBoson = WawGenInfoHelper::bosonDecayModes::kUndefined;
#else
  decayModeBoson = -1;
#endif

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


        float tauES = HTTAnalysis::getEnergyScale(pdg, mc_match, dm, type);
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
    }
    else if(lastSystEffect==type) return metCache;

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

float HTTPair::getPTTOT(HTTAnalysis::sysEffects defaultType) const
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
    currentUncert = JecUncertEnum::NONE;
    currentShift = true;

    p4.SetPtEtaPhiM(DEF,DEF,DEF,DEF);
    jecUncertSourceValuesUp.clear();
    jecUncertSourceValuesUp.reserve( (unsigned int)JecUncertEnum::NONE );

    jecUncertSourceValuesDown.clear();
    jecUncertSourceValuesDown.reserve( (unsigned int)JecUncertEnum::NONE );

}
////////////////////////////////////////////////
////////////////////////////////////////////////
const TLorentzVector & HTTJet::getP4()
{
    double scale = 0.;
    double shift = 1.;

    currentP4 = p4;
    if(currentUncert == JecUncertEnum::NONE) return currentP4;

    if(currentShift)
    {
        scale = jecUncertSourceValuesUp[(unsigned int)currentUncert];
        shift = 1.;

    }else
    {
        scale = jecUncertSourceValuesDown[(unsigned int)currentUncert];
        shift = -1.;
    }

    currentP4.SetPtEtaPhiM( p4.Pt()*(1 + scale*shift),
                            p4.Eta(),
                            p4.Phi(),
                            p4.M()*(1 + scale*shift)
                          );
    return currentP4;

}
////////////////////////////////////////////////
////////////////////////////////////////////////
double HTTJet::getMaxPt()
{
    double maxup = *std::max_element( jecUncertSourceValuesUp.begin(),jecUncertSourceValuesUp.end()  );
    double maxdown = *std::max_element( jecUncertSourceValuesDown.begin(),jecUncertSourceValuesDown.end()  ) ;

    return p4.Pt()*( 1 + std::max(maxup, maxdown*(-1) ) );
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJet::setJecUncertSourceValue(unsigned int uncert, double value, bool up)
{
    if(up) jecUncertSourceValuesUp[uncert] = value;
    else   jecUncertSourceValuesDown[uncert] = value;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HTTJetCollection::clear()
{
    currentUncert = JecUncertEnum::NONE;
    currentShift = true;

    dijet.SetPtEtaPhiM(DEF,DEF,DEF,DEF);
    jetCollection.clear();
    jetCurrentCollection.clear();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::setCurrentUncertShift(JecUncertEnum uncert, bool up)
{
    //Make sure that everything gets shifted here
    currentUncert = uncert;
    currentShift = up;
    fillCurrentCollection();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::addJet(HTTJet newJet)
{
    jetCollection.push_back(newJet);
    if(newJet.getP4().Pt() > 20 )
        jetCurrentCollection.push_back(newJet);

}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::fillCurrentCollection()
{
    jetCurrentCollection.clear();
    for(auto jet : jetCollection)
    {
        jet.setUncertShift(currentUncert, currentShift);
        if(jet.getP4().Pt() > 20)
            jetCurrentCollection.push_back(jet);
    }
    if(jetCurrentCollection.size() > 1)
        setDijetP4();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
void HTTJetCollection::setDijetP4()
{
    dijet = jetCurrentCollection.at(0).getP4()+jetCurrentCollection.at(1).getP4();
}
////////////////////////////////////////////////
////////////////////////////////////////////////
int HTTJetCollection::getNJets(double pt)
{
    int njets = 0;
    for(auto &jet : jetCurrentCollection)
        if(jet.getP4().Pt() > pt )njets++;
    
    return njets;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
int HTTJetCollection::getNJetInGap(double pt)
{
    int njetingap = 0;
    if(jetCurrentCollection.size() < 2) return 0;

    double j1_eta =  jetCurrentCollection.at(0).getP4().Eta();
    double j2_eta =  jetCurrentCollection.at(1).getP4().Eta();

    for (unsigned ij=2; ij<jetCurrentCollection.size(); ij++)
    {
        TLorentzVector addJet = jetCurrentCollection.at(ij).getP4();

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
int HTTJetCollection::getNBtag()
{
    int nbtags = 0;
    for(auto &jet : jetCurrentCollection)
        if( jet.isBtagJet() )nbtags++;

    return nbtags;
}

