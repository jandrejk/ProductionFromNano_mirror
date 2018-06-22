#ifdef PROJECT_NAME
#include "m2n/HTTDataFormats/interface/HTTEvent.h"
#include "m2n/HTT/interface/GenInfoHelper.h"
#else
#include "HTTEvent.h"
#include "EnergyScales.h"
#endif

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

  mcWeight = 1.0;
  ptReWeight = 1.0;
  ptReWeightSUSY = 1.0;
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
////////////////////////////////////////////////
////////////////////////////////////////////////
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
const TLorentzVector & HTTParticle::getNominalShiftedP4() const
{

    p4Cache = p4;

    if( std::abs(getPDGid())==15 )
    {
        int dm = getProperty(PropertyEnum::decayMode);
        float tauES = 0.0;

        if( getProperty(PropertyEnum::mc_match)==5 )
        {
            if(dm == 0 ) tauES = ES.Tau.oneProng0p0;
            if(dm == 1 ) tauES = ES.Tau.oneProng1p0;
            if(dm == 10) tauES = ES.Tau.threeProng0p0;

        }else if( getProperty(PropertyEnum::mc_match)==1 )
        {
            if(dm == 0 ) tauES = ES.Electron.oneProng0p0;
            if(dm == 1 ) tauES = ES.Electron.oneProng1p0;
            if(dm == 10) tauES = ES.Electron.threeProng0p0;

        }else if( getProperty(PropertyEnum::mc_match)==2 )
        {
            if(dm == 0 ) tauES = ES.Muon.oneProng0p0;
            if(dm == 1 ) tauES = ES.Muon.oneProng1p0;
            if(dm == 10) tauES = ES.Muon.threeProng0p0;
        }

        float tauES_mass = tauES;
        if (dm == 0) tauES_mass=0;

        p4Cache.SetPtEtaPhiM(p4.Pt() * (1.0+tauES),
                             p4.Eta(),
                             p4.Phi(),
                             p4.M() * (1.0+tauES_mass) );

    }

    return p4Cache;
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

    }else if(type==HTTAnalysis::NOMINAL){
        lastSystEffect = type;
        return getNominalShiftedP4();

    }else if( (unsigned int)type < (unsigned int)HTTAnalysis::DUMMY_SYS)
    {
        lastSystEffect = type;

        p4Cache = p4;

        if( std::abs(getPDGid())==15 )
        {
            int dm = getProperty(PropertyEnum::decayMode);
            float tauES = 0.0;

            if( getProperty(PropertyEnum::mc_match)==1 ) // Prompt Electrons
            {
                if(dm == 0 ) tauES = getShiftedES( ES.Electron.oneProng0p0,   ES.Electron.uncertaintyShift ,dm,type );
                if(dm == 1 ) tauES = getShiftedES( ES.Electron.oneProng1p0,   ES.Electron.uncertaintyShift ,dm,type );
                if(dm == 10) tauES = getShiftedES( ES.Electron.threeProng0p0, ES.Electron.uncertaintyShift ,dm,type );

            }else if( getProperty(PropertyEnum::mc_match)==2 ) // Prompt Muon
            {
                if(dm == 0 ) tauES = getShiftedES( ES.Muon.oneProng0p0,   ES.Muon.uncertaintyShift ,dm,type );
                if(dm == 1 ) tauES = getShiftedES( ES.Muon.oneProng1p0,   ES.Muon.uncertaintyShift ,dm,type );
                if(dm == 10) tauES = getShiftedES( ES.Muon.threeProng0p0, ES.Muon.uncertaintyShift ,dm,type );

            }else if( getProperty(PropertyEnum::mc_match)==5 ) // Genuine Tau
            {
                if(dm == 0 ) tauES = getShiftedES( ES.Tau.oneProng0p0,   ES.Tau.uncertaintyShift ,dm,type );
                if(dm == 1 ) tauES = getShiftedES( ES.Tau.oneProng1p0,   ES.Tau.uncertaintyShift ,dm,type );
                if(dm == 10) tauES = getShiftedES( ES.Tau.threeProng0p0, ES.Tau.uncertaintyShift ,dm,type );
            }

            float tauES_mass = tauES;
            if (dm == 0) tauES_mass=0;

            p4Cache.SetPtEtaPhiM(p4.Pt() * (1.0+tauES),
                                 p4.Eta(),
                                 p4.Phi(),
                                 p4.M() * (1.0+tauES_mass) );

        }

        return p4Cache;

    }


    p4Cache = p4;

    return p4;
}
////////////////////////////////////////////////
////////////////////////////////////////////////
float HTTParticle::getShiftedES( float ES, float uncert, int dm, HTTAnalysis::sysEffects type) const
{

    float shift = (unsigned int)type % 2 == 0 ? -1.0 : 1.0;

    if(dm == 0)
    {
        if(type == HTTAnalysis::TES0p0Up || type == HTTAnalysis::TES0p0Down
           || type == HTTAnalysis::MES0p0Up || type == HTTAnalysis::MES0p0Down
           || type == HTTAnalysis::EES0p0Up || type == HTTAnalysis::EES0p0Down)
        {
            return ES + shift*uncert;
        }
    }
    if(dm == 1)
    {
        if(type == HTTAnalysis::TES1p0Up || type == HTTAnalysis::TES1p0Down
           || type == HTTAnalysis::MES1p0Up || type == HTTAnalysis::MES1p0Down
           || type == HTTAnalysis::EES1p0Up || type == HTTAnalysis::EES1p0Down)
        {
            return ES + shift*uncert;
        }
    }
    if(dm == 10)
    {
        if(type == HTTAnalysis::TES3p0Up || type == HTTAnalysis::TES3p0Down
           || type == HTTAnalysis::MES3p0Up || type == HTTAnalysis::MES3p0Down
           || type == HTTAnalysis::EES3p0Up || type == HTTAnalysis::EES3p0Down)
        {
            return ES + shift*uncert;
        }
    }

    return ES;

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
////////////////////////////////////////////////
////////////////////////////////////////////////
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
        if( (std::abs(leg1.getPDGid())==15 && leg1.getProperty(PropertyEnum::mc_match)==5) 
            || (std::abs(leg2.getPDGid())==15 && leg2.getProperty(PropertyEnum::mc_match)==5) )
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
////////////////////////////////////////////////
////////////////////////////////////////////////
