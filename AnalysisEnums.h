#ifndef RootAnalysis_AnalysisEnums_H
#define RootAnalysis_AnalysisEnums_H

#include <map>
#include <vector>
#include <iostream>
#include "EnergyScales.h"

namespace HTTAnalysis {

///Copy from DataFormats/TauReco/interface/PFTauDecayMode.h
enum hadronicTauDecayModes
{
    tauDecay1ChargedPion0PiZero,
    tauDecay1ChargedPion1PiZero, // rho (770 MeV) mediated)
    tauDecay1ChargedPion2PiZero, // a1  (1.2 GeV) mediated
    tauDecay1ChargedPion3PiZero, // contaminated or unmerged photo
    tauDecay1ChargedPion4PiZero, // contaminated or unmerged photo
    tauDecay2ChargedPion0PiZero, // extra track or un-recod track
    tauDecay2ChargedPion1PiZero, // extra track or un-recod track
    tauDecay2ChargedPion2PiZero, // extra track or un-recod track
    tauDecay2ChargedPion3PiZero, // extra track or un-recod track
    tauDecay2ChargedPion4PiZero, // extra track or un-recod track
    tauDecay3ChargedPion0PiZero, // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion1PiZero, // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion2PiZero, // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion3PiZero, // a1  (1.2 GeV) mediated
    tauDecay3ChargedPion4PiZero, // a1  (1.2 GeV) mediated
    tauDecaysElectron,
    tauDecayMuon,
    tauDecayOther             // catch-all
};

class eventCategory {
public:

eventCategory(const std::string & aName, std::vector<const eventCategory*> & categoryRejester){
    myName = aName;

    myId = categoryRejester.size();
    categoryRejester.push_back(this);

    myWSF = 0;
    myQCDEstimate = 0;
    myQCDSFNum = 0;
    myQCDSFDenom = 0;
    if(aName.find("_W")==std::string::npos) myWSF = new eventCategory(aName+"_W",categoryRejester);
    if(aName.find("_QCD")==std::string::npos) {
        myQCDEstimate = new eventCategory(aName+"_QCD",categoryRejester);
        myQCDSFNum = new eventCategory(aName+"_QCDSFNum",categoryRejester);
        myQCDSFDenom = new eventCategory(aName+"_QCDSFDenom",categoryRejester);
    }
};

unsigned int id() const {
    return myId;
};

const eventCategory * wSF() const {
    if(myWSF) return myWSF;
    else return this;
};
const eventCategory * qcdEstimate() const {
    if(myQCDEstimate) return myQCDEstimate;
    else return this;
};
const eventCategory * qcdSFNumerator() const {
    if(myQCDSFNum) return myQCDSFNum;
    else return this;
};
const eventCategory * qcdSFDenominator() const {
    if(myQCDSFDenom) return myQCDSFDenom;
    else return this;
};

const std::string & name() const {
    return myName;
}

private:

std::string myName;
unsigned int myId;
eventCategory *myWSF;
eventCategory *myQCDEstimate;
eventCategory *myQCDSFNum;
eventCategory *myQCDSFDenom;


};

enum sysEffects 
{   NOMINAL,
    TES1p0p0Up, TES1p0p0Down,
    TES1p1p0Up, TES1p1p0Down,
    TES3p0p0Up, TES3p0p0Down,

    MES1p0p0Up, MES1p0p0Down,
    MES1p1p0Up, MES1p1p0Down,
    MES3p0p0Up, MES3p0p0Down,

    EES1p0p0Up, EES1p0p0Down,
    EES1p1p0Up, EES1p1p0Down,
    EES3p0p0Up, EES3p0p0Down,
    EESUp,      EESDown,
    DUMMY_SYS,
///Place systematic effects not affecting the SV calculation after DUMMY_SYS
///all quantities for following syst effects are calculated on fly, no need to rerun
///the ntuple making step.
    JESUp, JESDown,
    J2TUp, J2TDown,
    ZPtUp, ZPtDown,
    TTUp, TTDown,
    QCDSFUp, QCDSFDown,
    WSFUp, WSFDown,
    ggUp, ggDown,
    ZmumuUp, ZmumuDown
};

enum finalState
{
  EleTau = 0,
  MuTau = 1,
  TauTau = 2,
  NONE = 4
};

float getEnergyScale( int pdg, int mc_match, float eta, int dm, sysEffects type = sysEffects::NOMINAL )
{
    float shift = 0.0;
    float scale = 0.0;
    if(type > sysEffects::NOMINAL && type < sysEffects::DUMMY_SYS){
        shift = type % 2 == 0 ? -1.0 : 1.0;
    }
    if(pdg == 11)
    {
        if(mc_match == 1 || mc_match == 3)
        {
              if(type == sysEffects::EESUp || type == sysEffects::EESDown){
                if(fabs(eta) < 1.479 ) return shift*ES.Electron.uncertaintyBarrel;
                return shift*ES.Electron.uncertaintyEndcap;
              }
        }      
    }

    if(pdg == 15)
    {
        if(mc_match == 1 || mc_match == 3)
        {
            if(dm == 0)
            {
              scale = ES.Electron.oneProng0p0;
              if(type == sysEffects::EES1p0p0Up || type == sysEffects::EES1p0p0Down) scale += shift*ES.Electron.uncertaintyShift1p;
              return scale;
            }
            if(dm == 1)
            {
              scale = ES.Electron.oneProng1p0;
              if(type == sysEffects::EES1p1p0Up || type == sysEffects::EES1p1p0Down) scale += shift*ES.Electron.uncertaintyShift1p1p;
              return scale;
            }
            if(dm == 10)
            {
              scale = ES.Electron.threeProng0p0;
              if(type == sysEffects::EES3p0p0Up || type == sysEffects::EES3p0p0Down) scale += shift*ES.Electron.uncertaintyShift3p;
              return scale;
            }
        }
        if(mc_match == 2 || mc_match == 4)
        {
            if(dm == 0)
            {
              scale = ES.Muon.oneProng0p0;
              if(type == sysEffects::MES1p0p0Up || type == sysEffects::MES1p0p0Down) scale += shift*ES.Muon.uncertaintyShift1p;
              return scale;
            }
            if(dm == 1)
            {
              scale = ES.Muon.oneProng1p0;
              if(type == sysEffects::MES1p1p0Up || type == sysEffects::MES1p1p0Down) scale += shift*ES.Muon.uncertaintyShift1p1p;
              return scale;
            }
            if(dm == 10)
            {
              scale = ES.Muon.threeProng0p0;
              if(type == sysEffects::MES3p0p0Up || type == sysEffects::MES3p0p0Down) scale += shift*ES.Muon.uncertaintyShift3p;
              return scale;
            }
        }
        if(mc_match == 5)
        {
            if(dm == 0)
            {
              scale = ES.Tau.oneProng0p0;
              if(type == sysEffects::TES1p0p0Up || type == sysEffects::TES1p0p0Down) scale += shift*ES.Tau.uncertaintyShift1p;
              return scale;
            }
            if(dm == 1)
            {
              scale = ES.Tau.oneProng1p0;
              if(type == sysEffects::TES1p1p0Up || type == sysEffects::TES1p1p0Down) scale += shift*ES.Tau.uncertaintyShift1p1p;
              return scale;
            }
            if(dm == 10)
            {
              scale = ES.Tau.threeProng0p0;
              if(type == sysEffects::TES3p0p0Up || type == sysEffects::TES3p0p0Down) scale += shift*ES.Tau.uncertaintyShift3p;
              return scale;
            }
        }

    }
    return scale;
}

}
#endif
