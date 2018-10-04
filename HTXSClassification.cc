#include "HTXSClassification.h"
   
/// @brief Return bin index of x given the provided bin edges. 0=first bin, -1=underflow bin.
int getBin(double x, std::vector<double> bins)
{
  if (bins.size()==0||x<bins[0]) return -1; // should not happen!
  for (size_t i=1;i<bins.size();++i)
    if (x<bins[i]) return i-1;
  return bins.size()-1;
}

/// @brief VBF topolog selection
/// 0 = fail loose selction: m_jj > 400 GeV and Dy_jj > 2.8
/// 1 pass loose, but fail additional cut pT(Hjj)<25. 2 pass tight selection
int vbfTopology(HTTJetCollection *jets, const TLorentzVector &higgs)
{
  if (jets->getNJets(30)<2) return 0;
  TLorentzVector j1, j2;
  j1 = jets->getJet(0).P4();
  j2 = jets->getJet(1).P4();

  bool VBFtopo = (j1+j2).M() > 400.0 && std::abs(j1.Eta()-j2.Eta()) > 2.8;
  return VBFtopo ? (j1+j2+higgs).Pt()<25 ? 2 : 1 : 0;
}

/// @brief Stage-1 categorization
Stage1::Category getStage1Category( HiggsProdMode prodMode,
                                    const TLorentzVector &higgs,
                                    HTTJetCollection *jets)
{
  using namespace Stage1;
  int Njets=jets->getNJets(30);
  double pTj1 = Njets ? jets->getJet(0).Pt() : 0;
  int vbfTopo = vbfTopology(jets, higgs);

  // 1. GGF Stage 1 categories
  //    Following YR4 write-up: XXXXX
  if (prodMode==HiggsProdMode::GGF ) {
    if (Njets==0)        return Category::GG2H_0J;
    else if (Njets==1)   return Category(Category::GG2H_1J_PTH_0_60+getBin(higgs.Pt(),{0,60,120,200}));
    else if (Njets>=2) {
      // events with pT_H>200 get priority over VBF cuts
      if(higgs.Pt()<=200){
        if      (vbfTopo==2) return Category::GG2H_VBFTOPO_JET3VETO;
        else if (vbfTopo==1) return Category::GG2H_VBFTOPO_JET3;
      }
      // Njets >= 2jets without VBF topology
      return Category(Category::GG2H_GE2J_PTH_0_60+getBin(higgs.Pt(),{0,60,120,200}));
    }
  }
  // 2. Electroweak qq->Hqq Stage 1 categories
  else if (prodMode==HiggsProdMode::VBF  ) {
    if (pTj1>200) return Category::QQ2HQQ_PTJET1_GT200;
    if (vbfTopo==2) return Category::QQ2HQQ_VBFTOPO_JET3VETO;
    if (vbfTopo==1) return Category::QQ2HQQ_VBFTOPO_JET3;
    double mjj = Njets>1 ? (jets->getJet(0).P4()+jets->getJet(1).P4()).M():0;
    if ( 60 < mjj && mjj < 120 ) return Category::QQ2HQQ_VH2JET;
    return Category::QQ2HQQ_REST;
  }
  return Category::UNKNOWN;
}