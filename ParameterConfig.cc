struct Parameter{
  float defaultValue = -10.;
  struct Datasets{
    struct DY{

      struct Incl{
        float Xs = 4954 ;
        float NNLOXs = 5765.4 ;
        float Nevt = 145803217;
      } Incl;
      struct High{
        float Xs = 6.657  ;
        float Nevt = 6108651;
      } High;
      struct Low{
        float Xs = 18610 ;
        float Nevt = 35291566;
      } Low;
      struct Jet1{
        float Xs = 1012.5 ;
        float Nevt = 62627174;
      } Jet1;
      struct Jet2{
        float Xs =  332.8;
        float Nevt = 19970551;
      } Jet2;
      struct Jet3{
        float Xs = 101.8;
        float Nevt = 5856110;
      } Jet3;
      struct Jet4{
        float Xs = 54.8;
        float Nevt = 4197868;
      } Jet4;
      struct M150{
        float Xs = 6.657;
        float Nevt = 4099481;
      } M150;
    } DY;

    struct W{
      struct Incl{
        float Xs =  50380;
        float NNLOXs = 61526.7;
        float Nevt = 86731806;
      } Incl;
      struct Jet1{
        float Xs = 9644.5;
        float Nevt = 45367044;
      } Jet1;
      struct Jet2{
        float Xs = 3144.5;
        float Nevt = 60197766;
      } Jet2;
      struct Jet3{
        float Xs = 954.8;
        float Nevt = 59067548;
      } Jet3;
      struct Jet4{
        float Xs = 485.6;
        float Nevt = 29995313;
      } Jet4;
    } W;

  } Datasets;


///////////// Electron //////////////
  struct Electron{

    struct Baseline{
     float pt = 26.;
     float eta = 2.1;
     float dxy = 0.045;
     float dz = 0.2;
     float overlap = 0.5;
    } Baseline;

    struct TES{
      float one_prong_0p0   =  0.000;
      float one_prong_1p0   =  0.095;
      float three_prong_0p0 =  0.000;
      float uncert          =  0.03;
    }TES;

    struct DiElectronVeto{
      float pt = 15. ;
      float eta = 2.5 ;
      float dxy = 0.045 ;
      float dz = 0.2 ;
      float iso = 0.3;
      float overlap = 0.15;
    } DiElectronVeto;

    struct ThirdLeptonVeto{
      float pt = 10.;
      float eta = 2.5;
      float dxy = 0.045;
      float dz = 0.2;
      float iso = 0.3;
    } ThirdLeptonVeto;

    struct AdditionalElectrons{
      float pt = 13.;
      float eta = 2.1;
      float dxy = 0.045;
      float dz = 0.2 ;
      float iso = 0.3;
    } AdditionalElectrons;

    struct MVAId{
      bool WP90(float Eta, float Pt, float Id){
        return (( fabs(Eta)<0.8
                  && Pt > 10.
                  && Id < 0.837
                )
                ||( fabs(Eta)>=0.8 && fabs(Eta)<1.479
                    && Pt > 10.
                    && Id < 0.715
                )
                ||( fabs(Eta)>=1.479
                    && Id < 0.357
                )
               );
      }
     bool WP80(float Eta, float Pt, float Id){
        return (( fabs(Eta)<0.8
                   && Pt >10.  
                   && Id < 0.941
                )
                ||( fabs(Eta)>=0.8 && fabs(Eta)<1.479
                    && Pt > 10
                    && Id < 0.899
                ) 
                ||( fabs(Eta)>=1.479
                    && Id < 0.758
                )

               );
      }
    } MVAId;

  } Electron;

/////////////////////////////////  
///////////// Muon //////////////
  struct Muon{

    struct Baseline{
     float pt = 20.;
     float trailPt = 10;
     //float eta = 2.4;
     float eta = 2.1;
     float dxy = 0.045;
     float dz = 0.2;
     float overlap = 0.5;
     float trailOverlap = 0.3;
    } Baseline;

    struct TES{
      float one_prong_0p0   = -0.002;
      float one_prong_1p0   =  0.015;
      float three_prong_0p0 =  0.000;
      float uncert          =  0.015;
    }TES;

    struct DiMuonVeto{
      float pt = 15. ;
      float eta = 2.4 ;
      float dxy = 0.045 ;
      float dz = 0.2 ;
      float iso = 0.3;
      float overlap = 0.15;
    } DiMuonVeto;

    struct ThirdLeptonVeto{
      float pt = 10.;
      float eta = 2.4;
      float dxy = 0.045;
      float dz = 0.2;
      float iso = 0.3;
    } ThirdLeptonVeto;

    struct AdditionalMuon{
      float pt = 5.;
      float eta = 2.1;
      float dxy = 0.045;
      float dz = 0.2;
      float iso = 0.3;

    } AdditionalMuon;

  } Muon;

/////////////////////////////////  
///////////// Tau //////////////
  struct Tau{

    struct TES{
      float one_prong_0p0   = -0.018;
      float one_prong_1p0   =  0.010;
      float three_prong_0p0 =  0.004;
      float uncert          =  0.012;
    }TES;


    struct Baseline{

      float overlap = 0.5;

      float leadPt = 50;
      float pt(string channel){
        if(channel == "et") return 30.;
        if(channel == "mt") return 30.;
        if(channel == "tt") return 40.;

        return 0.;
      }

      float eta(string channel){
        if(channel == "et") return 2.3;
        if(channel == "mt") return 2.3;
        if(channel == "tt") return 2.1;

        return 0.;
      }

      float IdDecay(string channel){
        if(channel == "et") return 0.5;
        if(channel == "mt") return 0.5;
        if(channel == "tt") return 0.5;

        return 0.;
      }

      float LeadDz(string channel){
        if(channel == "et") return 0.2;
        if(channel == "mt") return 0.2;
        if(channel == "tt") return 0.2;

        return 0.;
      }
    }Baseline;

  } Tau;

/////////////////////////////////  
///////////// Jets //////////////

  struct Jets{

    float dR = 0.5;

    struct Jet{
      float pt = 20.;
      float eta = 4.7;
    } Jet;

    struct BJet{
      float pt = 20.;
      float eta = 2.4;
      float BCSV = 0.8484;
    } BJet;

  } Jets;

} Parameter;
