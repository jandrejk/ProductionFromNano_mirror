#define HTauTauTreeFromNanoBase_cxx

#include "HTauTauTreeFromNanoBase.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLeaf.h>

#include <iostream>
#include <fstream>
#include <algorithm>

//move these two to the configuration

HTauTauTreeFromNanoBase::HTauTauTreeFromNanoBase(TTree *tree, std::vector<edm::LuminosityBlockRange> lumiBlocks, std::string prefix) : NanoEventsSkeleton(tree)
{

    ifstream S("configBall.json");
    Settings = json::parse(S);

    isMC = Settings["isMC"].get<bool>();
    isSync = Settings["isSync"].get<bool>();

    httJetCollection.initCollection(isMC, isSync);

    ///Init HTT ntuple
    initHTTTree(tree, prefix);

    jsonVector = lumiBlocks;

    ///Initialization of SvFit
    if( Settings["svfit"].get<bool>() )
    {
        std::cout<<"[HTauTauTreeFromNanoBase]: Run w/ SVFit"<<std::endl;
        unsigned int verbosity = 0;//Set the debug level to 3 for testing
        svFitAlgo_ = std::unique_ptr<ClassicSVfit>(new ClassicSVfit(verbosity) );
    } else
    {
        std::cout<<"[HTauTauTreeFromNanoBase]: Run w/o SVFit"<<std::endl;
        svFitAlgo_=nullptr;
    }

    ///Initialization of RecoilCorrector
    if( applyRecoil )
    {
        std::cout<<"[HTauTauTreeFromNanoBase]: Apply MET recoil corrections"<<std::endl;
        std::string correctionFile = "HTT-utilities/RecoilCorrections/data/Type1_PFMET_2017.root";
        recoilCorrector_= std::unique_ptr<RecoilCorrector>( new RecoilCorrector(correctionFile) );
        metSys_         = std::unique_ptr<MEtSys>( new MEtSys("HTT-utilities/RecoilCorrections/data/MEtSys.root") );

    } else
    {
        std::cout<<"[HTauTauTreeFromNanoBase]: Do not apply MET recoil corrections"<<std::endl;
        recoilCorrector_=nullptr;
    }

    ///Get files with weights
    zPtReweightFile = std::unique_ptr<TFile>( new TFile("utils/zptweight/zpt_weights_2017.root") );  
    if(!zPtReweightFile) std::cout<<"Z pt reweight file zpt_weights.root is missing."<<std::endl;
    zptmass_histo = (TH2D*)zPtReweightFile->Get("zptmass_histo");

    puweights = std::unique_ptr<TFile>( new TFile("utils/puweight/puweights.root") );  
    if(!puweights) std::cout<<"puweights.root is missing."<<std::endl;
    // puweights_histo = (TH1D*)puweights->Get("#VBFHToTauTau_M125_13TeV_powheg_pythia8#RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1#MINIAODSIM");   
    puweights_histo = (TH1D*)puweights->Get( Settings["puTag"].get<string>().c_str() );   


    if(isMC)
    {
        ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
        std::cout<<"[HTauTauTreeFromNanoBase]: Instantiate JEC uncertainty sources"<<std::endl;
        initJecUnc("utils/jec_uncert/Fall17_17Nov2017_V6_MC_UncertaintySources_AK4PFchs.txt");

        std::cout<<"[HTauTauTreeFromNanoBase]: Load files and init for promote-demote"<<std::endl;
        httJetCollection.initForPromoteDemote();
    }

    if(httEvent->getSampleType() == HTTEvent::h)
    {
        nnlo_ggh_graphs = std::unique_ptr<TFile>( new TFile("utils/NNLO_ggH/NNLOPS_reweight.root") );
        NNLOPSratio_pt_mcatnlo_0jet = (TGraphErrors*)nnlo_ggh_graphs->Get("gr_NNLOPSratio_pt_mcatnlo_0jet");
        NNLOPSratio_pt_mcatnlo_1jet = (TGraphErrors*)nnlo_ggh_graphs->Get("gr_NNLOPSratio_pt_mcatnlo_1jet");
        NNLOPSratio_pt_mcatnlo_2jet = (TGraphErrors*)nnlo_ggh_graphs->Get("gr_NNLOPSratio_pt_mcatnlo_2jet");
        NNLOPSratio_pt_mcatnlo_3jet = (TGraphErrors*)nnlo_ggh_graphs->Get("gr_NNLOPSratio_pt_mcatnlo_3jet");
    }
    else
    {
        nnlo_ggh_graphs=nullptr;
    }
    

    firstWarningOccurence_=true;

}

HTauTauTreeFromNanoBase::~HTauTauTreeFromNanoBase()
{
  httFile->Write();
  cout << "Finished!" << endl;
}

/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::initHTTTree(const TTree *tree, std::string prefix)
{

    if(prefix=="") prefix="HTT";
    prefix += "_";

    std::string filePath(tree->GetCurrentFile()->GetName());
    size_t location = filePath.find_last_of("/");
    if(location==std::string::npos) location = 0;
    else location+=1;
    std::string fileName = prefix+filePath.substr(location,filePath.size());
    httFile = std::unique_ptr<TFile>( new TFile(fileName.c_str(),"RECREATE") );
    httEvent = std::unique_ptr<HTTEvent>(new HTTEvent() );
    httEvent->setSampleType( Settings["sample"].get<string>() );

    metShifts.clear();
    if( httEvent->getSampleType() == HTTEvent::TTbar 
        || httEvent->getSampleType() == HTTEvent::ST
        || httEvent->getSampleType() == HTTEvent::Diboson
        || !isMC )
    {
        applyRecoil = false;

    }else
    {
        applyRecoil = Settings["recoil"].get<bool>() ;       
    }

    if(applyRecoil)
    {
        // Dont judge me... Not very nice way to implement met uncerts
        // TODO: Put it in a smart container that calculates shifts and has knowledge of shift names
        metShifts.push_back( make_pair("",                                       make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Up ) ) ); // Dummy shifts for nominal
        metShifts.push_back( make_pair("CMS_htt_boson_reso_met_13TeVUp",         make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Up ) ) );
        metShifts.push_back( make_pair("CMS_htt_boson_reso_met_13TeVDown",       make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Down ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_reso_met_0Jet_13TeVUp",    make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Up ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_reso_met_0Jet_13TeVDown",  make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Down ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_reso_met_1Jet_13TeVUp",    make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Up ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_reso_met_1Jet_13TeVDown",  make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Down ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_reso_met_2Jet_13TeVUp",    make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Up ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_reso_met_2Jet_13TeVDown",  make_pair(MEtSys::SysType::Response,   MEtSys::SysShift::Down ) ) );
        metShifts.push_back( make_pair("CMS_htt_boson_scale_met_13TeVUp",        make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Up ) ) );
        metShifts.push_back( make_pair("CMS_htt_boson_scale_met_13TeVDown",      make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Down ) ) );        
        // metShifts.push_back( make_pair("CMS_htt_boson_scale_met_0Jet_13TeVUp",   make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Up ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_scale_met_0Jet_13TeVDown", make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Down ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_scale_met_1Jet_13TeVUp",   make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Up ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_scale_met_1Jet_13TeVDown", make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Down ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_scale_met_2Jet_13TeVUp",   make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Up ) ) );
        // metShifts.push_back( make_pair("CMS_htt_boson_scale_met_2Jet_13TeVDown", make_pair(MEtSys::SysType::Resolution, MEtSys::SysShift::Down ) ) ); 
    }
    

    //  httTree = new TTree("HTauTauTree","");
    //  httTree->SetDirectory(httFile);

    hStats = new TH1F("hStats","Bookkeeping histogram",11,-0.5,10.5);
    //  hStats->SetDirectory(httFile);

    t_TauCheck=new TTree("TauCheck","TauCheck");
    evtWriter = std::unique_ptr<EventWriter>( new EventWriter() );
    evtWriter->initTree(t_TauCheck, httJetCollection.getNeededJECShifts(), isMC, isSync, metShifts);
    
    leptonPropertiesList = leptonProperties; // Defined in PropertyEnum.h

    ////////////////////////////////////////////////////////////
    HTTEvent::usePropertyFor["electronIsolation"]  = PropertyEnum::pfRelIso03_all;
    HTTEvent::usePropertyFor["electronIDWP80"]     = PropertyEnum::mvaFall17noIso_WP80_v2;
    HTTEvent::usePropertyFor["electronIDWP90"]     = PropertyEnum::mvaFall17noIso_WP90_v2;
    HTTEvent::usePropertyFor["electronIDCutBased"] = PropertyEnum::cutBased_v2;
    HTTEvent::usePropertyFor["muonIsolation"]      = PropertyEnum::pfRelIso04_all;
    HTTEvent::usePropertyFor["muonID"]             = PropertyEnum::mediumId;
    HTTEvent::usePropertyFor["tauIsolation"]       = PropertyEnum::rawMVAoldDM2017v2;
    HTTEvent::usePropertyFor["tauID"]              = PropertyEnum::idMVAoldDM2017v2;

    ///Trigger bits to check
    /// Must be aligned with TriggerEnum.h
    ///FIXME: is there a nicer way to define trigger list, e.g. a cfg file?
    TriggerData aTrgData;

    // 2017 94X  Filter

    // Electron
    // 0 CaloIdL_TrackIdL_IsoVL                CaloIdLTrackIdLIsoVL*TrackIso*Filter
    // 1 WPTight                               hltEle*WPTight*TrackIsoFilter*
    // 2 WPLoose                               hltEle*WPLoose*TrackIsoFilter
    // 3 OverlapFilter PFTau                   *OverlapFilterIsoEle*PFTau*

    // Muon
    // 0 = TrkIsoVVL                           RelTrkIsoVVLFiltered0p4
    // 1 = Iso                                 hltL3crIso*Filtered0p07
    // 2 = OverlapFilter PFTau                 *OverlapFilterIsoMu*PFTau*

    // Tau
    // 0 = LooseChargedIso                     LooseChargedIso*
    // 1 = MediumChargedIso                    *MediumChargedIso*
    // 2 = TightChargedIso                     *TightChargedIso*
    // 3 = TightID OOSC photons                *TightOOSCPhotons*
    // 4 = L2p5 pixel iso                      hltL2TauIsoFilter
    // 5 = OverlapFilter IsoMu                 *OverlapFilterIsoMu*
    // 6 = OverlapFilter IsoEle                *OverlapFilterIsoEle*
    // 7 = L1-HLT matched                      *L1HLTMatched*
    // 8 = Dz                                  *Dz02*


    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu24";
    triggerBits_.back().leg1Id=13;
    triggerBits_.back().leg1BitMask=(1<<3);
    // triggerBits_.back().leg1Pt=24;
    triggerBits_.back().leg1L1Pt=22;
    // triggerBits_.back().leg1OfflinePt=25;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu27";
    triggerBits_.back().leg1Id=13;
    triggerBits_.back().leg1BitMask=(1<<3);
    // triggerBits_.back().leg1Pt=27;
    triggerBits_.back().leg1L1Pt=22; //22 or 25...
    // triggerBits_.back().leg1OfflinePt=28;

    // Mu tauh triggers
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1";
    triggerBits_.back().leg1Id=13;
    triggerBits_.back().leg1BitMask=(1<<1) + (1<<2); //iso+OL
    triggerBits_.back().leg1Pt=20;
    triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg1L1Pt=-1;
    //  triggerBits_.back().leg1L1Pt=18;
    triggerBits_.back().leg1OfflinePt=20;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<0) + (1<<8); //looseChargedIso+OL
    triggerBits_.back().leg2Pt=27;
    triggerBits_.back().leg2Eta=2.1;
    //  triggerBits_.back().leg2L1Pt=20;
    triggerBits_.back().leg2L1Pt=-1;

    //Single e triggers
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele27_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1);
    triggerBits_.back().leg1Pt=27;
    triggerBits_.back().leg1L1Pt=-1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele32_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1);
    triggerBits_.back().leg1Pt=32;
    triggerBits_.back().leg1L1Pt=-1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele35_WPTight_Gsf";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1);
    triggerBits_.back().leg1Pt=35;
    triggerBits_.back().leg1L1Pt=-1;


    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1";
    triggerBits_.back().leg1Id=11;
    triggerBits_.back().leg1BitMask=(1<<1)+(1<<3);
    triggerBits_.back().leg1Pt=24;
    triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<0) + (1<<7); 
    triggerBits_.back().leg2Pt=30;
    triggerBits_.back().leg2Eta=2.1;


    // Single tauh triggers
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=2;

    // tauh tauh triggers
    ///9th tau bit(1<<8) for di-tau dz filter  (should be OK 80X triggers)
    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask= (1<<2) + (1<<3) + (1<<6); //TightChargedIso+photons+dz
    triggerBits_.back().leg1Pt=35;
    triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask= (1<<2) + (1<<3) + (1<<6);
    triggerBits_.back().leg2Pt=35;
    triggerBits_.back().leg2Eta=2.1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=(1<<1) + (1<<3) + (1<<6); //MediumChargedIso+photons+dz
    triggerBits_.back().leg1Pt=40;
    triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<1) + (1<<3) + (1<<6);
    triggerBits_.back().leg2Pt=40;
    triggerBits_.back().leg2Eta=2.1;

    triggerBits_.push_back(aTrgData);
    triggerBits_.back().path_name="HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg";
    triggerBits_.back().leg1Id=15;
    triggerBits_.back().leg1BitMask=(1<<2) + (1<<6);; //TightChargedIso+dz
    triggerBits_.back().leg1Pt=40;
    triggerBits_.back().leg1Eta=2.1;
    triggerBits_.back().leg2Id=15;
    triggerBits_.back().leg2BitMask=(1<<2) + (1<<6);;
    triggerBits_.back().leg2Pt=40;
    triggerBits_.back().leg2Eta=2.1;
    ////////////////////////////////////////////////////////////
    ///Filter bits to check
    for(auto filter : FilterNames) // Defined in FilterEnum.h
        filterBits_.push_back(filter);
    
    return;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::Loop(Long64_t nentries_max, unsigned int sync_event)
{

    check_event_number = sync_event;

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();
    Long64_t nentries_use=nentries;
    if (nentries_max>0 && nentries_max < nentries)
    {
        nentries_use=nentries_max;
        check_event_number = 2038471; // Usefull event in mutau for debug in 2017
    }

    Long64_t nbytes = 0, nb = 0;
    int entry=0;
    float perc = 0.0;
    for (Long64_t jentry=0; jentry<nentries_use;jentry++)
    {
        Long64_t ientry = LoadTree(jentry);
       
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        // if (check_event_number>0 && event!=check_event_number) continue;
        httEvent->clear();
        evtWriter->setDefault();


        debugWayPoint("## FOUND EVENT ##", {},{(int)event});

        if(jentry%10000==0)
        {   
            perc = jentry > 0 ? (float)jentry/(float)nentries : 0.0;
            std::cout<<"Processing "<<jentry<<"th event  "<< perc << "%" <<std::endl;
        }

        if( !eventInJson() ) continue;
        debugWayPoint("[Loop] passes json");

        unsigned int bestPairIndex = Cut(ientry);
        debugWayPoint("[Loop] best pair index", {}, {(int)bestPairIndex});

        hStats->Fill(0);//Number of events analyzed
        // hStats->Fill(1,httEvent->getMCWeight());//Sum of weights

        if ( failsGlobalSelection() ) continue;
        debugWayPoint("[Loop] passes global selection");

        bestPairIndex_ = bestPairIndex;

        if(bestPairIndex<9999)
        {

            debugWayPoint("[Loop] good pair index found");

            ///Call pairSelection again to set selection bits for the selected pair.
            pairSelection(bestPairIndex);

            // Throw away a lot of events if not producing sync ntuples.
            if(!isSync && !httEvent->checkSelectionBit(SelectionBitsEnum::antiLeptonId) ) continue;

            fillJets(bestPairIndex);
            fillGenLeptons();
            fillPairs(bestPairIndex);
            fillEvent(bestPairIndex);

            HTTPair & bestPair = httPairCollection[0];
            applyMetRecoilCorrections(bestPair); // Adds met to pair 

            if( !httEvent->checkSelectionBit(SelectionBitsEnum::thirdLeptonVeto)
                && !httEvent->checkSelectionBit(SelectionBitsEnum::diLeptonVeto)
            ){
                computeSvFit(bestPair);
            }


            evtWriter->fill(httEvent.get(), &httJetCollection, httLeptonCollection, &bestPair);
            evtWriter->entry=entry++;
            evtWriter->fileEntry=jentry;
            t_TauCheck->Fill();

            hStats->Fill(2);//Number of events saved to ntuple
            hStats->Fill(3,httEvent->getMCWeight());//Sum of weights saved to ntuple
            if(firstWarningOccurence_) firstWarningOccurence_ = false; //stop to warn once the first pair is found and filled
        }
    }

    cout << "Found " << entry << " pairs" << endl;

   //if you change any of the lists, uncomment and execute in this directory (i.e. not running parallel in another) and then run again.
   //  Or, if you run in another directory, copy the *Enum*h over after running the first time, then run again.
   
   // writeJECSourceHeader(jecSources_);
   // writePropertiesHeader(leptonPropertiesList);
   // writeTriggersHeader(triggerBits_);
   // writeFiltersHeader(filterBits_);
   
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::failsGlobalSelection()
{

  //  if ( getMetFilterBits() != passMask_ ) return true;


  return false;
}

/////////////////////////////////////////////////
/*

Building pairs:
Cut
  1 fillLeptons
  2 buildPairs (fill info for all possible pairs, and sort them according to comparePairs)
  3 pairSelection [channel-specific, e.g. in HMuTau...] (filter out pairs not fulfilling kinematic/ID/deltaR requirements)
  4 bestPair (select pair with index 0)

*/
/////////////////////////////////////////////////

/////////////////////////////////////////////////
/////////////////////////////////////////////////
Int_t HTauTauTreeFromNanoBase::Cut(Long64_t entry)
{

    debugWayPoint("[Cut] ------ Begin -------");
    debugWayPoint("[Cut] Fill leptons");
    fillLeptons();

    if( !(httLeptonCollection.size()>1) ) return 9999;

    debugWayPoint("[Cut] Found at least two leptons");

    //build pairs
    if(!buildPairs()) return 9999;

    debugWayPoint("[Cut] Found pairs", {},{(int)httPairs_.size()});

    std::vector<unsigned int> pairIndices;
    for(unsigned int iPair=0;iPair<httPairs_.size();iPair++)
    {
        debugWayPoint("[Cut] Check pair", {},{(int)iPair});
        if(pairSelection(iPair))
        {
            pairIndices.push_back(iPair);
            debugWayPoint("[Cut] Pair passes pairSelection");
        }
    }
    debugWayPoint("[Cut] Pairs passing pairSelection", {}, {(int)pairIndices.size()  });
    
    return bestPair(pairIndices);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauTauTreeFromNanoBase::bestPair(std::vector<unsigned int> &pairIndices)
{
    ///Pair are already sorted during the ntuple creation
    if(!pairIndices.empty()) return pairIndices[0];
    else return 9999;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::pairSelection(unsigned int iPair)
{

    ///Requires channel specific implementation based on the following
    ///Baseline+post sync selection as on
    ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015#Baseline_mu_tau_h_AN1
    ///Indexes for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
    ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

    return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::extraMuonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin)
{
    return thirdLeptonVeto(signalLeg1Index,signalLeg2Index,13,dRmin);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::extraElectronVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin)
{
    return thirdLeptonVeto(signalLeg1Index,signalLeg2Index,11,dRmin);
}

bool HTauTauTreeFromNanoBase::thirdLeptonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, int leptonPdg, double dRmin)
{

    TLorentzVector leg1P4 = httLeptonCollection[signalLeg1Index].getP4();
    TLorentzVector leg2P4 = httLeptonCollection[signalLeg2Index].getP4();

    for(unsigned int iLepton=0;iLepton<httLeptonCollection.size();iLepton++)
    {
        if(iLepton==signalLeg1Index || iLepton==signalLeg2Index) continue;
        TLorentzVector leptonP4 = httLeptonCollection[iLepton].getP4();
        double dr = std::min(leg1P4.DeltaR(leptonP4),leg2P4.DeltaR(leptonP4));
        if(dr<dRmin) continue;
        if(leptonPdg == 13 && std::abs(httLeptonCollection[iLepton].getPDGid())==leptonPdg && httLeptonCollection[iLepton].isExtraLepton() ) return true;
        if(leptonPdg == 11 && std::abs(httLeptonCollection[iLepton].getPDGid())==leptonPdg && httLeptonCollection[iLepton].isExtraLepton() ) return true;
    }
    return false;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillEvent(unsigned int bestPairIndex)
{

    // httEvent->setSelectionBit(SelectionBitsEnum::diMuonVeto,1); //only set explicitly for mutau
    // httEvent->setSelectionBit(SelectionBitsEnum::diElectronVeto,1); //only set explicitly for etau

    httEvent->setRun(run);
    httEvent->setEvent(event);
    httEvent->setLS(luminosityBlock);
    httEvent->setNPV(PV_npvs);
    httEvent->setRho(fixedGridRhoFastjetAll);

    httEvent->setAODPV(TVector3(PV_x,PV_y,PV_z));
    httEvent->setRefittedPV(TVector3(PV_x,PV_y,PV_z));//FIXME
    httEvent->setIsRefit(false);//FIXME

    TVector2 metPF;
    metPF.SetMagPhi(MET_pt, MET_phi);
    httEvent->setMET(metPF);
    httEvent->setMET_uncorr(metPF);

    std::vector<Int_t> aFilters = getFilters(filterBits_);
    httEvent->setFilters(aFilters);

    httEvent->setMETFilterDecision(getMetFilterBits());

    if( isMC )//Assume that all those are filled for MC
    {

        httEvent->setStage1Cat( GenHiggs_stage1PtJet30 );

        httEvent->setMCWeight( sgn(genWeight) );
        httEvent->setXsec( Settings["xsec"].get<float>() );
        httEvent->setGenNEventsWeight( 1.0 / Settings["genNEvents"].get<float>() );

        httEvent->setMCatNLOWeight(LHEWeight_originalXWGTUP);//??
        httEvent->setLHE_Ht(LHE_HT);
        httEvent->setLHEnOutPartons(LHE_Njets);
        //FIXMEhttEvent->setGenPV(TVector3(pvGen_x,pvGen_y,pvGen_z));

        httEvent->setNPU(Pileup_nTrueInt); //??Pileup_nPU or Pileup_nTrueInt
        httEvent->setPUWeight( puweights_histo->GetBinContent( puweights_histo->GetXaxis()->FindBin(Pileup_nTrueInt) ) );

        // Zpt reweighting
        TLorentzVector genBosonP4, genBosonVisP4;
        float zPtReWeight = 1.;
        if( findBosonP4(genBosonP4,genBosonVisP4) )
        {
            httEvent->setGenBosonP4(genBosonP4,genBosonVisP4);

            if( ( httEvent->getSampleType() == HTTEvent::DY || httEvent->getSampleType() == HTTEvent::DYLowM ) )
            {
                zPtReWeight = getZPtReweight(genBosonP4);
            }
        }
        httEvent->setZPtReWeight(zPtReWeight);
        /////////////////////////////////////////////////////////////////////////////////////////////////////

        ///TT reweighting according to
        ///https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#pt_top_Reweighting
        TLorentzVector topP4, antitopP4;
        double topPtReWeight = 1.;
        double topPtReWeight_r1 = 1.;
        if( findTopP4(topP4, antitopP4)  )
        {
            httEvent->setTopP4(topP4, antitopP4);
            if( httEvent->getSampleType() == HTTEvent::TTbar)
            {
                double topPt     = topP4.Perp()      > 400 ? 400 : topP4.Perp() ;
                double antitopPt = antitopP4.Perp()  > 400 ? 400 : antitopP4.Perp();

                double weightTop = exp(0.0615-0.0005*topPt);
                double weightAntitop= exp(0.0615-0.0005*antitopPt);

                double weightTop_r1 = exp(0.156-0.00137*topPt);
                double weightAntitop_r1= exp(0.156-0.00137*antitopPt);

                topPtReWeight = sqrt(weightTop*weightAntitop);
                topPtReWeight_r1 = sqrt(weightTop_r1*weightAntitop_r1);
            }
        }
        httEvent->setTopPtReWeight(topPtReWeight);
        httEvent->setTopPtReWeightR1(topPtReWeight_r1);
        /////////////////////////////////////////////////////////////////////////////////////////////////////
        if(nnlo_ggh_graphs)
        {
            if      (GenHiggs_njets30==0)      httEvent->setNNLO_ggH_weight( NNLOPSratio_pt_mcatnlo_0jet->Eval( GenHiggs_pt > 125.0 ? 125.0 : GenHiggs_pt ) );
            else if (GenHiggs_njets30==1)      httEvent->setNNLO_ggH_weight( NNLOPSratio_pt_mcatnlo_1jet->Eval( GenHiggs_pt > 625.0 ? 625.0 : GenHiggs_pt ) );
            else if (GenHiggs_njets30==2)      httEvent->setNNLO_ggH_weight( NNLOPSratio_pt_mcatnlo_2jet->Eval( GenHiggs_pt > 800.0 ? 800.0 : GenHiggs_pt ) );
            else if (GenHiggs_njets30>=3)      httEvent->setNNLO_ggH_weight( NNLOPSratio_pt_mcatnlo_3jet->Eval( GenHiggs_pt > 925.0 ? 925.0 : GenHiggs_pt ) );
            else                               httEvent->setNNLO_ggH_weight( 1.0 );
        }
    }

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::initJecUnc(std::string correctionFile)
{

    for(unsigned int isrc = 0; isrc < (unsigned int)JecUncertEnum::NONE; isrc++)
    {
        JetCorrectorParameters *p = new JetCorrectorParameters(correctionFile, JecUncertNames[isrc]);
        JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
        jecUncerts.push_back(unc);
        // outputFile<<jecSources_[isrc]<<" = "<<isrc<<", "<<std::endl;
    }
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double HTauTauTreeFromNanoBase::getJecUnc(unsigned int index, unsigned int isrc, bool up)
{
    if(b_nGenPart==nullptr) return 0;//MB: do not check it for data
    double result = 0;
    double jetpt = Jet_pt[index];
    double jeteta =  Jet_eta[index];

    JetCorrectionUncertainty *unc = jecUncerts[isrc];
    unc->setJetPt(jetpt);
    unc->setJetEta(jeteta);
    return unc->getUncertainty(up);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
map<string,double> HTauTauTreeFromNanoBase::getValuesAfterJecSplitting(unsigned int iJet)
{   
    map<string,double> values = { {"", 0.} };
    for(auto uncert : JecAfterSplitting)
    {   
        double uncertSum = 0.;
        for(auto source : uncert.second)
        {
            uncertSum = sqrt( uncertSum*uncertSum +  pow( getJecUnc(iJet, (unsigned int)source, true), 2 )  );
        }
        values[uncert.first] = uncertSum ;
    }
    return values;
}

bool HTauTauTreeFromNanoBase::jetSelection(unsigned int index, unsigned int bestPairIndex)
{

    TLorentzVector aP4;
    aP4.SetPtEtaPhiM(Jet_pt[index],
                     Jet_eta[index],
                     Jet_phi[index],
                     Jet_mass[index]);

    bool passSelection = std::abs(aP4.Eta())<4.7
                         && Jet_jetId[index]>=1;//it means at least loose
   
    if(bestPairIndex<9999)
    {
        TLorentzVector leg1P4 = httPairs_[bestPairIndex].getLeg1().getP4();
        TLorentzVector leg2P4 = httPairs_[bestPairIndex].getLeg2().getP4();

        passSelection &= aP4.DeltaR(leg1P4) > 0.5
                         && aP4.DeltaR(leg2P4) > 0.5;
    }

    if( 2.65 < std::abs( aP4.Eta() ) && std::abs( aP4.Eta() ) <  3.139 && aP4.Pt() < 50) return false; // removal of jets from EE noise
    return passSelection;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillJets(unsigned int bestPairIndex)
{

    debugWayPoint("[fillJets] ------ Begin -------");
    httJetCollection.clear();

    for(unsigned int iJet=0;iJet<nJet;++iJet)
    {

        if(!jetSelection(iJet, bestPairIndex)) continue;
        debugWayPoint("[fillJets] passing jet selection",{(double)Jet_pt[iJet], (double)Jet_eta[iJet]},{},{"pt","eta"});


        HTTJet aJet;
        aJet.SetPtEtaPhiM(Jet_pt[iJet],
                        Jet_eta[iJet],
                        Jet_phi[iJet],
                        Jet_mass[iJet]);


        ///JEC uncertaintes
        aJet.setJecUncertValues( getValuesAfterJecSplitting(iJet) );

        std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iJet, aJet.P4(), "Jet");
        ///Set jet PDG id by hand
        aProperties[(unsigned int)PropertyEnum::pdgId] = 98.0;

        aJet.setProperties(aProperties);
        httJetCollection.addJet(aJet);
        
    }
    //Set Jet collection to unshifted jets
    httJetCollection.fillCurrentCollections();

}
int HTauTauTreeFromNanoBase::muonSelection(HTTParticle aLepton)
{
    int bitmask = 0;
    float muonPt = aLepton.getP4().Pt();
    float muonEta = std::abs( aLepton.getP4().Eta() );
 

    if(std::abs(aLepton.getProperty(PropertyEnum::dz))<0.2
       && std::abs(aLepton.getProperty(PropertyEnum::dxy))<0.045)
    {
        float muonIso = aLepton.getProperty( HTTEvent::usePropertyFor.at("muonIsolation") ) < 0.3;
        float muonID = aLepton.getProperty(  HTTEvent::usePropertyFor.at("muonID")) > 0.5;

        // Passes Baseline cuts 
        if(muonPt > LeptonCuts::Baseline.Muon.pt
           && muonEta < LeptonCuts::Baseline.Muon.eta
           && muonID
        ) bitmask |= LeptonCuts::Baseline.bitmask;

        
        if(muonIso)
        {
            // Passes ilepton cuts
            if(muonPt > LeptonCuts::Di.Muon.pt
               && muonEta < LeptonCuts::Di.Muon.eta
            ) bitmask |= LeptonCuts::Di.bitmask; 

            if(muonID)
            {
                // Passes extra lepton cuts
                if(muonPt > LeptonCuts::Extra.Muon.pt
                   && muonEta < LeptonCuts::Extra.Muon.eta
                ) bitmask |= LeptonCuts::Extra.bitmask;

                // Passes additional lepton cuts
                if(muonPt > LeptonCuts::Additional.Muon.pt
                   && muonEta < LeptonCuts::Additional.Muon.eta
                ) bitmask |= LeptonCuts::Additional.bitmask;
            }
        }

    }
    debugWayPoint("[muonSelection]",
        {(double)muonPt,
         (double)muonEta,
         (double)aLepton.getProperty(PropertyEnum::dz) ,
         (double)aLepton.getProperty(PropertyEnum::dxy),         
         (double)aLepton.getProperty( HTTEvent::usePropertyFor.at("muonIsolation") ),
         (double)aLepton.getProperty( HTTEvent::usePropertyFor.at("muonID")),
         (double)bitmask
        },
        {},
        {"pt","eta","dz","dxy","iso","id","mask"}
    ); 
    return bitmask;
}

int HTauTauTreeFromNanoBase::electronSelection(HTTParticle aLepton)
{
    int bitmask = 0;
    float elePt = aLepton.getP4().Pt();
    float eleEta = std::abs( aLepton.getP4().Eta() );


    if(std::abs(aLepton.getProperty(PropertyEnum::dz))<0.2
       && std::abs(aLepton.getProperty(PropertyEnum::dxy))<0.045)
    {

        bool eleIso =    aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIsolation") ) < 0.3;
        bool eleIDWP80 = aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIDWP80")) > 0.5;
        bool eleIDWP90 = aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIDWP90")) > 0.5;
        bool eleIDCB   = aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIDCutBased") ) > 0;
        bool convVeto  = aLepton.getProperty(PropertyEnum::convVeto)>0.5;
        bool lostHits  = aLepton.getProperty(PropertyEnum::lostHits)<1.5; //0 or 1


        // Passes dilepton cuts
        if(elePt >  LeptonCuts::Di.Electron.pt
            && eleEta < LeptonCuts::Di.Electron.eta
            && eleIso
            && eleIDCB
        ) bitmask |= LeptonCuts::Di.bitmask;

        if(convVeto && lostHits)
        {
            // Passes baseline cuts
            if(elePt >  LeptonCuts::Baseline.Electron.pt
                && eleEta < LeptonCuts::Baseline.Electron.eta
                && eleIDWP90
            ) bitmask |= LeptonCuts::Baseline.bitmask;

            // Passes additional lepton cuts
            if(elePt >  LeptonCuts::Additional.Electron.pt
                && eleEta < LeptonCuts::Additional.Electron.eta
                && eleIDWP90
                && eleIso
            ) bitmask |= LeptonCuts::Additional.bitmask;

            // Passes extra lepton cuts
            if(elePt >  LeptonCuts::Extra.Electron.pt
                && eleEta < LeptonCuts::Extra.Electron.eta
                && eleIDWP90
                && eleIso
            ) bitmask |= LeptonCuts::Extra.bitmask;
        }

    }

    debugWayPoint("[electronSelection]",
        {(double)elePt,
         (double)eleEta,
         (double)aLepton.getProperty(PropertyEnum::dz) ,
         (double)aLepton.getProperty(PropertyEnum::dxy),
         (double)aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIsolation") ),
         (double)aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIDWP80")),
         (double)aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIDWP90")),
         (double)aLepton.getProperty(HTTEvent::usePropertyFor.at("electronIDCutBased")),
         (double)aLepton.getProperty(PropertyEnum::convVeto),
         (double)aLepton.getProperty(PropertyEnum::lostHits),
         (double)bitmask
        },
        {},
        {"pt","eta","dz","dxy","iso","id80","id90","idcb","coVe","loHi","mask"}
    );

    return bitmask;
}

int HTauTauTreeFromNanoBase::tauSelection(HTTParticle aLepton)
{
    int bitmask = 0;

    float tauPt = aLepton.getP4().Pt();
    float tauEta = std::abs( aLepton.getP4().Eta() );

    if(tauPt > LeptonCuts::Baseline.Tau.SemiLep.pt
       && tauEta < LeptonCuts::Baseline.Tau.SemiLep.eta
    ) bitmask |= (LeptonCuts::Baseline.Tau.SemiLep.bitmask + LeptonCuts::Baseline.Tau.Additional);

    if(tauEta < LeptonCuts::Baseline.Tau.FullHad.eta)
    {
        if(tauPt  > LeptonCuts::Baseline.Tau.FullHad.lead_pt)
            bitmask |= (LeptonCuts::Baseline.Tau.FullHad.lead_bitmask + LeptonCuts::Baseline.Tau.Additional);

        if(tauPt  > LeptonCuts::Baseline.Tau.FullHad.sublead_pt)
            bitmask |= (LeptonCuts::Baseline.Tau.FullHad.sub_bitmask + LeptonCuts::Baseline.Tau.Additional);
    }

    debugWayPoint("[tauSelection]",
        {(double)tauPt,
         (double)tauEta,
         (double)aLepton.getProperty(PropertyEnum::dz) ,
         (double)aLepton.getProperty(HTTEvent::usePropertyFor.at("tauIsolation") ),
         (double)aLepton.getProperty(HTTEvent::usePropertyFor.at("tauID") ),
         (double)bitmask
        },
        {},
        {"pt","eta","dz","iso","id","mask"}
    );

    return bitmask;
}

void HTauTauTreeFromNanoBase::fillLeptons()
{

    debugWayPoint("[fillLeptons] ------ Begin -------");
    httLeptonCollection.clear();

    // Only interested in events with hadronic taus
    if(nTau == 0) return;

    //Muons //////////////////////////////////////////////////////////////////////////////////
    for(unsigned int iMu=0; iMu<nMuon; ++iMu)
    {
        float loosestMuonPtCut = min( { LeptonCuts::Baseline.Muon.pt,
                                        LeptonCuts::Additional.Muon.pt,
                                        LeptonCuts::Extra.Muon.pt, 
                                        LeptonCuts::Di.Muon.pt
                                      } );

        if( !(Muon_pt[iMu] > loosestMuonPtCut ) ) continue;
        debugWayPoint("[fillLeptons] Muon passes loosest pt cut",{(double)Muon_pt[iMu], (double)loosestMuonPtCut },{}, {"pt","cut"} );
        HTTParticle aLepton;
        TLorentzVector p4;

        p4.SetPtEtaPhiM(Muon_pt[iMu],
                        Muon_eta[iMu],
                        Muon_phi[iMu],
                        0.10566); //muon mass
        
        TVector3 pca;//FIXME: can partly recover with ip3d and momentum?
        std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iMu, p4, "Muon");
        aLepton.setProperties(aProperties);
        aLepton.setP4(p4);
        aLepton.setChargedP4(p4);//same as p4 for muon
        //aLepton.setNeutralP4(p4Neutral); not defined for muon
        aLepton.setPCA(pca);

        aLepton.setCutBitmask( muonSelection(aLepton) );

        httLeptonCollection.push_back(aLepton);
    }//Muons

    //Electrons //////////////////////////////////////////////////////////////////////////////////
    for(unsigned int iEl=0; iEl<nElectron; ++iEl)
    {
        float e_pt=Electron_pt[iEl];
        float loosestElectronPtCut = min( { LeptonCuts::Baseline.Electron.pt,
                                            LeptonCuts::Additional.Electron.pt,
                                            LeptonCuts::Extra.Electron.pt, 
                                            LeptonCuts::Di.Electron.pt

                                          } );

        if (Electron_eCorr[iEl]>0) e_pt/=Electron_eCorr[iEl];
        if( !(e_pt>loosestElectronPtCut) ) continue;
        debugWayPoint("[fillLeptons] Electron passes loosest pt cut",{(double)e_pt,(double)loosestElectronPtCut },{}, {"pt","cut"} );
        HTTParticle aLepton;
        TLorentzVector p4;

        p4.SetPtEtaPhiM(e_pt,
                        Electron_eta[iEl],
                        Electron_phi[iEl],
                        0.51100e-3); //electron mass

        TVector3 pca;//FIXME: can partly recover with ip3d and momentum?

        std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iEl, p4, "Electron");
        aLepton.setProperties(aProperties);
        aLepton.setP4(p4);
        aLepton.setChargedP4(p4);//same as p4 for electron
        //aLepton.setNeutralP4(p4Neutral); not defined for electron
        aLepton.setPCA(pca);

        aLepton.setCutBitmask( electronSelection(aLepton) );
        httLeptonCollection.push_back(aLepton);
    }//Electrons

    //Taus //////////////////////////////////////////////////////////////////////////////////
    for(unsigned int iTau=0; iTau<nTau; ++iTau)
    {
        if( std::abs(Tau_eta[iTau])>2.3 ) continue;
        if( Tau_idDecayMode[iTau]<0.5 ) continue; //oldDMs
        debugWayPoint("[fillLeptons] Tau passes eta cut and DMId");

        HTTParticle aLepton;

        TLorentzVector newp4;
        newp4.SetPtEtaPhiM(Tau_pt[iTau],
                           Tau_eta[iTau],
                           Tau_phi[iTau],
                           Tau_mass[iTau]);

        std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iTau, newp4, "Tau");
        aLepton.setProperties(aProperties); //Set properties to allow calculation of TES in HTTParticle (mc_match needed)
        aLepton.setP4(newp4);
        aLepton.setCutBitmask( tauSelection(aLepton) );

        if( aLepton.getP4().Pt() < 20 ) continue;
        if( std::abs(aLepton.getProperty(PropertyEnum::dz)) > 0.2 ) continue;
        if( (int)std::abs(aLepton.getProperty(PropertyEnum::charge)) != 1 )continue;
        debugWayPoint("[fillLeptons] Tau passes loosest pt cut after ES",{(double)aLepton.getP4().Pt()},{(int)aLepton.getProperty(PropertyEnum::decayMode)},{"pt","dm"});

        UChar_t bitmask=aLepton.getProperty( HTTEvent::usePropertyFor.at("tauID") ); //byIsolationMVArun2v1DBoldDMwLTraw
        if ( !(bitmask & 0x1 ) ) continue; //require at least very loose tau (in NanoAOD, only OR of loosest WP of all discriminators is stored)

        TLorentzVector chargedP4;//approximate by leadTrack p4
        double leadTkPhi = TVector2::Phi_mpi_pi(Tau_phi[iTau]+Tau_leadTkDeltaPhi[iTau]);

        chargedP4.SetPtEtaPhiM(Tau_pt[iTau]*Tau_leadTkPtOverTauPt[iTau],
                               Tau_eta[iTau]+Tau_leadTkDeltaEta[iTau],
                               leadTkPhi,
                               0.13957); //pi+/- mass

        aLepton.setChargedP4(chargedP4);
        aLepton.setNeutralP4( aLepton.getP4() -chargedP4);

        TVector3 pca;//FIXME: can partly recover with dxy,dz and momentum?
        aLepton.setPCA(pca);

        httLeptonCollection.push_back(aLepton);
    }//Taus

    //Sort leptons
    std::sort(httLeptonCollection.begin(),httLeptonCollection.end(),compareLeptons);

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillGenLeptons()
{

    httGenLeptonCollection.clear();

    if(!fChain->FindBranch("nGenPart")) return;

    for(unsigned int iGenPart=0;iGenPart<nGenPart;++iGenPart)
    {
        if(std::abs(GenPart_pdgId[iGenPart])!=15) continue;
        TLorentzVector p4;
        p4.SetPtEtaPhiM(GenPart_pt[iGenPart],
                        GenPart_eta[iGenPart],
                        GenPart_phi[iGenPart],
                        1.777);//should use pdg mass as masses below 10GeV are zeroed
        //do not consider low momentum candidates??
        if( !(p4.P()>10) ) continue;
        //find direct daughters
        std::vector<unsigned int> daughterIndexes;
        bool isFinalTau=getDirectDaughterIndexes(daughterIndexes,(int)iGenPart);
        if(!isFinalTau) continue;

        HTTParticle aLepton;
        aLepton.setP4(p4);
        aLepton.setChargedP4(getGenComponentP4(daughterIndexes,1));
        aLepton.setNeutralP4(getGenComponentP4(daughterIndexes,0));
        //TVector3 pca(genpart_pca_x->at(iGenPart), genpart_pca_y->at(iGenPart), genpart_pca_z->at(iGenPart));
        //aLepton.setPCA(pca);

        //std::vector<Double_t> aProperties = getProperties(genLeptonPropertiesList, iGenPart);
        //set properties by hand (keep correct order)
        std::vector<Double_t> aProperties;
        aProperties.push_back(GenPart_pdgId[iGenPart]);
        aProperties.push_back(genTauDecayMode(daughterIndexes));    
        aLepton.setProperties(aProperties);

        httGenLeptonCollection.push_back(aLepton);

    }
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector HTauTauTreeFromNanoBase::getGenComponentP4(std::vector<unsigned int> &indexes, unsigned int iAbsCharge)
{

    TLorentzVector aNeutralP4, aChargedP4, aLeptonP4;

    for(unsigned int idx=0;idx<indexes.size();++idx){
        unsigned int iGenPart = indexes[idx];
        unsigned int pdg_id = std::abs(GenPart_pdgId[iGenPart]);
        if(pdg_id == 11 || pdg_id == 13){
            aLeptonP4.SetPtEtaPhiM(GenPart_pt[iGenPart],
                                   GenPart_eta[iGenPart],
                                   GenPart_phi[iGenPart],
                                   (pdg_id==11?0.51100e-3:0.10566));//set mass

        }else if(pdg_id == 211 || pdg_id == 321 ){
            aChargedP4.SetPtEtaPhiM(GenPart_pt[iGenPart],
                                    GenPart_eta[iGenPart],
                                    GenPart_phi[iGenPart],
                                    (pdg_id==211?0.1396:0.4937));//set mass

        }else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 || pdg_id == 311 ){
            aNeutralP4.SetPtEtaPhiM(GenPart_pt[iGenPart],
                                    GenPart_eta[iGenPart],
                                    GenPart_phi[iGenPart],
                                    (pdg_id==111?0.1350:0.4976));//set mass
        }
    }

    TLorentzVector aP4;
    if(iAbsCharge==0) aP4 = aNeutralP4;
    else if(aChargedP4.E()>0) aP4 = aChargedP4;
    else if(aLeptonP4.E()>0) aP4 = aLeptonP4;

    return aP4;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getGenMatch(unsigned int index, std::string colType) //overloaded
{ 
    if(!isMC) return -999;
    if(b_nGenPart==nullptr) return -999;
    if(nGenPart==0) return -999;

    TLorentzVector p4;
    if(colType=="Muon")
    {
        if(index>=nMuon) return -999;
        p4.SetPtEtaPhiM(Muon_pt[index],
                        Muon_eta[index],
                        Muon_phi[index],
                        0.10566); //muon mass

    }else if(colType=="Electron"){
        if(index>=nElectron) return -999;
        p4.SetPtEtaPhiM(Electron_pt[index],
                        Electron_eta[index],
                        Electron_phi[index],
                        0.51100e-3); //electron mass

    }else if(colType=="Tau"){
        if(index>=nTau) return -999;
        p4.SetPtEtaPhiM(Tau_pt[index],
                        Tau_eta[index],
                        Tau_phi[index],
                        Tau_mass[index]);
    }
    else
        return -999;

    return getGenMatch(p4);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getGenMatch(TLorentzVector selObj)
{

    float dRTmp = 1.;
    float matchings[5] = {15.,15.,15.,15.,15.};
    debugWayPoint("[getGenMatch] ------ Begin -------");
    debugWayPoint("[getGenMatch] particle with",
                  {(double)selObj.Pt(),(double)selObj.Eta(),(double)selObj.Phi()},
                  {},
                  {"pt","eta","phi"});

    for(unsigned int iGen=0;iGen<nGenPart;++iGen)
    {

        if( GenPart_pt[iGen] > 8)
        {
            dRTmp = evtWriter->calcDR( GenPart_eta[iGen],GenPart_phi[iGen],selObj.Eta(),selObj.Phi() );

            int statusFlags=GenPart_statusFlags[iGen];
            bool GenPart_isPrompt=(statusFlags & (1<<0) ) == (1<<0);
            bool GenPart_isDirectPromptTauDecayProduct=(statusFlags & (1<<5) ) == (1<<5);

            //electron
            if(std::abs(GenPart_pdgId[iGen]) == 11)
            {
                if( dRTmp < matchings[0] && GenPart_isPrompt){
                    matchings[0] = dRTmp;
                    debugWayPoint("[getGenMatch] is prompt electron with",{(double)GenPart_pt[iGen], (double)dRTmp},{},{"pt","dR"});
                }
                if( dRTmp < matchings[2] && GenPart_isDirectPromptTauDecayProduct){
                    matchings[2] = dRTmp;
                    debugWayPoint("[getGenMatch] is direct prompt electron with",{(double)GenPart_pt[iGen], (double)dRTmp},{},{"pt","dR"});
                }
            }

            //muon
            if(fabs(GenPart_pdgId[iGen]) == 13)
            {
                if( dRTmp < matchings[1] && GenPart_isPrompt){
                    matchings[1] = dRTmp;
                    debugWayPoint("[getGenMatch] is prompt muon with",{(double)GenPart_pt[iGen], (double)dRTmp},{},{"pt","dR"});
                }
                if( dRTmp < matchings[3] && GenPart_isDirectPromptTauDecayProduct){
                    matchings[3] = dRTmp;
                    debugWayPoint("[getGenMatch] is direct prompt muon with",{(double)GenPart_pt[iGen], (double)dRTmp},{},{"pt","dR"});
                }
            }

            //tauhad
            if( fabs(GenPart_pdgId[iGen]) == 15 && GenPart_isPrompt)
            {
                TLorentzVector tmpLVec;
                TLorentzVector remParticles;
                remParticles.SetPtEtaPhiM(0.,0.,0.,0.);
                int nr_neutrinos = 0;
                int nr_gammas = 0;
                bool vetoLep = false;

                for(unsigned int iDau=0;iDau<nGenPart;++iDau)
                {
                    if( ( fabs(GenPart_pdgId[iDau]) == 11 || fabs(GenPart_pdgId[iDau]) == 13 )
                        && (int)GenPart_genPartIdxMother[iDau] == (int)iGen
                        ) vetoLep = true;

                    if( fabs(GenPart_pdgId[iDau]) == 16 && (int)GenPart_genPartIdxMother[iDau] ==  (int)iGen)
                    {

                        tmpLVec.SetPtEtaPhiM(GenPart_pt[iDau],
                                             GenPart_eta[iDau],
                                             GenPart_phi[iDau],
                                             GenPart_mass[iDau]
                                             );
                        remParticles += tmpLVec;
                        if(fabs(GenPart_pdgId[iDau]) == 16) nr_neutrinos++;
                    }

                }
              
                if(vetoLep==false && nr_neutrinos == 1 )
                {
                    TLorentzVector tau;
                    tau.SetPtEtaPhiM(GenPart_pt[iGen],
                                     GenPart_eta[iGen],
                                     GenPart_phi[iGen],
                                     GenPart_mass[iGen]
                                     );

                    if((tau-remParticles).Pt() > 15)
                    {
                        dRTmp = evtWriter->calcDR( (tau-remParticles).Eta(), (tau-remParticles).Phi(), selObj.Eta(), selObj.Phi() );
                        if( dRTmp < matchings[4] ){
                            matchings[4] = dRTmp;
                            debugWayPoint("[getGenMatch] genuine tau after radiation",{(double)dRTmp},{},{"dR"});
                        }
                    }
                }
            }//end tauhad
        }
    }

    int whichObj = 1;

    float smallestObj = matchings[0];
    for(int i=1; i<5; i++)
    {
        if(matchings[i] < smallestObj){
            smallestObj = matchings[i];
            whichObj = i+1;
        }
    }
    
    if(whichObj < 6 && smallestObj < 0.2)
    {
        debugWayPoint("[getGenMatch] tightest match",{(double)smallestObj},{(int)whichObj},{"dR","=> gen_match"});
        return whichObj;
    }
    debugWayPoint("[getGenMatch] no match => gen_match = 6");
    return 6;

}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillPairs(unsigned int bestPairIndex)
{
    httPairCollection.clear();
    if(bestPairIndex<9999 &&  bestPairIndex<httPairs_.size()){
        httPairCollection.push_back(httPairs_[bestPairIndex]);
    }
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::buildPairs()
{
    httPairs_.clear();
    for(unsigned int iL1=0; iL1<httLeptonCollection.size()-1; ++iL1)
    {
        for(unsigned int iL2=iL1+1; iL2<httLeptonCollection.size(); ++iL2)
        {

            debugWayPoint("[buildPairs] Try pair",{},{(int)iL1,(int)iL2});
            if( !(httLeptonCollection[iL1].getP4().DeltaR(httLeptonCollection[iL2].getP4())>0.3) ) continue;
            debugWayPoint("[buildPairs] No overlap");


            HTTPair aHTTpair;

            aHTTpair.setLeg1(httLeptonCollection.at(iL1),iL1);
            aHTTpair.setLeg2(httLeptonCollection.at(iL2),iL2);
            
            httPairs_.push_back(aHTTpair); 
        }
    }
    std::sort(httPairs_.begin(),httPairs_.end(),comparePairs);
    return !httPairs_.empty();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
template<class T> T HTauTauTreeFromNanoBase::getBranchValue(const char *branchAddress, unsigned int index)
{
    std::vector<T> *aVector = *(std::vector<T> **)(branchAddress);
    if(aVector->size()<=index){
        return 0;
    }
    return aVector->at(index);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
Int_t  HTauTauTreeFromNanoBase::getFilter(std::string name)
{
    TBranch *branch = fChain->GetBranch(name.c_str());
    if(!branch)
    {
        if(firstWarningOccurence_)
          std::cout<<"Branch: "<<name<<" not found in the TTree."<<std::endl;
        return -999;
    }else{
        TLeaf *leaf = branch->FindLeaf(name.c_str());
        bool decision = leaf!=nullptr ? leaf->GetValue() : false;
        return decision;
    }
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
std::vector<Double_t>  HTauTauTreeFromNanoBase::getProperties(const std::vector<std::string> & propertiesList,
                                                              unsigned int index,
                                                              std::string colType)
{

    TLorentzVector obj;
    obj.SetPtEtaPhiM(0,0,0,999999);
    return getProperties(propertiesList,index,obj,colType);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
std::vector<Double_t>  HTauTauTreeFromNanoBase::getProperties(const std::vector<std::string> & propertiesList,
                                                              unsigned int index,
                                                              TLorentzVector obj,
                                                              std::string colType)
{

    std::vector<Double_t> aProperties;
    for(auto propertyName:propertiesList){
        aProperties.push_back(getProperty(propertyName,index,obj,colType));
    }
    return aProperties;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
Double_t  HTauTauTreeFromNanoBase::getProperty(std::string name, unsigned int index, std::string colType)
{
    TLorentzVector obj;
    obj.SetPtEtaPhiM(0,0,0,999999);
    return getProperty(name,index,obj,colType);
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
Double_t  HTauTauTreeFromNanoBase::getProperty(std::string name, unsigned int index, TLorentzVector obj, std::string colType)
{

    if(name=="mc_match") return getGenMatch(index,colType);
    if(name=="isGoodTriggerType") return getTriggerMatching(index,obj,false,colType);
    if(name=="FilterFired") return getTriggerMatching(index,obj,true,colType);//some overhead due to calling it again with a different option, but kept for backward compatibility (and debug)

    if(colType=="Electron")
    {
      if(name.find("Muon_")!=std::string::npos
           || name.find("Tau_")!=std::string::npos
           || name.find("Jet_")!=std::string::npos
           || false) return 0;

    }else if(colType=="Muon"){
        if(name.find("Electron_")!=std::string::npos
           || name.find("Tau_")!=std::string::npos
           || name.find("Jet_")!=std::string::npos
           || false) return 0;

    }else if(colType=="Tau"){
        if(name.find("Electron_")!=std::string::npos
           || name.find("Muon_")!=std::string::npos
           || name.find("Jet_")!=std::string::npos
           || name.find("pfRelIso03_all")!=std::string::npos
           || name.find("sip3d")!=std::string::npos
           || false) return 0;

        if(name.find("pdgId")!=std::string::npos)
        {
            TBranch *branch = fChain->GetBranch("Tau_charge");
            if(!branch)
            {
                if(firstWarningOccurence_) std::cout<<"Branch: Tau_charge not found in the TTree, return pdgId=-15"<<std::endl;
                return -15;
            }
            TLeaf *leaf = branch->FindLeaf("Tau_charge");
            int charge = leaf!=nullptr ? leaf->GetValue(index) : 0;
            return charge>0 ? -15 : 15; 
        }      
    }
    else if(colType=="Jet"){
        if(name.find("Jet_")==std::string::npos) return 0;
    }

    if(colType!="" && name.find(colType+"_")==std::string::npos){
        name = colType+"_"+name;
    }

    TBranch *branch = fChain->GetBranch(name.c_str());
    if(!branch)
    {
        if(firstWarningOccurence_) std::cout<<"Branch: "<<name<<" not found in the TTree."<<std::endl;
        return 0;
    }

    const char *branchAddress = branch->GetAddress();
    std::string branchClass(branch->GetClassName());
    if(branchClass=="vector<double>") return getBranchValue<double>(branchAddress, index);
    if(branchClass=="vector<int>") return getBranchValue<int>(branchAddress, index);
    if(branchClass=="vector<bool>") return getBranchValue<bool>(branchAddress, index);
    if(branchClass=="vector<Long64_t>") return getBranchValue<Long64_t>(branchAddress, index);
    //otherwise assume that index value makes sense
    if(branchClass==""){//check leafs
        TLeaf *leaf = branch->FindLeaf(name.c_str());
        return leaf!=nullptr ? leaf->GetValue(index) : 0;
    }

    return 0;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// void  HTauTauTreeFromNanoBase::writeJECSourceHeader(const std::vector<string> &jecSources)
// {
//     ofstream outputFile("JecUncEnum.h");
//     outputFile<<"enum class JecUncEnum { ";

//     for(unsigned int isrc = 0; isrc < jecSources_.size(); isrc++)
//     {
//         outputFile<<jecSources_[isrc]<<" = "<<isrc<<", "<<std::endl;
//     }
//     outputFile<<"NONE"<<" = "<<jecSources_.size()<<std::endl;
//     outputFile<<"};"<<std::endl;
//     outputFile.close();
// }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void  HTauTauTreeFromNanoBase::writeTriggersHeader(const std::vector<TriggerData> &triggerBits)
{
    ofstream outputFile("TriggerEnum.h");
    outputFile<<"enum class TriggerEnum { ";

    for(unsigned int iBit=0;iBit<triggerBits.size();++iBit)
    {
        std::string name = triggerBits[iBit].path_name;
        outputFile<<name<<" = "<<iBit<<", "<<std::endl;
    }
    outputFile<<"NONE"<<" = "<<triggerBits.size()<<std::endl;
    outputFile<<"};"<<std::endl;
    outputFile.close();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void  HTauTauTreeFromNanoBase::writeFiltersHeader(const std::vector<std::string> &filterBits)
{
    ofstream outputFile("FilterEnum.h");
    outputFile<<"enum class FilterEnum { ";

    for(unsigned int iItem=0;iItem<filterBits.size();++iItem){
        std::string name = filterBits[iItem];
        outputFile<<name<<" = "<<iItem<<", "<<std::endl;
    }
    outputFile<<"NONE"<<" = "<<filterBits.size()<<std::endl;
    outputFile<<"};"<<std::endl;
    outputFile.close();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void  HTauTauTreeFromNanoBase::writePropertiesHeader(const std::vector<std::string> & propertiesList)
{

    ofstream outputFile("PropertyEnum.h");
    outputFile<<"enum class PropertyEnum { ";

    for(unsigned int iItem=0;iItem<propertiesList.size();++iItem)
    {
        std::string name = propertiesList[iItem];
        std::string pattern = "daughters_";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        pattern = "Daughters";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        pattern = "jets_";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        pattern = "Jet_";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        pattern = "Tau_";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        pattern = "Muon_";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        pattern = "Electron_";
        if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
        outputFile<<name<<" = "<<iItem<<", "<<std::endl;
    }
    outputFile<<"NONE"<<" = "<<propertiesList.size()<<std::endl;
    outputFile<<"};"<<std::endl;
    outputFile.close();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////

std::vector<Int_t>  HTauTauTreeFromNanoBase::getFilters(const std::vector<std::string> & filtersList)
{
    std::vector<Int_t> aFilters;
    for(auto propertyName:filtersList){
        aFilters.push_back(getFilter(propertyName));
    }
    return aFilters;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::isGoodToMatch(unsigned int ind)
{  
    //MB: This method not mandatory as selection of good to match is applied upper in the flow. Anyway, an trivial implementation kept.
    return true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getTriggerMatching(unsigned int index, TLorentzVector p4_1, bool checkBit, std::string colType)
{
    //  TLorentzVector p4_1;
    unsigned int particleId=0;
    double dRmax=0.5;//was 0.25

    if(checkBit) debugWayPoint("[getTriggerMatching] ------ Begin -------");

    if(colType=="Muon")            particleId=13;
    else if(colType=="Electron")   particleId=11;
    else if(colType=="Tau")        particleId=15;
    else return 0;

    int firedBits = 0;
    for(unsigned int iTrg=0; iTrg<triggerBits_.size(); ++iTrg)
    {
        bool decision = false;
        if(checkBit) debugWayPoint("[getTriggerMatching] " + triggerBits_[iTrg].path_name);
        if (p4_1.Pt()<triggerBits_[iTrg].leg1OfflinePt) continue;

        //check if trigger is fired
        TBranch *branch = fChain->GetBranch(triggerBits_[iTrg].path_name.c_str());
        if(branch!=nullptr)
        {
            TLeaf *leaf = branch->FindLeaf(triggerBits_[iTrg].path_name.c_str());
            decision = leaf!=nullptr ? leaf->GetValue() : false;
        }
        if(!decision) continue; // do not check rest if trigger is not fired
        if(checkBit) debugWayPoint("[getTriggerMatching] trigger fired");

        //////////////////////// check legs //////////////////////////////////
        //////  first leg   //////// 
        decision = false;
        if(particleId==triggerBits_[iTrg].leg1Id)
        {
            for(unsigned int iObj=0; iObj<nTrigObj; ++iObj)
            {
                if(TrigObj_id[iObj]!=(int)particleId) continue;
                TLorentzVector p4_trg;
                p4_trg.SetPtEtaPhiM(TrigObj_pt[iObj],
                                    TrigObj_eta[iObj],
                                    TrigObj_phi[iObj],
                                    0.);

                if( !(p4_1.DeltaR(p4_trg)<dRmax) ) continue;
                if(checkBit) debugWayPoint("[getTriggerMatching] matched with trigger object",{(double)TrigObj_pt[iObj],
                                                                                               (double)TrigObj_eta[iObj],
                                                                                               (double)TrigObj_phi[iObj] }, {},
                                                                                               {"pt","eta","phi"});

                if( triggerBits_[iTrg].leg1Pt>0   && !( TrigObj_pt[iObj]            > triggerBits_[iTrg].leg1Pt) ) continue;
                if( triggerBits_[iTrg].leg1Eta>0  && !( std::abs(TrigObj_eta[iObj]) < triggerBits_[iTrg].leg1Eta) ) continue;
                if( triggerBits_[iTrg].leg1L1Pt>0 && !( TrigObj_l1pt[iObj]          > triggerBits_[iTrg].leg1L1Pt) ) continue;
                if(checkBit) debugWayPoint("[getTriggerMatching] object passes l1pt cut",{(double)triggerBits_[iTrg].leg1L1Pt, (double)TrigObj_l1pt[iObj] },{},{"cut","l1pt"});

                if( checkBit && !( ((int)TrigObj_filterBits[iObj] & triggerBits_[iTrg].leg1BitMask)==triggerBits_[iTrg].leg1BitMask) ) continue;
                if(checkBit) debugWayPoint("[getTriggerMatching] passes filter");
                decision = true;
                // break;
            }
        }
        if(decision){
            firedBits |= (1<<iTrg);
            continue;
        }
        //////////////////////////////////////////////////////////////////////
        //////  second leg   //////// 
        decision = false;
        if(particleId==triggerBits_[iTrg].leg2Id)
        {
            for(unsigned int iObj=0; iObj<nTrigObj; ++iObj)
            {
                if(TrigObj_id[iObj]!=(int)particleId) continue;
                TLorentzVector p4_trg;
                p4_trg.SetPtEtaPhiM(TrigObj_pt[iObj],
                                    TrigObj_eta[iObj],
                                    TrigObj_phi[iObj],
                                    0.);

                if( !(p4_1.DeltaR(p4_trg)<dRmax) ) continue;
                if( triggerBits_[iTrg].leg2Pt>0   && !( TrigObj_pt[iObj]            > triggerBits_[iTrg].leg2Pt) ) continue;
                if( triggerBits_[iTrg].leg2Eta>0  && !( std::abs(TrigObj_eta[iObj]) < triggerBits_[iTrg].leg2Eta) ) continue;
                if( triggerBits_[iTrg].leg2L1Pt>0 && !( TrigObj_l1pt[iObj]          > triggerBits_[iTrg].leg2L1Pt) ) continue;

                if( checkBit && !( ((int)TrigObj_filterBits[iObj] & triggerBits_[iTrg].leg2BitMask)==triggerBits_[iTrg].leg2BitMask) ) continue;
                decision = true;
                break;
            }
        }
        if(decision){
            firedBits |= (1<<iTrg);
            continue;
        }
        ////////////////////////////////////////////////////////////////////// 
    }
    return firedBits;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getMetFilterBits()
{
    int firedBits=0;
    for(unsigned int iFlt=0; iFlt<filterBits_.size(); ++iFlt)
    {
        bool decision = false;
        //check if trigger is fired
        TBranch *branch = fChain->GetBranch(filterBits_[iFlt].c_str());
        if(branch!=nullptr)
        {
            TLeaf *leaf = branch->FindLeaf(filterBits_[iFlt].c_str());
            decision = leaf!=nullptr ? leaf->GetValue() : false;
        }
        if(decision) firedBits |= (1<<iFlt);
    }
    return firedBits;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::findBosonP4(TLorentzVector &bosonP4, TLorentzVector &visBosonP4)
{
    bosonP4.SetXYZM(0,0,0,0);
    visBosonP4.SetXYZM(0,0,0,0);

    if(b_nGenPart==nullptr) return false;

    for(unsigned int iGen = 0; iGen < nGenPart; ++iGen)
    {
        unsigned int absPdgId = std::abs(GenPart_pdgId[iGen]);
        if(absPdgId != 11 && absPdgId != 13 && absPdgId != 15 //charged leptons
           && absPdgId != 12 && absPdgId != 14 && absPdgId != 16 ) //neutrinos
          continue;

        bool fromHardProcessFinalState = ( (GenPart_statusFlags[iGen] & (1<<8)) == (1<<8) );
        if ( !(fromHardProcessFinalState ) ) continue;

        //mass is stored only m>10GeV and photons m>1GeV
        double mass = GenPart_mass[iGen]>0 ? GenPart_mass[iGen] :
          absPdgId==15 ? 1.7776 :
          absPdgId==13 ? 0.10566 :
          absPdgId==11 ? 0.51100e-3 : 0.;

        TLorentzVector p4;
        p4.SetPtEtaPhiM(GenPart_pt[iGen],
                        GenPart_eta[iGen],
                        GenPart_phi[iGen],
                        mass);

        std::vector<unsigned int> daughterIndexes;
        bool isFinal=getDirectDaughterIndexes(daughterIndexes,(int)iGen,false);//store neutrinos for further use
        if(!isFinal) continue;

        bool isElectron = (absPdgId == 11);
        bool isMuon = (absPdgId == 13);
        bool isNeutrino = (absPdgId == 12 || absPdgId == 14 || absPdgId == 16);
        bool isTau = (absPdgId == 15);
        /*
        bool isDirectHardProcessTauDecayProduct = true;//FIXME(genFlags & (1<<10)) == (1<<10);
        if ( (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct){
          genBosonP4 += p4;
        }
        ///This GenParticle list is missing pions, so we have
        ///to add hadronic tau, subtract neutral component
        if(absPdgId == 66615) genBosonP4 += p4;
        if(absPdgId == 77715) genBosonP4 -= p4;
        ///
        */
        //std::cout<<"\tpdgId="<<absPdgId<<", fromHardProcess:"<<fromHardProcessFinalState<<std::endl;
       //is the following correct? When undecayed particles are used (incl. taus) there should not be double counting
        if ( fromHardProcessFinalState && (isMuon || isElectron || isNeutrino || isTau) )
        {
            bosonP4 += p4;
            if(!isNeutrino) visBosonP4 += p4;

            if(isTau)
            {
              //subtract p4 of tau neutrinos from visBosonP4
                for(unsigned int iDau=0; iDau<daughterIndexes.size(); ++iDau)
                {
                    unsigned int absPdgIdDau = std::abs(GenPart_pdgId[daughterIndexes[iDau]]);
                    if(absPdgIdDau != 12 && absPdgIdDau != 14 && absPdgIdDau != 16 ) continue; //neutrinos
                      
                    TLorentzVector p4Dau;
                    p4Dau.SetPtEtaPhiM(GenPart_pt[daughterIndexes[iDau]],
                                       GenPart_eta[daughterIndexes[iDau]],
                                       GenPart_phi[daughterIndexes[iDau]],
                                       0.);
                    visBosonP4 -= p4Dau;
                }
            }
        }
    }
    //std::cout<<"boson mass="<<bosonP4.M()<<std::endl;

    return ( bosonP4.M()>1E-3 || bosonP4.P()>1E-3 );
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::findTopP4(TLorentzVector &topP4, TLorentzVector &antiTopP4)
{
    topP4.SetXYZM(0,0,0,0);
    antiTopP4.SetXYZM(0,0,0,0);

    if(b_nGenPart==nullptr) return false;

    for(unsigned int iGen = 0; iGen < nGenPart; ++iGen)
    {
        if(std::abs(GenPart_pdgId[iGen])!=6) continue;

        //mass is stored only m>10GeV and photons m>1GeV (fine for top)
        TLorentzVector p4;
        p4.SetPtEtaPhiM(GenPart_pt[iGen],
                        GenPart_eta[iGen],
                        GenPart_phi[iGen],
                        GenPart_mass[iGen]);

        std::vector<unsigned int> daughterIndexes;
        bool isFinal=getDirectDaughterIndexes(daughterIndexes,(int)iGen);
        if(!isFinal) continue;

        if(GenPart_pdgId[iGen]==6) topP4 = p4;
        if(GenPart_pdgId[iGen]==-6) antiTopP4 = p4;
    }

    return (topP4.M()>1E-3 && antiTopP4.M()>1E-3);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double HTauTauTreeFromNanoBase::getZPtReweight(const TLorentzVector &genBosonP4, bool doSUSY)
{

    double weight = 1.0;

    //Z pt reweighting
    
    if(genBosonP4.M()>1E-3)
    {
        TH2D *hWeight = zptmass_histo;
        double mass = genBosonP4.M();
        double pt = genBosonP4.Perp();    
        int massBin = hWeight->GetXaxis()->FindBin(mass);
        int ptBin = hWeight->GetYaxis()->FindBin(pt);
        weight = hWeight->GetBinContent(massBin,ptBin);
    }
    return weight;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::computeSvFit(HTTPair &aPair)
{
    if(svFitAlgo_==nullptr) return;

    //Legs
    HTTParticle leg1 = aPair.getLeg1();
    double mass1;
    int decay1 = -1;
    classic_svFit::MeasuredTauLepton::kDecayType type1;
    if(std::abs(leg1.getPDGid())==11)
    {
        mass1 = 0.51100e-3; //electron mass
        type1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;

    }else if(std::abs(leg1.getPDGid())==13){
        mass1 = 0.10566; //muon mass
        type1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;

    }else{//tau->hadrs.
        decay1 = leg1.getProperty(PropertyEnum::decayMode);
        mass1 = leg1.getP4().M();
        if(decay1==0) mass1 = 0.13957; //pi+/- mass
        type1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    }

    HTTParticle leg2 = aPair.getLeg2();
    double mass2;
    int decay2 = -1;
    classic_svFit::MeasuredTauLepton::kDecayType type2;
    if(std::abs(leg2.getPDGid())==11){
        mass2 = 0.51100e-3; //electron mass
        type2 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;

    }else if(std::abs(leg2.getPDGid())==13){
        mass2 = 0.10566; //muon mass
        type2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;

    }else{//tau->hadrs.
        decay2 = leg2.getProperty(PropertyEnum::decayMode);
        mass2 = leg2.getP4().M();
        if(decay2==0) mass2 = 0.13957; //pi+/- mass
        type2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
    }
    //Leptons for SvFit
    std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type1, leg1.getP4().Pt(), leg1.getP4().Eta(), leg1.getP4().Phi(), mass1, decay1) );
    measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type2, leg2.getP4().Pt(), leg2.getP4().Eta(), leg2.getP4().Phi(), mass2, decay2) );
    //MET

    TMatrixD covMET(2, 2);
    covMET[0][0] = aPair.getMETMatrix().at(0);
    covMET[0][1] = aPair.getMETMatrix().at(1);
    covMET[1][0] = aPair.getMETMatrix().at(2);
    covMET[1][1] = aPair.getMETMatrix().at(3);
    if(covMET[0][0]==0 && covMET[1][0]==0 && covMET[0][1]==0 && covMET[1][1]==0) return; //singular covariance matrix    

    TLorentzVector p4SVFit;
    vector<string> shifts;

    if(applyRecoil)
    {
        for(auto shift : metShifts)                              shifts.push_back( shift.first );
    }else
    {
        for(auto shift : httJetCollection.getNeededJECShifts() ) shifts.push_back( shift.first );
    }

    for(auto shift : shifts )
    {
        p4SVFit.SetPtEtaPhiM(-10,-.10,-10,-10);

        aPair.setCurrentMETShift(shift);

        //Only calculate svfit for shapes that are in SR
        if( ( strcmp(shift.c_str(),"") != 0 && aPair.isInLooseSR() )
            || ( strcmp(shift.c_str(),"") == 0  )
        ){
            p4SVFit = runSVFitAlgo(measuredTauLeptons, aPair.getMET(), covMET);
        }

        aPair.setP4(p4SVFit,shift);
    }
    aPair.setCurrentMETShift("");
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector HTauTauTreeFromNanoBase::runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons, const TVector2 &aMET, const TMatrixD &covMET)
{

    TLorentzVector p4SVFit;
    if(measuredTauLeptons.size()!=2 || svFitAlgo_==nullptr) return p4SVFit;

    //set logM regularization term which is final state dependent
    double kappa = 4;
    if(measuredTauLeptons[0].type()==classic_svFit::MeasuredTauLepton::kTauToElecDecay 
       || measuredTauLeptons[0].type()==classic_svFit::MeasuredTauLepton::kTauToMuDecay)
    { //1st tau is lepton
        if(measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToElecDecay 
           || measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToMuDecay)
        {
            kappa = 3; //ll decay
        } else{
            kappa = 4; //lt decay
        }

    } else
    { //1st tau is hadron
      if(measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToElecDecay
         || measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToMuDecay)
      {
          kappa = 4; //ltau decay
      } else{
          kappa = 5; //tt decay
      }
    }    
    svFitAlgo_->addLogM_fixed(true, kappa);
    svFitAlgo_->integrate(measuredTauLeptons, aMET.X(), aMET.Y(), covMET);

    if(svFitAlgo_->isValidSolution() )//Get solution
    {

        classic_svFit::DiTauSystemHistogramAdapter* aHistogramAdapter = static_cast< classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_->getHistogramAdapter());
        p4SVFit.SetPtEtaPhiM(aHistogramAdapter->getPt(),
                             aHistogramAdapter->getEta(),
                             aHistogramAdapter->getPhi(),
                             aHistogramAdapter->getMass());

    }else{
        p4SVFit.SetPtEtaPhiM(0,0,0,0);
    }

    return p4SVFit;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::applyMetRecoilCorrections(HTTPair &aPair)
{

    // Shift met by jec if there is not best pair or recoilCorrector is not initialized
    TVector2 met; met.SetMagPhi(MET_pt, MET_phi);
    aPair.setMETMatrix(MET_covXX, MET_covXY, MET_covXY, MET_covYY);

    if( recoilCorrector_==nullptr
        || httPairCollection.empty()
        || !applyRecoil
    )
    {
        for(auto shift : httJetCollection.getNeededJECShifts() )
        {
            aPair.setMET( met - httJetCollection.getTotalJetShift(shift.second.first, shift.second.second), shift.first );
        }
        aPair.setCurrentMETShift(""); //Explicitly set current met to unshifted
        return;
    }

    debugWayPoint("[applyMetRecoilCorrections] original met   ",{(double)met.Px(),
                                                              (double)met.Py(),
                                                              (double)met.Mod()},{},
                                                              {"met_px","met_py","met_pt"} ) ;
    aPair.setMET(met,""); // Set Met to propagate tes

    TLorentzVector genBosonP4 = httEvent->getGenBosonP4();
    TLorentzVector genBosonVisP4 = httEvent->getGenBosonP4(true);

    float corrMEtPx, corrMEtPy;
    float MEtPx = aPair.getMET().Px();
    float MEtPy = aPair.getMET().Py();
    float gen_ll_px = genBosonP4.Px();
    float gen_ll_py = genBosonP4.Py();
    float gen_ll_vis_px = genBosonVisP4.Px();
    float gen_ll_vis_py = genBosonVisP4.Py();

    int nJets = httJetCollection.getNJets(20);
    if(httEvent->getDecayModeBoson()>=10) nJets++; //W, add jet for fake tau
      

    debugWayPoint("[applyMetRecoilCorrections] met after applying tes   ",{(double)MEtPx,
                                                                           (double)MEtPy,
                                                                           (double)TVector2(MEtPx,MEtPy).Mod()},{},
                                                                           {"met_px","met_py","met_pt"} );

    debugWayPoint("[applyMetRecoilCorrections] gen boson   ",{(double)gen_ll_px,
                                                           (double)gen_ll_py,
                                                           (double)gen_ll_vis_px,
                                                           (double)gen_ll_vis_py},{},
                                                           {"gen_px","gen_py","genvis_px","genvis_py"} ) ;

    recoilCorrector_->CorrectByMeanResolution(
    //recoilCorrector_->Correct( //Quantile correction works better for MVA MET
      MEtPx,
      MEtPy,
      gen_ll_px,
      gen_ll_py,
      gen_ll_vis_px,
      gen_ll_vis_py,
      nJets,
      corrMEtPx,
      corrMEtPy
    );

    debugWayPoint("[applyMetRecoilCorrections] met after applying recoil corrections   ",{(double)corrMEtPx,
                                                                                          (double)corrMEtPy,
                                                                                          (double)TVector2(corrMEtPx,corrMEtPy).Mod()},{},
                                                                                          {"met_px","met_py","met_pt"} ) ;


    float met_scale_x, met_scale_y;
    for(auto shift : metShifts )
    {
        met_scale_x = corrMEtPx;
        met_scale_y = corrMEtPy;
        if(strcmp(shift.first.c_str(),"") != 0)
        {
            metSys_->ApplyMEtSys(
                corrMEtPx, corrMEtPy,        
                gen_ll_px,gen_ll_py,         
                gen_ll_vis_px,gen_ll_vis_py, 
                nJets,                       
                MEtSys::ProcessType::BOSON,  
                shift.second.first,   
                shift.second.second,        
                met_scale_x,met_scale_y      
            );
        }     
        aPair.setMET( TVector2(met_scale_x,met_scale_y), shift.first, false); // Set met without applying tes a second time
    }
    aPair.setCurrentMETShift("");


    return ;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::compareLeptons(const HTTParticle& i, const HTTParticle& j)
{
    unsigned int i_type=2, j_type=2;
    if(std::abs(i.getPDGid())==11)      i_type=1;//ele
    else if(std::abs(i.getPDGid())==13) i_type=0;//mu
      
    if(std::abs(j.getPDGid())==11)      j_type=1;//ele
    else if(std::abs(j.getPDGid())==13) j_type=0;//mu
      
    if(i_type > j_type) return false;
    if(i_type == j_type && i.getP4().Pt() < j.getP4().Pt() ) return false;

    return true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::comparePairs(const HTTPair& i, const HTTPair& j)
{
    double i_iso=999, j_iso=999;
    unsigned int i_type=2, j_type=2;

    ////////////////////////////////////// Leg 1 ////////////////////////////////////////////////
    //step 0.5, leg 1 type: 2: tau, 1: e, 0: mu

    if( std::abs(i.getLeg1().getPDGid())==11 ) i_type = 1;
    else if( std::abs(i.getLeg1().getPDGid())==13 ) i_type = 0;

    if( std::abs(j.getLeg1().getPDGid())==11 ) j_type = 1;
    else if( std::abs(j.getLeg1().getPDGid())==13 ) j_type = 0;

    if (i_type<j_type) return true;
    else if(i_type>j_type) return false;

    //step 1: two pairs with same leg1 lepton: leg 1 ISO
    i_iso = std::abs(i.getLeg1().getPDGid())==15 ? -i.getLeg1().getProperty(HTTEvent::usePropertyFor.at("tauIsolation")) :
            std::abs(i.getLeg1().getPDGid())==11 ? i.getLeg1().getProperty(HTTEvent::usePropertyFor.at("electronIsolation")):
            i.getLeg1().getProperty(HTTEvent::usePropertyFor.at("muonIsolation"));
    if(i_iso<-1) i_iso=999; //something went wrong

    j_iso = std::abs(j.getLeg1().getPDGid())==15 ? -j.getLeg1().getProperty(HTTEvent::usePropertyFor.at("tauIsolation")) :
            std::abs(j.getLeg1().getPDGid())==11 ? j.getLeg1().getProperty(HTTEvent::usePropertyFor.at("electronIsolation")):
            j.getLeg1().getProperty(HTTEvent::usePropertyFor.at("muonIsolation"));
    if(j_iso<-1) j_iso=999; //something went wrong

    if (i_iso<j_iso) return true;
    else if(i_iso>j_iso) return false;

    //step 2: two pairs with same leg1 lepton that have same isolation: leg1 pt
    if(i.getLeg1().getP4().Pt()>j.getLeg1().getP4().Pt()) return true;
    else if(i.getLeg1().getP4().Pt()<j.getLeg1().getP4().Pt()) return false;

    ////////////////////////////////////// Leg 2 ////////////////////////////////////////////////
    //step 2.5, leg 2 type
    if( std::abs(i.getLeg2().getPDGid())==11 ) i_type = 1;
    else if( std::abs(i.getLeg2().getPDGid())==13 ) i_type = 0;

    if( std::abs(j.getLeg2().getPDGid())==11 ) j_type = 1;
    else if( std::abs(j.getLeg2().getPDGid())==13 ) j_type = 0;

    if (i_type<j_type) return true;
    else if(i_type>j_type) return false;

    //step 3: two pairs with same leg2 lepton: leg 2 ISO
    i_iso = std::abs(i.getLeg2().getPDGid())==15 ? -i.getLeg2().getProperty(HTTEvent::usePropertyFor.at("tauIsolation")) :
            std::abs(i.getLeg2().getPDGid())==11 ? i.getLeg2().getProperty(HTTEvent::usePropertyFor.at("electronIsolation")):
            i.getLeg2().getProperty(HTTEvent::usePropertyFor.at("muonIsolation"));
    if(i_iso<-1) i_iso=999; //something went wrong

    j_iso = std::abs(j.getLeg2().getPDGid())==15 ? -j.getLeg2().getProperty(HTTEvent::usePropertyFor.at("tauIsolation")) :
            std::abs(j.getLeg2().getPDGid())==11 ? j.getLeg2().getProperty(HTTEvent::usePropertyFor.at("electronIsolation")):
            j.getLeg2().getProperty(HTTEvent::usePropertyFor.at("muonIsolation"));
    if(j_iso<-1) j_iso=999; //something went wrong

    if (i_iso<j_iso) return true;
    else if(i_iso>j_iso) return false;

    //step 4: two pairs with same leg2 lepton that have same isolation: leg2 pt
    if(i.getLeg2().getP4().Pt()>j.getLeg2().getP4().Pt()) return true;

    return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// int HTauTauTreeFromNanoBase::isGenPartDaughterPdgId(int index, unsigned int absPdgId){
//   if(index>nGenPart || index<0 || GenPart_genPartIdxMother[index]<0) return -1;
//   if(std::abs(GenPart_pdgId[GenPart_genPartIdxMother[index]])==absPdgId)
//     return GenPart_genPartIdxMother[index];
//   else
//     return isGenPartDaughterPdgId(GenPart_genPartIdxMother[index],absPdgId);
// }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// bool HTauTauTreeFromNanoBase::isGenPartDaughterIdx(int index, int mother){
//   if(index>nGenPart || index<0 || GenPart_genPartIdxMother[index]<0) return false;
//   if(GenPart_genPartIdxMother[index]==mother)
//     return true;
//   else
//     return isGenPartDaughterIdx(GenPart_genPartIdxMother[index],mother);
// }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::getDirectDaughterIndexes(std::vector<unsigned int> &indexes, unsigned int motherIndex, bool ignoreNeutrinos)
{  
    indexes.clear();
    bool isFinal=true;
    for(unsigned int iDau=0;iDau<nGenPart;iDau++)
    {
        if(GenPart_genPartIdxMother[iDau]==(int)motherIndex)
        {
            int aPdgId = std::abs(GenPart_pdgId[iDau]);
            if(aPdgId==std::abs(GenPart_pdgId[motherIndex])){
                isFinal=false;
            }
            if((aPdgId==12 || aPdgId==14 || aPdgId==16) && ignoreNeutrinos) continue;
            indexes.push_back(iDau);
        }
    }
    return isFinal;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauTauTreeFromNanoBase::findFinalCopy(unsigned int index)
{  
    std::vector<unsigned int> daughterIndexes;
    if(getDirectDaughterIndexes(daughterIndexes,index)) return index;

    unsigned int newIndex=index;
    for(unsigned int iDau=0;iDau<daughterIndexes.size();iDau++)
    {
        int pdgId = GenPart_pdgId[daughterIndexes[iDau]];
        if(pdgId==GenPart_pdgId[index]){
            newIndex=daughterIndexes[iDau];
            break;
        }
    }
    //check if found copy is decaying (should not happen for unstable particles)
    getDirectDaughterIndexes(daughterIndexes,newIndex,false);
    if(daughterIndexes.empty()) return index;

    return findFinalCopy(newIndex);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauTauTreeFromNanoBase::findFirstCopy(unsigned int index)
{
    int mo_idx = GenPart_genPartIdxMother[index];
    if(mo_idx<0) return index;

    int mo_pdgId = GenPart_pdgId[mo_idx];
    if(mo_pdgId!=GenPart_pdgId[index]) return index;
    
    return findFirstCopy((unsigned int)mo_idx);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::genTauDecayMode(std::vector<unsigned int> &daughterIndexes)
{

    if(daughterIndexes.empty()) return -99;
   
    int numElectrons      = 0;
    int numMuons          = 0;
    int numChargedPions   = 0;
    int numNeutralPions   = 0;
    int numPhotons        = 0;
    int numNeutrinos      = 0;
    int numOtherParticles = 0;
    
    for(unsigned int idx=0; idx<daughterIndexes.size(); idx++) 
    {
        int pdg_id = std::abs(GenPart_pdgId[daughterIndexes[idx]]);
        if(pdg_id == 11) numElectrons++;
        else if(pdg_id == 13) numMuons++;
        else if(pdg_id == 211 || pdg_id == 321 ) numChargedPions++; //Count both pi+ and K+
        else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 || pdg_id == 311 ) numNeutralPions++; //Count both pi0 and K0 {_L/S}
        else if(pdg_id == 12 || 
                pdg_id == 14 || 
                pdg_id == 16) {
            numNeutrinos++;
        }
        else if(pdg_id == 22) numPhotons++;
        else {
            numOtherParticles++;
        }
    }
    if(numElectrons>1)//sometimes there are gamma->ee conversions
    { 
        numPhotons += numElectrons/2;
        numElectrons -= 2*(numElectrons/2);
    }
    //"convert" photons to piZeros
    numNeutralPions += numPhotons/2;
    numPhotons -= 2*(numPhotons/2);
    int tauDecayMode = HTTAnalysis::tauDecayOther;

    if( numOtherParticles == 0 )
    {
        if( numElectrons == 1 )
        {
          //--- tau decays into electrons
          tauDecayMode = HTTAnalysis::tauDecaysElectron;
        }else if( numMuons == 1 ){
          //--- tau decays into muons
          tauDecayMode = HTTAnalysis::tauDecayMuon;
        }else {
            //--- hadronic tau decays
            switch ( numChargedPions )
            {
              case 1 :
                  if( numPhotons != 0 )
                  {
                      tauDecayMode =  HTTAnalysis::tauDecayOther;
                      break;
                  }
                  switch ( numNeutralPions )
                  {
                      case 0:
                          tauDecayMode = HTTAnalysis::tauDecay1ChargedPion0PiZero;
                          break;
                      case 1:
                          tauDecayMode = HTTAnalysis::tauDecay1ChargedPion1PiZero;
                          break;
                      case 2:
                          tauDecayMode = HTTAnalysis::tauDecay1ChargedPion2PiZero;
                          break;
                      case 3:
                          tauDecayMode = HTTAnalysis::tauDecay1ChargedPion3PiZero;
                          break;
                      case 4:
                          tauDecayMode = HTTAnalysis::tauDecay1ChargedPion4PiZero;
                          break;
                      default:
                          tauDecayMode = HTTAnalysis::tauDecayOther;
                          break;
                  }
                  break;
              case 3 : 
                  if( numPhotons != 0 )
                  {
                      tauDecayMode = HTTAnalysis::tauDecayOther;
                      break;
                  }
                  switch ( numNeutralPions )
                  {
                      case 0 : 
                          tauDecayMode = HTTAnalysis::tauDecay3ChargedPion0PiZero;
                          break;
                      case 1 : 
                          tauDecayMode = HTTAnalysis::tauDecay3ChargedPion1PiZero;
                          break;
                      case 2 : 
                          tauDecayMode = HTTAnalysis::tauDecay3ChargedPion2PiZero;
                          break;
                      case 3 : 
                          tauDecayMode = HTTAnalysis::tauDecay3ChargedPion3PiZero;
                          break;
                      case 4 : 
                          tauDecayMode = HTTAnalysis::tauDecay3ChargedPion4PiZero;
                          break;
                      default:
                          tauDecayMode = HTTAnalysis::tauDecayOther;
                          break;
                  }
                  break;
            }
        }
    }
    return tauDecayMode;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::eventInJson()
{
    //If the jsonVec is empty, then no JSON file was provided so all events should be accepted
    if( jsonVector.empty() ) return true;

    edm::LuminosityBlockID lumiID(run, luminosityBlock);
    for(std::vector<edm::LuminosityBlockRange>::const_iterator itr = jsonVector.begin(); itr != jsonVector.end(); itr++){
        if( edm::contains(*itr,lumiID) ) return true;
    }
    return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::debugWayPoint(string description, vector<double> dbls, vector<int> ints, vector<string> descr)
{
    if (event==check_event_number)
    {   int ind = 0;
        bool what = dbls.size() + ints.size() == descr.size();

        cout << description << ": ";
        for( auto &d : dbls )
        {   
            cout << descr[ind] << "(" << d << ") ";
            if(what) ind++;
        }
        for( auto &i : ints )
        {
            cout << descr[ind] << "(" << i << ") ";
            if(what) ind++;
        }
        cout << endl;
    }
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
