////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//                  Define all needed JEC uncertainty sources here            //
////////////////////////////////////////////////////////////////////////////////
enum class JecUncertEnum {
AbsoluteStat = 0,
AbsoluteScale,
AbsoluteFlavMap,
AbsoluteMPFBias,
Fragmentation,
SinglePionECAL,
SinglePionHCAL,
FlavorQCD,
TimePtEta,
RelativeJEREC1,
RelativeJEREC2,
RelativeJERHF,
RelativePtBB,
RelativePtEC1,
RelativePtEC2,
RelativePtHF,
RelativeBal,
RelativeSample,
RelativeFSR,
RelativeStatFSR,
RelativeStatEC,
RelativeStatHF,
PileUpDataMC,
PileUpPtRef,
PileUpPtBB,
PileUpPtEC1,
PileUpPtEC2,
PileUpPtHF,
Total,
NONE
};

const vector<string> JecUncertNames =
{
"AbsoluteStat",
"AbsoluteScale",
"AbsoluteFlavMap",
"AbsoluteMPFBias",
"Fragmentation",
"SinglePionECAL",
"SinglePionHCAL",
"FlavorQCD",
"TimePtEta",
"RelativeJEREC1",
"RelativeJEREC2",
"RelativeJERHF",
"RelativePtBB",
"RelativePtEC1",
"RelativePtEC2",
"RelativePtHF",
"RelativeBal",  
"RelativeSample",
"RelativeFSR",
"RelativeStatFSR",
"RelativeStatEC",
"RelativeStatHF",
"PileUpDataMC",
"PileUpPtRef",
"PileUpPtBB",
"PileUpPtEC1",
"PileUpPtEC2",
"PileUpPtHF",
"Total",
""
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

const map< string, std::vector<JecUncertEnum> > JecAfterSplitting = {
    { "Total", {JecUncertEnum::Total} },
    {"CMS_scale_j_eta0to5_13Tev", {JecUncertEnum::SinglePionECAL, JecUncertEnum::SinglePionHCAL, JecUncertEnum::AbsoluteFlavMap,
                                   JecUncertEnum::AbsoluteMPFBias, JecUncertEnum::AbsoluteScale, JecUncertEnum::AbsoluteStat, 
                                   JecUncertEnum::Fragmentation, JecUncertEnum::FlavorQCD, JecUncertEnum::TimePtEta,
                                   JecUncertEnum::PileUpDataMC, JecUncertEnum::RelativeFSR, JecUncertEnum::RelativeStatFSR,
                                   JecUncertEnum::PileUpPtRef }
    },

    {"CMS_scale_j_eta0to3_13TeV", {JecUncertEnum::PileUpPtEC1, JecUncertEnum::PileUpPtEC2, JecUncertEnum::PileUpPtBB,
                                   JecUncertEnum::RelativeJEREC1, JecUncertEnum::RelativeJEREC2, JecUncertEnum::RelativePtEC1,
                                   JecUncertEnum::RelativePtEC2, JecUncertEnum::RelativeStatEC, JecUncertEnum::RelativePtBB}
    },

    {"CMS_scale_j_eta3to5_13TeV",{JecUncertEnum::RelativeStatHF, JecUncertEnum::RelativePtHF, JecUncertEnum::PileUpPtHF, JecUncertEnum::RelativeJERHF} },

    {"CMS_scale_j_RelativeBal_13TeV", {JecUncertEnum::RelativeBal} },

    {"CMS_scale_j_RelativeSample_13TeV", {JecUncertEnum::RelativeSample} },
};

