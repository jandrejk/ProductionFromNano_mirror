import sys, os
from glob import glob

base_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/"


#file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]
#MC_hist_list = [mc[112:] for mc in glob(base_dir+'/MCpuHistograms/*')]
#print len (MC_hist_list)



# lsit(set()) structure to remove dublicates in case there are some ext for some samples
#file_names = list(set(['_'.join(f.split('/')[-1].split('_')[0:3]) for f in file_list]))
#file_names= file_names[1:]



    
if not os.path.exists("./MCpuHistograms/merged") :
    os.makedirs("./MCpuHistograms/merged")
if not os.path.exists("./MCpuHistograms/usedHistos") :
    os.makedirs("./MCpuHistograms/usedHistos")
"""
counter = 0
for name in file_names :
    save_name = './MCpuHistograms/merged/'+name+'.root'
    
    #print 'hadd '+save_name+' MCpuHistograms/store_mc_RunIIFall17NanoAOD_'+name+'*'
    command = 'hadd '+save_name+' MCpuHistograms/store_mc_RunIIFall17NanoAOD_'+name+'*'
    #os.system(command)
print 'merging completed'
"""

this_list = [
    ("WW", "WW"),
    ("WZ", "WZ"),
    ("ZZ", "ZZ"),
    ("DY1JetsToLL_ext1", "DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_"),
    ("DY2JetsToLL_ext1", "DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_"),
    ("DY3JetsToLL_ext1", "DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_"),
    ("DYJetsToLL_ext1" , "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1_"),
    ("DY1JetsToLL"      , "DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ("DY2JetsToLL"      , "DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ("DY3JetsToLL"      , "DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ("DY4JetsToLL"      , "DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ("DY5JetsToLL"      , "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ("EWKWMinus2Jets_WToLNu" , "EWKWMinus2Jets_WToLNu"),
    ("EWKWPlus2Jets_WToLNu"  , "EWKWPlus2Jets_WToLNu"),
    ("EWKZ2Jets_ZToLL"       , "EWKZ2Jets_ZToLL"),
    ("EWKZ2Jets_ZToNuNu"     , "EWKZ2Jets_ZToNuNu"),
    ("VBFHToTauTau"    , "VBFHToTauTau"),
    ("GluGluHToTauTau" , "GluGluHToTauTau"),
    ("ST_t-channel_antitop_4f" , "ST_t-channel_antitop_4f"),
    ("ST_t-channel_top_4f"     , "ST_t-channel_top_4f"),
    ("ST_tW_antitop_5f"        , "ST_tW_antitop_5f"),
    ("ST_tW_top_5f"            , "ST_tW_top_5f"),
    ("TTTo2L2Nu"        , "TTTo2L2Nu"),
    ("TTToHadronic"     , "TTToHadronic"),
    ("TTToSemiLeptonic" , "TTToSemiLeptonic"),
    ("WJetsToLNu_ext1" , "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2"),
    ("WJetsToLNu"      , "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_NANOAODSIM_PU2017_12Apr2018_94X_mc2017_realistic_v14"),
    ("W1JetsToLNu"     ,  "W1JetsToLNu"),
    ("W2JetsToLNu"     , "W2JetsToLNu"),
    ("W3JetsToLNu"     , "W3JetsToLNu"),
    ("W4JetsToLNu"     , "W4JetsToLNu"),
]


for it in this_list :
    key = it[0]
    value = it[1]
    #print key	
    save_name = './MCpuHistograms/merged/'+key+'.root'
    cmd = 'hadd '+save_name+' MCpuHistograms/store_mc_RunIIFall17NanoAOD_'+value+'*'
    mvcmd = 'mv ./MCpuHistograms/*'+value+'*root ./MCpuHistograms/usedHistos/'
    #print mvcmd
    os.system(cmd)
    os.system(mvcmd)






