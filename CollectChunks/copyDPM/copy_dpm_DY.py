import os

import sys
directory = sys.argv[1]

names = [
'DY1JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v1',
'DY1JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_v3_94X_mc2017_realistic_v14_ext1_v2',
'DY2JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1_v1',
'DY2JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'DY3JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1_v1',
'DY3JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'DY4JetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_v2_94X_mc2017_realistic_v14_v2',
'DYJetsToLL_M_10to50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1_v1',
'DYJetsToLL_M_50_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_v1'
]

for name in names :
	cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/dy/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/dy/{0}'.format(name,directory)
	print cmd1
	print cmd2
	os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
