import os

import sys
directory = sys.argv[1]

names = [
'W1JetsToLNu_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v3',
'W2JetsToLNu_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v4',
'W3JetsToLNu_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'W4JetsToLNu_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'WJetsToLNu_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1_v2',
'WJetsToLNu_TuneCP5_13TeV_madgraphMLM_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v2'
]

for name in names :
	cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/w/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/w/{0}'.format(name,directory)
	print cmd1
	print cmd2
	os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
