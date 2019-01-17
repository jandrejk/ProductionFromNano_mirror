import os

import sys
directory = sys.argv[1]

names = [
'EWKWMinus2Jets_WToLNu_M_50_TuneCP5_13TeV_madgraph_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'EWKWPlus2Jets_WToLNu_M_50_TuneCP5_13TeV_madgraph_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'EWKZ2Jets_ZToLL_M_50_TuneCP5_13TeV_madgraph_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'EWKZ2Jets_ZToNuNu_TuneCP5_13TeV_madgraph_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1'
]

for name in names :
	cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/ewk/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/ewk/{0}'.format(name,directory)
	print cmd1
	print cmd2
	os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
