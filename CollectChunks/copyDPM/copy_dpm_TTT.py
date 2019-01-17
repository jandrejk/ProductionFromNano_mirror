import os

import sys
directory = sys.argv[1]

names = [
'TTTo2L2Nu_TuneCP5_13TeV_powheg_pythia8_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v1',
'TTToHadronic_TuneCP5_13TeV_powheg_pythia8_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v2',
'TTToSemiLeptonic_TuneCP5_13TeV_powheg_pythia8_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v1'
]

for name in names :
	cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/tt/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/tt/{0}'.format(name,directory)
	print cmd1
	print cmd2
	os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
