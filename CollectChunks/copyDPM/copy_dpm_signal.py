import os

import sys
directory = sys.argv[1]

names = [
#'GluGluHToTauTau_M125_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1_v1',
#'GluGluHToTauTau_M125_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
#'VBFHToTauTau_M125_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
#'WminusHToTauTau_M125_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
#'WplusHToTauTau_M125_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
#'ZHToTauTau_M125_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1',
'VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v1'
]

for name in names :
	#cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/signal/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/signal/{0}'.format(name,directory)
	#print cmd1
	print cmd2
	#os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
