import os

import sys
directory = sys.argv[1]

names = [
'SingleElectron_Run2017B_31Mar2018_v1',
'SingleElectron_Run2017C_31Mar2018_v1',
'SingleElectron_Run2017D_31Mar2018_v1',
'SingleElectron_Run2017E_31Mar2018_v1',
'SingleElectron_Run2017F_31Mar2018_v1'
]

for name in names :
	cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/data/SingleElectron/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/data/SingleElectron/{0}'.format(name,directory)
	print cmd1
	print cmd2
	os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
