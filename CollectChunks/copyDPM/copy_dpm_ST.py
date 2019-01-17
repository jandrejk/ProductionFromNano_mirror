import os

import sys
directory = sys.argv[1]

names = [
'ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v2',
'ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV_powheg_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v2',
'ST_t_channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8_PU2017_12Apr2018_94X_mc2017_realistic_v14_v2',
'ST_t_channel_top_4f_inclusiveDecays_TuneCP5_13TeV_powhegV2_madspin_pythia8_PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_v1'
]

for name in names :
	cmd1 = 'gfal-mkdir srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/st/{0}'.format(name)
	cmd2 = 'gfal-copy -r -t 100000 -n 8 {1}/{0} srm://hephyse.oeaw.ac.at//dpm/oeaw.ac.at/home/cms/store/user/jaandrej/nanoAOD/20190107/samples/mc/st/{0}'.format(name,directory)
	print cmd1
	print cmd2
	os.system(cmd1)
	os.system(cmd2)
print 'copying to dpm finished'
