import subprocess as sp
import shlex
import json
import shutil
import os 
from runUtils import checkProxy
import getFilelists as getFL
import sys, os
from glob import glob

# create the diretory merged within MCpuHistograms (if it does not exist)
if not os.path.exists("./MCpuHistograms/merged") :
    os.makedirs("./MCpuHistograms/merged")


with open("sample_collection.json","r") as FSO:
    config = json.load(FSO)

base_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/"
# get the list of all available MC histograms to be merged
MC_hist_list = [mc[112:] for mc in glob(base_dir+'/MCpuHistograms/*') if '.root' in mc]


# use only MC samples as this is used to produce a list 
# that will be used later for MC PU distribution merging
frmt = 'mc'
f = config[frmt]
r = f["run"]["RunIIFall17NanoAOD"]

	
no_of_available_histos = len(MC_hist_list)
no_of_merged_histos = 0
counter = 0
for sample in r["samples"]:
    #print sample
    
    links = r["samples"][sample]
    for link in links:
	#print link
	first_parser = link.split('/')[1]
    	if 'ext' in link :
	    hadd_list = [mc for mc in MC_hist_list if 'ext' in mc and first_parser in mc]
	    #print link
	    #print first_parser
	    #print len(hadd_list)
    	    #print '-'*30
        else :
	    hadd_list = [mc for mc in MC_hist_list if 'ext' not in mc and first_parser in mc]
	no_of_merged_histos += len(hadd_list)
    	hadd_list_withPath = ['MCpuHistograms/'+s for s in hadd_list ]
	
        hadd_cmd = 'hadd ./MCpuHistograms/merged/' + link.replace('/','_')+'.root ' +' '.join(hadd_list_withPath)
	if link == '/DYJetsToLL_M-5to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/NANOAODSIM' :
	     print hadd_cmd			
	#print hadd_cmd
	#os.system(hadd_cmd)
        counter += 1
	#print link
print 'number of output files ', counter
if (no_of_merged_histos - no_of_available_histos) == 0 :
    print 'merging was successful!'
else :
    print 'not all files were merged, something went wrong'
