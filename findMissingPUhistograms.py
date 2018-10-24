import sys, os
from glob import glob

base_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/"

file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]
root_files = []

MC_hist_list = [mc[112:] for mc in glob(base_dir+'/MCpuHistograms/*')]
counter = 0

for file_name in file_list :
    for sample in open(file_name,'r') :
        sample = sample.replace('/','_')[31:].rstrip()
        root_files.append(sample)
	if not any(sample in x for x in MC_hist_list) :
	    print 'missing file \n', sample, '\n'
	    counter += 1
if counter == 0 :
    print 'all Pileup histograms are there!' 
else :
    print counter, 'files are missing' 


