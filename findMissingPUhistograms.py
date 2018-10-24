import sys, os
from glob import glob

base_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/"


file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]
root_files = []

for file_name in file_list :
    for sample in open(file_name,'r') :
        sample = sample.replace('/','_')
        root_files.append(sample)

MC_hist_list = [mc for mc in glob(base_dir+'/MCpuHistograms/*')]

for m in MC_hist_list :
    if m not in root_files :
        print('The following histogram is missing ', m)