#!/usr/bin/env python
import os
import sys
import ROOT as R

"""
This program takes the url of a sample and reads out the "Pileup_nTrueInt" branch of that 
sample. The output is a root file named after the DAS query of the sample where / is replaced
by #. The output contains a histogram named h1 showing the pile-up distribution.
"""

print 'running python file'
sample = sys.argv[1]
print sample
output_name = sample.replace('/','_')[31:]
output_dir = "{CMSSW_BASE}/src/WawTools/NanoAODTools/MCpuHistograms/".format(**os.environ)+ output_name
print output_dir

R.gROOT.SetBatch(True) #used to prevent from default canvas c1 popping out
aROOTFile = R.TFile.Open(sample)
aTree = aROOTFile.Get("Events")
entries = aTree.GetEntries()
print entries
f = R.TFile.Open(output_dir,"recreate")
h1 = R.TH1F("h1","pile-up",200,0,200)
aTree.Draw("Pileup_nTrueInt>>h1")
	
h1.Write()
f.Close()
exit(0)
