#!/usr/bin/env python
import os
import sys
import threading
from ROOT import gROOT, gSystem, TChain, TSystem, TFile, TString, vector, TFileCollection, edm, TH1F

"""
This program takes the url of a sample and reads out the "Pileup_nTrueInt" branch of that 
sample. The output is a root file named after the DAS query of the sample where / is replaced
by #. The output contains a histogram named h1 showing the pile-up distribution.
"""

print 'running python file'
sample = sys.argv[1]

output_name = sample.replace('/','_')[31:]
output_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/MCpuHistograms/"+output_name

gROOT.SetBatch(True) #used to prevent from default canvas c1 popping out
aROOTFile = TFile.Open(sample)
aTree = aROOTFile.Get("Events")
entries = aTree.GetEntries()
f = TFile.Open(output_dir,"recreate")
h1 = TH1F("h1","pile-up",200,0,200)
aTree.Draw("Pileup_nTrueInt>>h1")
	
h1.Write()
f.Close()


exit(0)
