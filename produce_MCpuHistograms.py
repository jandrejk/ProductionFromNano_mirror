#!/usr/bin/env python

import os
import sys
import threading

from ROOT import gROOT, gSystem, TChain, TSystem, TFile, TString, vector, TFileCollection, edm, TH1F



print 'running python file'
sample = sys.argv[1]

#output_name = sample.replace('/','_').replace('.','_')[31:]
output_name = sample.replace('/','_')[31:]
output_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/MCpuHistograms/"+output_name
#print output_name
#print output_dir


gROOT.SetBatch(True) #used to prevent from default canvas c1 popping out
aROOTFile = TFile.Open(sample)
aTree = aROOTFile.Get("Events")
entries = aTree.GetEntries()
f = TFile.Open(output_dir,"recreate")
h1 = TH1F("h1","pile-up",200,0,200)
aTree.Draw("Pileup_nPU>>h1")
	
h1.Write()
f.Close()


exit(0)
