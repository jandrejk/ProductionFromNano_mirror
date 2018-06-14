#!/usr/bin/env python

import os
import sys
import threading

from ROOT import gSystem, TChain, TSystem, TFile, TString, vector, TFileCollection

from PSet import process

#dir = "/data/higgs/nanonaod_2016/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/VBFHToTauTau_M125_13TeV_powheg_pythia8/"
channel =          sys.argv[1]
aFile =            sys.argv[2]
doSvFit =          int(sys.argv[3])
applyRecoil=       int(sys.argv[4])

if not "root://" in aFile: aFile = "file://" + aFile
print sys.argv
#sync_event=850381
sync_event=0
#doSvFit = True

#applyRecoil=False
nevents=-1      #all
#nevents=5000
vlumis = vector('string')()

print 'Channel: ',channel

if doSvFit :
    print "Run with SVFit computation"
if applyRecoil :
    print "Apply MET recoil corrections"

print "Using file: ",aFile
aROOTFile = TFile.Open(aFile)
aTree = aROOTFile.Get("Events")

print "TTree entries: ",aTree.GetEntries()
print "Compiling...."

#Some system have problem runnig compilation (missing glibc-static library?).
#First we try to compile, and only then we start time consuming cmssw
status = 1
gSystem.CompileMacro('HTTEvent.cxx','k')
status *= gSystem.CompileMacro('syncDATA.C','k')
#status *= gSystem.CompileMacro('NanoEventsSkeleton.C') #RECOMPILE IF IT CHANGES!
status *= gSystem.CompileMacro('NanoEventsSkeleton.C','k')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisClassicSVfit.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisSVfitTF.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHTT-utilitiesRecoilCorrections.so')


status *= gSystem.CompileMacro('HTauTauTreeFromNanoBase.C','k')
if channel=='mt' or channel=='all': status *= gSystem.CompileMacro('HMuTauhTreeFromNano.C','k')
if channel=='et' or channel=='all': status *= gSystem.CompileMacro('HElTauhTreeFromNano.C','k')
if channel=='tt' or channel=='all': status *= gSystem.CompileMacro('HTauhTauhTreeFromNano.C','k')

print "Compilation status: ",status
if status==0:
    exit(-1)

if channel=='mt' or channel=='all': from ROOT import HMuTauhTreeFromNano
if channel=='et' or channel=='all': from ROOT import HElTauhTreeFromNano
if channel=='tt' or channel=='all': from ROOT import HTauhTauhTreeFromNano

lumisToProcess = process.source.lumisToProcess
#import FWCore.ParameterSet.Config as cms
#lumisToProcess = cms.untracked.VLuminosityBlockRange( ("1:2047-1:2047", "1:2048-1:2048", "1:6145-1:6145", "1:4098-1:4098", "1:3-1:7", "1:6152-1:6152", "1:9-1:11", "1:273-1:273", "1:4109-1:4109", "1:4112-1:4112", "1:4115-1:4116") )
for lumi in lumisToProcess:
    vlumis.push_back(lumi)

# for name in fileNames:

    # aFile = "file://"+name



if channel=='mt' or channel=='all': HMuTauhTreeFromNano(  aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)
if channel=='et' or channel=='all': HElTauhTreeFromNano(  aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)
if channel=='tt' or channel=='all': HTauhTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)


exit(0)
