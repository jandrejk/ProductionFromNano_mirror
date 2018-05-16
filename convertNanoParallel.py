#!/usr/bin/env python

import os
import sys
import threading

from ROOT import gSystem, TChain, TSystem, TFile, TString, vector

from PSet import process

dir = "/data/higgs/nanonaod_2016/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/VBFHToTauTau_M125_13TeV_powheg_pythia8/"
channel = sys.argv[1]
fileNames = [ sys.argv[2] ]

sync_event=1321942
#sync_event=0
#doSvFit = True
doSvFit = False
applyRecoil=True
#applyRecoil=False
nevents=-1      #all
#nevents=5000
vlumis = vector('string')()
nthreads = 6

print 'Channel: ',channel

if doSvFit :
    print "Run with SVFit computation"
if applyRecoil :
    print "Apply MET recoil corrections"

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

# xline=gSystem.GetMakeSharedLib()+' -Wattributes'
# xline=xline.replace(' -W ',' -W -Wattributes ')
# gSystem.SetMakeSharedLib(xline)
# print 'MMM ',gSystem.GetMakeSharedLib()

# yline=gSystem.GetMakeExe()+' -Wattributes'
# yline=yline.replace(' -W ',' -W -Wattributes ')
# gSystem.SetMakeExe(yline)
# print 'NNN ',gSystem.GetMakeExe()

stdout = sys.stdout
sys.stdout = open('/tmp/pstd', 'w')
stderr = sys.stderr
sys.stderr = open('/tmp/perr', 'w')
status *= gSystem.CompileMacro('HTauTauTreeFromNanoBase.C','k')
if channel=='mt' or channel=='all': status *= gSystem.CompileMacro('HMuTauhTreeFromNano.C','k')
if channel=='et' or channel=='all': status *= gSystem.CompileMacro('HElTauhTreeFromNano.C','k')
if channel=='tt' or channel=='all': status *= gSystem.CompileMacro('HTauhTauhTreeFromNano.C','k')
sys.stdout=stdout
sys.stderr=stderr

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

threads = []
ctr=0
for name in fileNames:
#    aFile = "file:///home/mbluj/work/data/NanoAOD/80X_with944/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISummer16NanoAOD_PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/"+name
    aFile = "file://"+dir+name

    print "Using file: ",aFile
    aROOTFile = TFile.Open(aFile)
    aTree = aROOTFile.Get("Events")
    print "TTree entries: ",aTree.GetEntries()
    if channel=='mt' or channel=='all': HMuTauhTreeFromNano(  aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)
    if channel=='et' or channel=='all': HElTauhTreeFromNano(  aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)
    if channel=='tt' or channel=='all': HTauhTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)

#    print 'A',name,threading.active_count()
#    t = threading.Thread(target=runFile, args=(aFile,) )
#    threads += [t]
#    print 'B',name,threading.active_count()
#    while threading.active_count()>nthreads:  #pause until thread slots become available                                                                                
#        pass
#    print 'C',name,threading.active_count()
#    t.start()





#    print "Making the MuTau tree"
#    aROOTFile = TFile.Open(aFile)
#    aTree = aROOTFile.Get("Events")
#    print "TTree entries: ",aTree.GetEntries()
#    HMuTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop(nevents,sync_event)
##    HMuTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop()

#    print "Making the TauTau tree"
#    aROOTFile = TFile.Open(aFile)
#    aTree = aROOTFile.Get("Events")
#    HTauhTauhTreeFromNano(aTree,doSvFit,applyRecoil,vlumis).Loop(nevents)

    if nevents>0: break

exit(0)
