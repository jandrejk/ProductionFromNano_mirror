#!/usr/bin/env python

import os
import sys
import threading

from ROOT import gSystem, TChain, TSystem, TFile, TString, vector, TFileCollection, edm

# from PSet import process
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Config as cms
import json
from glob import glob

def getLumisToRun(JSON):
    if JSON == "": return vector( 'edm::LuminosityBlockRange' )()

    vlumis = vector( 'edm::LuminosityBlockRange' )()
    myList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')
    lumisToProcess = cms.untracked.VLuminosityBlockRange( myList )


    for BlockRange in lumisToProcess:
        Block = BlockRange.split('-')

        startRun =       int(Block[0].split(':')[0])
        startLumiBlock = int(Block[0].split(':')[1])

        if len(Block) > 1:
            endRun =       int(Block[1].split(':')[0])
            endLumiBlock = int(Block[1].split(':')[1])
        else:
            endRun = startRun
            endLumiBlock = endLumiBlock

        vlumis.push_back( edm.LuminosityBlockRange( edm.LuminosityBlockID(startRun, startLumiBlock),
                                                    edm.LuminosityBlockID(endRun, endLumiBlock) ) )
    return vlumis
#######################################################################################################

#dir = "/data/higgs/nanonaod_2016/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/VBFHToTauTau_M125_13TeV_powheg_pythia8/"

aFile =            sys.argv[1]
channel =          sys.argv[2]
systShift =        sys.argv[3]
doSvFit =          int(sys.argv[4])
applyRecoil=       int(sys.argv[5])
nevents =          int(sys.argv[6])

if not "root://" in aFile: aFile = "file://" + aFile
print sys.argv

sync_event=0

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


assert gSystem.CompileMacro('HTTEvent.cxx','k') 
assert gSystem.CompileMacro('syncDATA.C','k')
#status *= gSystem.CompileMacro('NanoEventsSkeleton.C') #RECOMPILE IF IT CHANGES!
assert gSystem.CompileMacro('NanoEventsSkeleton.C','k')

gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisClassicSVfit.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisSVfitTF.so')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libHTT-utilitiesRecoilCorrections.so')

assert gSystem.CompileMacro('HTauTauTreeFromNanoBase.C','k')
from ROOT import HTTParticle, HTTAnalysis

if channel=='mt':
    assert gSystem.CompileMacro('HMuTauhTreeFromNano.C','k')
    from ROOT import HMuTauhTreeFromNano as Ntuplizer

if channel=='et':
    assert gSystem.CompileMacro('HElTauhTreeFromNano.C','k')
    from ROOT import HElTauhTreeFromNano as Ntuplizer

if channel=='tt':
    assert gSystem.CompileMacro('HTauhTauhTreeFromNano.C','k')
    from ROOT import HTauhTauhTreeFromNano as Ntuplizer


JSONfile = ""
# JSONfile = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
vlumis = getLumisToRun(JSONfile)

HTTParticle.corrType = getattr(HTTAnalysis, systShift )

prefix = "-".join([channel, systShift])
Ntuplizer(  aTree,doSvFit,applyRecoil,vlumis, prefix).Loop(nevents,sync_event)


for f in glob('*'):
    if not f in glob(prefix + '*.root') and not f in glob("log*.txt"):
        os.remove(f)


exit(0)
