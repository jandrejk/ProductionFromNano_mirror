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

with open("configBall.json","r") as FSO:
    configBall = json.load(FSO)

aFile =            configBall["file"]
channel =          str(configBall["channel"])
systShift =        str(configBall["systShift"])
JSONfile =         str(configBall["certJson"])
doSvFit =          int(configBall["svfit"])
applyRecoil =      int(configBall["recoil"])
nevents =          int(configBall["nevents"])
check_event =      int(configBall["check_event"])
isMC =             0 if configBall["format"] == "data" else 1

if not "root://" in aFile: aFile = "file://" + aFile
print sys.argv

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



vlumis = getLumisToRun(JSONfile)

HTTParticle.corrType = getattr(HTTAnalysis, systShift )

prefix = "-".join([channel, systShift])
Ntuplizer(  aTree,doSvFit,applyRecoil, isMC, vlumis, prefix).Loop(nevents,check_event)


# for f in glob('*'):
#     if not f in glob(prefix + '*.root') and not f in glob("log*.txt"):
#         os.remove(f)


exit(0)
