#!/usr/bin/env python

import sys, os
import threading

channel='tt'

nthreads = 6
dir = '/afs/hephy.at/work/m/mflechl/cmssw/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/'

def runOneFile(idx,file):
#    os.system('cp -p *so *pcm '+dir+'rundir_'+str(idx))
    os.system('cp -p *h *cxx *C *cc PSet.py zpt*root '+dir+'rundir_'+str(idx))
    os.system('cp -p '+dir+'convertNanoParallel.py '+dir+'rundir_'+str(idx))
    os.chdir(dir+'rundir_'+str(idx))
    os.system('./convertNanoParallel.py '+channel+' '+file+' &>log.txt')
    os.chdir(dir)

fileNames = [
    "DEBF5F61-CC12-E811-B47A-0CC47AA9943A.root",
    "5A038C2A-CC12-E811-B729-7845C4FC3B8D.root",
    "0E6F4B78-CC12-E811-B37D-FA163EA12C78.root",
    "50BE09DD-CC12-E811-869D-F04DA27542B9.root",
    "844BE355-CD12-E811-8871-FA163ED9B872.root",
]                                

#os.system('rm -rf rundir_*')
threads = []
ctr=0
for idx,file in enumerate(fileNames):
    if not os.path.exists(dir+'rundir_'+str(idx)):
        os.makedirs(dir+'rundir_'+str(idx))
    t = threading.Thread(target=runOneFile, args=(idx,file,) )
    threads += [t]

    while threading.active_count()>nthreads:  #pause until thread slots become available
        pass

    t.start()


for x in threads:
    x.join()

#hadd -O ntuples/v15_mt.root rundir_*/HTTM*root
if channel=='mt' or channel=='all': os.system('hadd -O ntuple_mt.root rundir_*/HTTM*root')
if channel=='tt' or channel=='all': os.system('hadd -O ntuple_tt.root rundir_*/HTTT*root')
