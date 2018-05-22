#!/usr/bin/env python

import sys, os
import threading
from time import sleep

channel='mt'

nthreads = 6
#dir = '/afs/hephy.at/work/m/mflechl/cmssw/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/'
dir = os.getcwd()+'/'

def runOneFile(idx,file):
#    os.system('cp -p *so *pcm '+dir+'rundir_'+channel+'_'+str(idx))
    print str(idx),file,dir+'rundir_'+channel+'_'+str(idx)
    os.system('cp -p *h *cxx *C *cc PSet.py zpt*root '+dir+'rundir_'+channel+'_'+str(idx))
    os.system('cp -p '+dir+'convertNanoParallel.py '+dir+'rundir_'+channel+'_'+str(idx))
    os.chdir(dir+'rundir_'+channel+'_'+str(idx))
    os.system('rm -f HTT*root')
    os.system('./convertNanoParallel.py '+channel+' '+file+' &>log.txt')
    os.chdir(dir)
    print str(idx)+' done'

print 'Channel:',channel

#2016 sync
#dirName = '/data/higgs/nanoaod_2016/PUMoriond17_05Feb2018_94X_mcRun2_asymptotic_v2-v1/VBFHToTauTau_M125_13TeV_powheg_pythia8/'
#fileNames = [
#    "DEBF5F61-CC12-E811-B47A-0CC47AA9943A.root",
#    "5A038C2A-CC12-E811-B729-7845C4FC3B8D.root",
#    "0E6F4B78-CC12-E811-B37D-FA163EA12C78.root",
#    "50BE09DD-CC12-E811-869D-F04DA27542B9.root",
#    "844BE355-CD12-E811-8871-FA163ED9B872.root",
#]

#2017 sync
dirName = '/data/higgs/nanoaod_2017/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/'
fileNames = [
    "4E3EF595-DD43-E811-834A-002590E7DFEE.root"
]
#fileNames = [
#    "107D8998-A543-E811-8F6A-0CC47A4D7628.root",
#    "1E41103A-E743-E811-AC30-A0369FD0B122.root",
#    "42CA2AF3-3E43-E811-8C57-00248C55CC3C.root",
#    "4E3EF595-DD43-E811-834A-002590E7DFEE.root",
#    "4E84E643-DF43-E811-A2DB-0025905B857A.root",
#    "52DEBA60-9043-E811-A4ED-002590E7DFFE.root",
#    "54B7F462-DD43-E811-A4D4-48D539D33331.root",
#    "88E41764-DF43-E811-BB8E-002590E7D7DE.root",
#    "A433D7DB-DD43-E811-BEFC-0025905A6056.root",
#    "BA8DF356-9043-E811-BF25-A0369FD0B192.root",
#    "E0B523DB-4C43-E811-A9D9-002590E7DE26.root",
#    "FA4F6F5B-9043-E811-8DB2-0CC47A4D76C8.root",
#]

#commented out: when rerunning, will only recompile what is needed.
#os.system('rm -rf rundir_'+channel+'_*')

threads = []
ctr=0
for idx,file in enumerate(fileNames):
    if not os.path.exists(dir+'rundir_'+channel+'_'+str(idx)):
        os.makedirs(dir+'rundir_'+channel+'_'+str(idx))
    t = threading.Thread(target=runOneFile, args=(idx,dirName+file,) )
    threads += [t]

    while threading.active_count()>nthreads:  #pause until thread slots become available
        pass

    #avoid problems when all start at the same time
    sleep(5)

    t.start()


for x in threads:
    x.join()

#hadd -O ntuples/v15_mt.root rundir_'+channel+'_*/HTTM*root
if channel=='et' or channel=='all': os.system('hadd -f -O ntuple_et.root rundir_'+channel+'_*/HTTE*root')
if channel=='mt' or channel=='all': os.system('hadd -f -O ntuple_mt.root rundir_'+channel+'_*/HTTM*root')
if channel=='tt' or channel=='all': os.system('hadd -f -O ntuple_tt.root rundir_'+channel+'_*/HTTT*root')
