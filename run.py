#!/usr/bin/env python

import sys, os
import shutil
import threading
from glob import glob
from time import sleep
import subprocess as sp
import multiprocessing as mp
import shlex
import argparse
import string
import json

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sample', help='Sample to run over', type=str, metavar = 'SAMPLE', default = "")
    parser.add_argument('-v', dest='version', help='Version of ntuples.', type=str, metavar = 'VERSION', default = 'v1')
    parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt'], default = 'mt')
    parser.add_argument('-e', dest='shift', help='Uncert shift of energy scale',choices = ['t0u','t1u','t10u','t0d','t1d','t10d','m0u','m1u','m10u','m0d','m1d','m10d','e0u','e1u','e10u','e0d','e1d','e10d'], default = '')
    parser.add_argument('-j', dest='jobs', help='If set to NJOBS > 0: Run NJOBS in parallel on heplx. Otherwise submit to batch.', type=int, default = 0)
    parser.add_argument('-o', dest='outdir', help='Where to write output when running on batch.', type=str, default = '/afs/hephy.at/data/higgs01')
    parser.add_argument('-m', dest='merge', help='Merge sample in outputfolder of batch', action = "store_true")
    parser.add_argument('-d', dest='debug', help='Debug', action = "store_true")
    parser.add_argument('--cert', dest='cert', help='Cert when running over data.', type=str, default =  'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')
    parser.add_argument('--event', dest='event', help='Debug', default = 0)


    args = parser.parse_args()

    #2017 sync
    # sample = 'VBFHToTauTau_M125_13TeV_powheg_pythia8'

    if not args.merge:
        print 'Channel:',args.channel
        SNP = SteerNanoProduction(args.channel, args.shift, args.outdir, args.jobs, args.debug, args.event, args.cert)
        # print SNP.makeConfigBalls(args.sample)
        SNP.runOneSample(args.sample)

    else:
        mergeSample(args.version, args.outdir,args.sample)

class SteerNanoProduction():

    def __init__(self, channel, shift, outdir='', nthreads=6, debug=False, event = 0, cert = ""):

        self.basedir = os.getcwd()
        self.outdir = outdir
        self.channel = channel
        self.svfit = False
        self.recoil = False


        self.debug = debug
        self.event = event

        if debug:
            self.nthreads = 1
            self.nevents = 10001 if not event else -1

        else:
            self.nthreads = nthreads
            self.nevents = -1

        shifts={'t0u' : 'TES1p0p0Up',
                't1u' : 'TES1p1p0Up',
                't10u': 'TES3p0p0Up',
                't0d' : 'TES1p0p0Down',
                't1d' : 'TES1p1p0Down',
                't10d': 'TES3p0p0Down',
                'm0u' : 'MES1p0p0Up',
                'm1u' : 'MES1p1p0Up',
                'm10u': 'MES3p0p0Up',
                'm0d' : 'MES1p0p0Down',
                'm1d' : 'MES1p1p0Down',
                'm10d': 'MES3p0p0Down',
                'e0u' : 'EES1p0p0Up',
                'e1u' : 'EES1p1p0Up',
                'e10u': 'EES3p0p0Up',
                'e0d' : 'EES1p0p0Down',
                'e1d' : 'EES1p1p0Down',
                'e10d': 'EES3p0p0Down'
        }

        if shift:
            self.systShift = shifts.get(shift,'NOMINAL')
        else:
            self.systShift = "NOMINAL"

        self.certJson = cert

    def runOneSample(self, sample, version="v1"):
        threads = []

        assert sample
        sample = sample.split("/")[-1].replace(".txt","")

        with open("submit_on_batch.sh") as FSO:
            templ = string.Template( FSO.read() )

        runpath = "/".join([self.basedir,"out", version ,sample, 'rundir_'+self.channel+'_' ])
        outdir = "/".join([self.outdir, version, sample])
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        print "Submitting: ", sample

        for idx,configBall in enumerate( self.makeConfigBalls(sample) ):
            file = configBall["file"]

            if self.debug and idx > 0: break

            rundir = runpath+str(idx+1)

            ##### Create rundir. Overwrite if it already exists
            if not os.path.exists(rundir):
                os.makedirs(rundir)
            else:
                shutil.rmtree(rundir)
                os.makedirs(rundir)


            os.chdir(self.basedir)
            self.prepareRunDir(rundir)
            self.throwConfigBall(configBall, rundir)

            # Run local
            if self.nthreads > 0:
                t = threading.Thread(target=self.runOneFileLocal, args=(rundir,file,) )
                threads.append(t)

                while threading.active_count()>self.nthreads:  #pause until thread slots become available
                    pass
                t.start()
                sleep(5)
            # Run on batch system
            else:
                runscript = templ.substitute(rundir = rundir,
                                             outdir = outdir,
                                             sleeping = idx*10,
                                             channel = self.channel)

                os.chdir(rundir )
                with open("submit.sh","w") as FSO:
                    FSO.write( runscript )
                os.system( "sbatch submit.sh" )

        if self.nthreads > 0:
            for x in threads:
                x.join()

            if not self.debug:
                os.chdir( "/".join([ self.basedir, "out", version ,sample ]) )
                os.system('hadd -f -O '+'{0}_all.root rundir_{1}_*/{0}_*root'.format("-".join([self.channel, self.systShift]), self.channel ) )

    def runOneFileLocal(self, rundir ,file):


        os.chdir(rundir )
        njob = int(rundir.split("_")[-1])
        # Run local

        print "\033[93mrun\033[0m  Job {0}: ".format(njob) + file.split("/")[-1]

        if self.debug:
                p = sp.Popen(shlex.split( './convertNanoParallel.py' ),
                             stdout = sys.__stdout__,
                             stderr = sys.__stderr__,
                             shell=False)

                p.communicate()
        else:
            with open("log.txt", 'w') as log:
                p = sp.Popen(shlex.split( './convertNanoParallel.py' ),
                             stdout = log,
                             stderr = log,
                             shell=False)

                p.communicate()


        print "\033[92mdone\033[0m Job {0}: ".format(njob)  + file.split("/")[-1]


    def prepareRunDir(self, rundir):

        headerfiles = glob("*.h")
        Cfiles = glob("*.cxx") + glob("*.C") + glob("*.cc")
        addFiles =['convertNanoParallel.py']

        shutil.copytree("utils", "/".join([rundir,"utils"]))
        for f in headerfiles + Cfiles + addFiles:
            shutil.copyfile("/".join([self.basedir,f]), "/".join([rundir,f]) )

            os.chmod("/".join([rundir,f]), 0777)

    def makeConfigBalls(self,sample):
        configBalls = []
        samples_avail = glob("samples/*/*/*")
        for sa in samples_avail:
            if sample in sa:
                files = self.getFiles(sa)
                parts = sa.split("/")

        for file in files:
            configBall = {}
            configBall["file"]        = file
            configBall["format"]      = parts[1]
            configBall["sample"]      = parts[2]
            configBall["channel"]     = self.channel
            configBall["systShift"]   = self.systShift
            configBall["svfit"]       = int(self.svfit)
            configBall["recoil"]      = int(self.recoil)
            configBall["nevents"]     = int(self.nevents)
            configBall["check_event"] = int(self.event)

            if parts[1] == "data":
                configBall["certJson"] = self.certJson
            else:
                configBall["certJson"] = ""

            configBalls.append(configBall)

        return configBalls

    def throwConfigBall(self, ball, where):

        with open("/".join([where, "configBall.json"]),"w") as FSO:
            json.dump(ball, FSO, indent=4)


    def getFiles(self,sample):

        with open( "{0}/{1}".format(self.basedir, sample) ) as FSO:
            buf = FSO.read()
        return buf.splitlines()


def mergeSample(version, outdir, single_sample=""):

    import ROOT as R

    path = "/".join([ outdir, version])

    if single_sample: samples = [ "/".join([path, single_sample.split("/")[-1].replace(".txt","")] ) ]
    else: samples = glob(path + "/*")


    for sample in samples:
        files = glob( sample + "/*" )

        types = {}
        for f in files:
            t = f.split("/")[-1].split("_")[0]
            if "_all.root" in f or "rundir" in f: continue

            if not t in types:
                types[t] = {"chain":R.TChain("TauCheck"), "N":1}
                types[t]["chain"].Add(f)
            else:
                types[t]["chain"].Add(f)
                types[t]["N"] += 1

        for t in types:
            print "Merging {0} for sample {1}: total {2} files".format(t, sample.split("/")[-1], types[t]["N"] )

            types[t]["chain"].Merge( "/".join([sample, t+"_all.root"]) )

    



if __name__ == '__main__':
  main()

