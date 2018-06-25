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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sample', help='Sample to run over', type=str, metavar = 'SAMPLE')
    parser.add_argument('-v', dest='version', help='Version of ntuples.', type=str, metavar = 'VERSION', default = 'v1')
    parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt'], default = 'mt')
    parser.add_argument('-e', dest='shift', help='Uncert shift of energy scale',choices = ['t0u','t1u','t10u','t0d','t1d','t10d','m0u','m1u','m10u','m0d','m1d','m10d','e0u','e1u','e10u','e0d','e1d','e10d'], default = '')
    parser.add_argument('-j', dest='jobs', help='If set to NJOBS > 0: Run NJOBS in parallel on heplx. Otherwise submit to batch.', type=int, default = 0)
    parser.add_argument('-d', dest='debug', help='Debug', action = "store_true")

    args = parser.parse_args()

    print 'Channel:',args.channel

    #2017 sync
    # sample = 'VBFHToTauTau_M125_13TeV_powheg_pythia8'

    SNP = SteerNanoProduction(args.channel, args.shift, args.jobs, args.debug)
    SNP.runOneSample(args.sample)

class SteerNanoProduction():

    def __init__(self, channel, shift, nthreads=6, debug=False):

        self.basedir = os.getcwd()
        self.outdir = '/afs/hephy.at/data/higgs01'
        self.channel = channel
        self.svfit = False
        self.recoil = False

        self.debug = debug

        if debug:
            self.nthreads = 1
            self.nevents = 10001
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

    def runOneSample(self, sample, version="v1"):
        threads = []

        with open("submit_on_batch.sh") as FSO:
            templ = string.Template( FSO.read() )

        runpath = "/".join([self.basedir,"out" ,sample, 'rundir_'+self.channel+'_' ])
        outdir = "/".join([self.outdir, version, sample])

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        print "Submitting: ", sample

        for idx,file in enumerate( self.getFiles(sample) ):
            if self.debug and idx > 0: break

            rundir = runpath+str(idx)

            ##### Create rundir. Overwrite if it already exists
            if not os.path.exists(rundir):
                os.makedirs(rundir)
            else:
                shutil.rmtree(rundir)
                os.makedirs(rundir)

            os.chdir(self.basedir)
            self.prepareRunDir(rundir)

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
                                             file = file,
                                             channel = self.channel,
                                             systShift = self.systShift,
                                             svfit = int(self.svfit),
                                             recoil = int(self.recoil),
                                             nevents = int(self.nevents) )

                os.chdir(rundir )
                with open("submit.sh","w") as FSO:
                    FSO.write( runscript )
                os.system( "sbatch submit.sh" )
                sleep(10)

        if self.nthreads > 0:
            for x in threads:
                x.join()

            if not self.debug:
                os.chdir( "/".join([ self.basedir, "out" ,sample ]) )
                os.system('hadd -f -O '+'ntuple_{0}.root rundir_{0}_*/{1}_*root'.format(self.channel, "-".join([self.channel, self.systShift]) ) )

    def runOneFileLocal(self, rundir ,file):


        os.chdir(rundir )
        njob = int(rundir.split("_")[-1]) + 1
        # Run local

        print "\033[93mrun\033[0m  Job {0}: ".format(njob) + file.split("/")[-1]

        runcmd  = './convertNanoParallel.py {0} {1} {2} {3} {4} {5}'.format( file,
                                                                         self.channel,
                                                                         self.systShift,
                                                                         int(self.svfit),
                                                                         int(self.recoil),
                                                                         int(self.nevents) )
        if self.debug:
                p = sp.Popen(shlex.split( runcmd ),
                             stdout = sys.__stdout__,
                             stderr = sys.__stderr__,
                             shell=False)

                p.communicate()
        else:
            with open("log.txt", 'w') as log:
                p = sp.Popen(shlex.split( runcmd ),
                             stdout = log,
                             stderr = log,
                             shell=False)

                p.communicate()


        print "\033[92mdone\033[0m Job {0}: ".format(njob)  + file.split("/")[-1]


    def prepareRunDir(self, rundir):

        headerfiles = glob("*.h")
        Cfiles = glob("*.cxx") + glob("*.C") + glob("*.cc")
        addFiles = glob("zpt*root")
        addFiles.append('PSet.py')
        addFiles.append('convertNanoParallel.py')
        addFiles.append('Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')


        for f in headerfiles + Cfiles + addFiles:
            shutil.copyfile("/".join([self.basedir,f]), "/".join([rundir,f]) )

            os.chmod("/".join([rundir,f]), 0777)


    def getFiles(self,sample):

        with open( "{0}/samples/{1}.txt".format(self.basedir, sample) ) as FSO:
            buf = FSO.read()
        return buf.splitlines()


if __name__ == '__main__':
  main()

