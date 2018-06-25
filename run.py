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
    parser.add_argument('-j', dest='jobs', help='If set to NJOBS > 0: Run NJOBS in parallel on heplx. Otherwise submit to batch.', type=int, default = 0)
    parser.add_argument('-d', dest='debug', help='Debug', action = "store_true")

    args = parser.parse_args()

    print 'Channel:',args.channel

    #2017 sync
    # sample = 'VBFHToTauTau_M125_13TeV_powheg_pythia8'

    SNP = SteerNanoProduction(args.channel, args.jobs, args.debug)
    SNP.runOneSample(args.sample)

class SteerNanoProduction():

    def __init__(self, channel, nthreads=6, debug=False):

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

            if not os.path.exists(rundir):
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
                os.system('hadd -f -O '+'ntuple_{0}.root rundir_{0}_*/HTT{1}*root'.format(self.channel, self.channel.upper() ) )

    def runOneFileLocal(self, rundir ,file):


        os.chdir(rundir )
        njob = int(rundir.split("_")[-1]) + 1
        # Run local

        print "\033[93mrun\033[0m  Job {0}: ".format(njob) + file.split("/")[-1]

        runcmd  = './convertNanoParallel.py {0} {1} {2} {3} {4}'.format( file,
                                                                         self.channel,                                                                   
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

