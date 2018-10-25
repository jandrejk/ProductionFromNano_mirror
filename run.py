#!/usr/bin/env python

import sys, os
import shutil
import threading
from glob import glob
from time import sleep, time
import subprocess as sp
import multiprocessing as mp
import shlex
import argparse
import string
import json
from runUtils import checkProxy, checkTokens, useToken, getSystem, getHeplxPublicFolder

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='sample', help='Sample to run over', type=str, metavar = 'SAMPLE', default = "")
    parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt','all'], default = 'mt')
    parser.add_argument('-e', dest='shift', help='Uncert shift of energy scale',choices = ['t0u','t1u','t10u','t0d','t1d','t10d','t0','t1','t10','t',
                                                                                           'm0u','m1u','m10u','m0d','m1d','m10d','m0','m1','m10','m',
                                                                                           'e0u','e1u','e10u','e0d','e1d','e10d','e0','e1','e10','e'], default = '')
    parser.add_argument('-t', dest='submit', help='Where to submit the job',choices = ['condor','batch','local'], default = 'local')
    parser.add_argument('-j', dest='jobs', help='If set to NJOBS > 0: Run NJOBS in parallel on heplx. Otherwise submit to batch.', type=int, default = 8)
    parser.add_argument('-o', dest='outdir', help='Where to write output when running on batch.', type=str, default = '/afs/hephy.at/data/higgs01')
    parser.add_argument('-m', dest='merge', help='Merge sample in outputfolder of batch', action = "store_true")
    parser.add_argument('-d', dest='debug', help='Debug', action = "store_true")
    parser.add_argument('-f', dest='force', help="Forces submission to batch when status in submit_log is 'NEW'", action = "store_true")
    parser.add_argument('--cert', dest='cert', help='Cert when running over data.', type=str, default =  'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')
    parser.add_argument('--event', dest='event', help='Debug', default = 0)
    parser.add_argument('--svfit', dest='svfit', help='Calculate svfit mass', action = "store_true")
    parser.add_argument('--sync', dest='sync', help='Produce sync ntuple', action = "store_true")     


    args = parser.parse_args()

    #2017 sync
    # sample = 'VBFHToTauTau_M125_13TeV_powheg_pythia8'
    if not checkProxy(): sys.exit()
    if not checkTokens(): sys.exit()
    if not os.environ.get("CMSSW_BASE", False):
        print "You forgot to source cmssw"
        sys.exit()


    if not args.merge:

        SNP = SteerNanoProduction(args.outdir, args.submit, args.svfit, args.jobs, args.debug, args.sync, args.force, args.event, args.cert)
        for cmd in makeSubmitList(args.sample, args.channel, args.shift):

            print "\033[1m" +  "_"*100 + "\033[0m"
            print 'Channel:',cmd["channel"]
            print 'Sample:', cmd["sample"]
            print 'Shift:', cmd["shift"]
            print "\033[1m" +  "_"*100 + "\033[0m"

            SNP.runOneSample( **cmd )

    else:
        mergeSample(args.version, args.outdir,args.sample)

def makeSubmitList( sample, channel, shift):
    #Oh god this is ugly...

    if channel == "all": channels = ["et","mt","tt"]
    else: channels = [channel]

    if shift == "t": shifts = ['t0u','t1u','t10u','t0d','t1d','t10d']
    elif shift == "m": shifts = ['m0u','m1u','m10u','m0d','m1d','m10d']
    elif shift == "e": shifts = ['e0u','e1u','e10u','e0d','e1d','e10d']
    elif shift == "t0": shifts = ['t0u','t0d']
    elif shift == "m0": shifts = ['m0u','m0d']
    elif shift == "e0": shifts = ['e0u','e0d']
    elif shift == "t1": shifts = ['t1u','t1d']
    elif shift == "m1": shifts = ['m1u','m1d']
    elif shift == "e1": shifts = ['e1u','e1d']
    elif shift == "t10": shifts = ['t10u','t10d']
    elif shift == "m10": shifts = ['m10u','m10d']
    elif shift == "e10": shifts = ['e10u','e10d']
    else: shifts = [shift]

    if os.path.exists(sample): samples = [sample]
    else:
        if sample == "all": samples = glob("samples/*/*/*")
        if sample == "mc": samples = glob("samples/mc/*/*")
        if sample == "data": samples = glob("samples/data/*/*")
        if sample == "dy": samples = glob("samples/mc/dy/*")
        if sample == "diboson": samples = glob("samples/mc/diboson/*")
        if sample == "ewk": samples = glob("samples/mc/ewk/*")
        if sample == "signal": samples = glob("samples/mc/signal/*")
        if sample == "st": samples = glob("samples/mc/st/*")
        if sample == "tt": samples = glob("samples/mc/tt/*")
        if sample == "w": samples = glob("samples/mc/w/*")

    submitlist = []
    for c in channels:
        for s in samples:
            for sh in shifts:
                if "/data/" in s and sh != "": continue
                if c == "et" and not "/mc/" in s and not "samples/data/SingleElectron" in s: continue
                if c == "mt" and not "/mc/" in s and not "samples/data/SingleMuon" in s: continue
                if c == "tt" and not "/mc/" in s and not "samples/data/Tau" in s: continue

                submitlist.append( { "sample":s, "channel":c, "shift":sh } )

    return submitlist



class SteerNanoProduction():

    def __init__(self, outdir='', submit='local', svfit = False, nthreads=8, debug=False, sync = False, force = False, event = 0, cert = ""):

        self.basedir = os.getcwd()
        self.outdir = outdir
        self.recoil = False
        self.svfit = svfit

        cell = getSystem()
        if submit == "batch":
            if cell == "cern.ch": self.submit = "lxplus"
            if cell == "hephy.at": self.submit = "hephybatch"
        elif submit == "condor":
            if cell == "cern.ch": self.submit = "condor"
            if cell == "hephy.at":
                print "Sorry! Condor not available on heplx. Submitting to batch."
                self.submit = "hephybatch"

        else: self.submit = submit

        self.debug = debug
        self.sync = sync
        self.force = force
        self.event = event

        if debug:
            self.nthreads = 1
            self.nevents = 10001 if not event else -1
        else:
            self.nthreads = nthreads
            self.nevents = -1

        self.shifts={'t0u' : 'TES1p0p0Up',
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



        self.certJson = cert

    def runOneSample(self, sample, channel, shift):
        useToken("hephy")
        threads = []
        os.chdir(self.basedir)

        self.channel = channel
        if shift:
            self.systShift = self.shifts.get(shift,'NOMINAL')
        else:
            self.systShift = "NOMINAL"        

        assert sample
        sample = sample.split("/")[-1].replace(".txt","")

        runpath = "/".join([self.basedir,"out", sample, 'rundir_{0}_{1}'.format(self.channel,self.systShift) ])
        outdir = "/".join([self.outdir, sample])
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        if not self.submit == "local":
            if self.submit == "condor":
                templ = "lxplus"
                condor_jobs = {}
            else:
                templ = self.submit
            with open("submit_on_{0}.sh".format(templ) ) as FSO:
                templ = string.Template( FSO.read() )

        if not self.debug and not self.submit == "local":
            if not self.writeSubmitLog( sample, self.systShift ):
                print "Sample seems to be running... Make sure you know what you are doing"
                return 

        sleeping = 0
        for idx,configBall in enumerate( self.makeConfigBalls(sample) ):
            file = configBall["file"]

            if self.debug and idx > 0: break

            rundir = "_".join([runpath, str(idx+1)])

            ##### Create rundir. Overwrite if it already exists... Well, actually delete it and make it new. 
            if not os.path.exists(rundir):
                os.makedirs(rundir)
            else:
                shutil.rmtree(rundir)
                os.makedirs(rundir)


            os.chdir(self.basedir)
            self.prepareRunDir(rundir, configBall)

            # Run local
            if self.submit == "local":
                t = threading.Thread(target=self.runOneFileLocal, args=(rundir,file,) )
                threads.append(t)

                while threading.active_count()>self.nthreads:  #pause until thread slots become available
                    pass
                t.start()
                if not self.event: sleep(5)
            # Run on batch system
            else:
                if sleeping > 10: sleeping = 0
                jobname = "{0}+{1}-{2}_{3}".format( sample, self.channel, configBall["systShift"], file.split("/")[-1] )

                runscript = templ.substitute(samplename = jobname,
                                             rundir = rundir,
                                             outdir = outdir,
                                             cell = getSystem(inverse = True),
                                             sleeping = sleeping*10,
                                             channel = self.channel)
                sleeping += 1

                os.chdir(rundir )
                with open("submit.sh","w") as FSO:
                    FSO.write( runscript )
                os.chmod("submit.sh", 0777)

                if self.submit == "hephybatch":
                    os.system( "sbatch submit.sh" )

                if self.submit == "lxplus":
                    os.system( " bsub -q 1nd -J {0} submit.sh".format( jobname ) )

                # if self.submit == "condor":
                    # condor_jobs[ "_".join([self.channel, self.systShift, str(idx+1)]) ] = "/".join([rundir,"condor.sub"])

        if self.submit == "condor":
            os.chdir(self.basedir)
            self.submitToCondor(runpath)


        if self.submit == "local":
            for x in threads:
                x.join()

            if not self.debug:
                os.chdir( "/".join([ self.basedir, "out", sample ]) )
                os.system('hadd -f -O '+'{0}_all.root rundir_{1}_*/{0}_*root'.format("-".join([self.channel, self.systShift]), self.channel ) )

    def submitToCondor(self, runpath):

        run_file = runpath.replace("rundir_","") + ".sub"
        for df in glob(run_file + "*"):
            os.remove(df)

        with open("condor_template.sub","r") as FSO:
            condor_templ = string.Template(FSO.read())

        with open(run_file,"w") as FSO:
            FSO.write(condor_templ.substitute(rundir=runpath+"*" ))
        useToken("cern")
        os.system("condor_submit {0}".format(run_file))


        
            

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

    def prepareRunDir(self, rundir, configBall):

        self.throwConfigBall(configBall, rundir)
        shutil.copytree("proxy", "/".join([rundir,"proxy"]))
        shutil.copytree("kerberos", "/".join([rundir,"kerberos"]))
                
        if not self.submit == "lxplus" and not self.submit == "condor":
            headerfiles = glob("*.h*")
            Cfiles = glob("*.c*") + glob("*.C")
            addFiles =['convertNanoParallel.py','validateAndCopy.py']

            shutil.copytree("utils", "/".join([rundir,"utils"]))
            for f in headerfiles + Cfiles + addFiles:
                shutil.copyfile("/".join([self.basedir,f]), "/".join([rundir,f]) )
                os.chmod("/".join([rundir,f]), 0777)

    def writeSubmitLog(self, sample, shift):

        log = {}

        logpath = "/".join([ getHeplxPublicFolder(),"submit_log.log" ])
        if os.path.exists(logpath):
            with open(logpath,"r") as FSO:
                log = json.load(FSO)

        if not log.get(sample,False): log[sample] = {}
        if not log[sample].get(self.channel,False): log[sample][self.channel] = {}
        if not log[sample][self.channel].get(shift,False): log[sample][self.channel][shift] =  {"submit_time": None, "status":"" }

        if not log[sample][self.channel][shift]["status"] == "NEW" or self.force:
            log[sample][self.channel][shift]["submit_time"] = int(time())
            log[sample][self.channel][shift]["status"] = "NEW"
            log[sample][self.channel][shift]["site"] = self.submit

            with open(logpath,"w") as FSO:
                json.dump(log, FSO, indent=4)
            return True

        return False

    def makeConfigBalls(self,sample):
        configBalls = []
        samples_avail = glob("samples/*/*/*")
        for sa in samples_avail:
            if sample + ".txt" in sa:
                files = self.getFiles(sa)
                parts = sa.split("/")

        with open("tagMapping.json","r") as FSO:
            puTag = json.load(FSO)

        for file in files:
            configBall = {}
            configBall["file"]        = file

            configBall["sample"]      = parts[2]
            configBall["samplename"]  = sample
            configBall["channel"]     = self.channel
            configBall["systShift"]   = self.systShift
            configBall["svfit"]       = self.svfit
            configBall["recoil"]      = self.recoil
            configBall["system"]      = self.submit
            configBall["nevents"]     = int(self.nevents)
            configBall["check_event"] = int(self.event)
            configBall["isSync"]      = self.sync
            if  parts[1] == "mc":
                configBall["isMC"]        = True
                configBall["certJson"] = ""
                configBall["puTag"]    = puTag[ sample ]["putag"]
                configBall["xsec"]     = puTag[ sample ]["xsec"]
                configBall["genNEvents"]  = puTag[ sample ]["nevents"]

            else:
                configBall["isMC"]        = False
                configBall["certJson"] = self.certJson
                configBall["puTag"]   = "pileup"


            configBalls.append(configBall)

        return configBalls

    def throwConfigBall(self, ball, where):

        with open("/".join([where, "configBall.json"]),"w") as FSO:
            json.dump(ball, FSO, indent=4)


    def getFiles(self,sample):

        with open( "{0}/{1}".format(self.basedir, sample) ) as FSO:
            buf = FSO.read()
        return buf.splitlines()    


def checkProxy():
    proxy_path = glob("/tmp/x509*_u{0}".format( os.getuid() ) )
    if len(proxy_path) == 1 and os.path.exists( proxy_path[0] ):

        p = sp.Popen( shlex.split("voms-proxy-info --timeleft"), stdout=sp.PIPE, stderr=sp.PIPE )
        (out,err) =  p.communicate()


        if err:
            print err
            return False
        if out:
            if int(out) > 0:
                if not os.path.exists("proxy"):
                    os.mkdir("proxy")
                shutil.copyfile(proxy_path[0], "proxy/x509_proxy")
                return True
            else:
                print "Proxy not valid! Get a new one with 'voms-proxy-init --voms cms'"
                return False

    
    print "No proxy found! Get a new one with 'voms-proxy-init --voms cms'" 
    return False
    



if __name__ == '__main__':
  main()

