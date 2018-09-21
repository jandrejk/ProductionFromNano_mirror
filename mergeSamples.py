import ROOT as R
import json
import os
import glob
import shutil
import subprocess as sp
import multiprocessing as mp
import shlex
import sys
import numpy as np
import argparse
from runUtils import checkProxy, checkTokens, getSystem, getHeplxPublicFolder

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', dest='version', help='Version of the merged samples', default = "v1")
    parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt'], default = 'mt')
    parser.add_argument('-s', dest='stitch', help='Stitch samples',action="store_true")


    args = parser.parse_args()

    if not checkTokens(): sys.exit()

    M = Merger(version=args.version, channel = args.channel)
    M.createSamples()
    if args.stitch: M.mergeSamples()

class Merger():

    def __init__(self, force=False, version="v1", channel = "mt"):

        self.force = force
        self.version = version
        self.channel = channel

        if "cern" in getSystem():
            print "Merging from lxplus - You should switch to heplx to run faster."
  

        self.logpath = "/".join([ getHeplxPublicFolder(),"submit_log.log" ])
        with open(self.logpath,"r") as FSO:
            self.log = json.load(FSO)

        self.outdir = "/afs/hephy.at/data/higgs01/"

        self.samples = self.collectFiles()
        self.mergekeys = self.mapCompletedJobs()

    def __del__(self):

        with open(self.logpath,"w") as FSO:
            json.dump(self.log, FSO, indent=2)

    def createSamples(self):

        for m in self.mergekeys:
            mergedir = "/".join([self.outdir,m[0],self.version])
            if not os.path.exists(mergedir):
                os.makedirs(mergedir)

            outfile = "/".join([mergedir, "{1}-{2}_{0}.root".format( *m ) ])
            if os.path.exists(outfile):
                shutil.rmtree(outfile, ignore_errors=True)

            print "\033[1m{0}\033[0m: {1}-{2}".format(*m)
            # Fallback
            # os.system("hadd -f {0} {1}".format(outfile, " ".join(self.samples[m[0]][m[1]][m[2]]["files"]) ) ) 

            self.combineFiles(outfile, mergekey=m)

    def combineFiles(self, outfile, mergekey):
            m = mergekey
            FM = R.TFileMerger()
            FM.OutputFile(outfile)
            success = True
            print "Adding files for merging... "
            for f in self.samples[m[0]][m[1]][m[2]]["files"]:
                if not FM.AddFile(f,False):
                    success = False
                    break
            if success:
                if not FM.Merge():
                    print "\033[0;31mProblem during merge...\033[0m\n"
                else:
                    print "\033[0;32mMerge successful!\033[0m\n"
                    self.log[m[0]][m[1]][m[2]]["status"] = "DONE"

    def mergeSamples(self):

        with open("stitchConfig.json","r") as FSO:
            stitch_config = json.load(FSO)

        cmd_list = {}

        mergedir = "/".join([self.outdir,self.version])
        if not os.path.exists(mergedir):
            os.makedirs(mergedir)

        for mergename in stitch_config:
            # if mergename == "BASIS_ntuple_VBF":
                tomerge = []
                for sample in stitch_config[mergename]["samples"]:
                    for ext in glob.glob( "/".join([ self.outdir, sample + "*", self.version, self.channel + "-*" ]) ):
                        tomerge.append(ext)

                mergecmds = self.getMergeCmds(mergename, tomerge) 
                if mergecmds:
                    for mergecmd in mergecmds:
                        os.system(mergecmd)
        #             cmd_list[mergename] = mergecmd 

        # self.applyCmdMulti(cmd_list)
  
    def getMergeCmds(self, name, parts):

        mergedir = "/".join([self.outdir,self.version])
        shifts = ["NOMINAL"]
        for es in ["TES","MES","EES"]:        
            for dm in ["1p0p0","1p1p0","3p0p0"]:
                for sh in ["Up","Down"]:
                    shifts.append( es + dm + sh )

        mergeCmds = []

        for shift in shifts:

            filename = "-".join([ self.channel, name.replace("BASIS",shift) ]) + ".root"
            if "_SingleElectron" in filename: filename = filename.replace("SingleElectron","Data")
            if "_SingleMuon" in filename: filename = filename.replace("SingleMuon","Data")
            if "_Tau" in filename: filename = filename.replace("Tau","Data")
            outfile =  "/".join([mergedir, filename ])
            addfiles = []

            for i,file in enumerate(parts):
                if "{0}-{1}".format( self.channel, shift ) in file:
                    addfiles.append(file)
            if addfiles:
                mergeCmds.append("hadd -f {0} {1}".format(outfile, " ".join( addfiles ) ) )
        return mergeCmds



    def applyCmdMulti(self, cmd_list, max_proc=8):

        done_queue = mp.Queue()

        for i, mergename in enumerate(cmd_list):
            print "Stitchhing: " + mergename
            if i >= max_proc:
                done_queue.get(block=True)
            proc = mp.Process(target=self.exec_cmd, args=(cmd_list[mergename], done_queue ))
            proc.start()

    def exec_cmd(self, cmd, q):
        shlCmd = shlex.split(cmd)
        # print shlCmd
        p = sp.Popen(shlCmd,stdout = sp.PIPE, stderr = sys.__stderr__, shell=False)
        p.communicate()
        q.put(object)

    def mapCompletedJobs(self):

        mergekeys = []
        samples = self.log.keys()
        samples.sort()
        for sample in samples:
            for channel in self.log[sample]:
                for shift in self.log[sample][channel]:
                    if self.log[sample][channel][shift]["status"] == "MERGE" or (self.force and self.log[sample][channel][shift]["status"] == "DONE"):
                        mergekeys.append((sample,channel,shift))
        return np.array(mergekeys)

    def collectFiles(self):

        samples = {}
        for sample in glob.glob(self.outdir + "*"):
            samplename = sample.replace(self.outdir,"")
            samples[samplename] = {}

            for file in glob.glob( "/".join([ sample, "*"]) ):
                root_file = file.split("/")[-1]
                channel = root_file.split("-")[0]
                shift = root_file.split("_")[0].replace(channel+"-","")

                if not samples[samplename].get(channel,False): samples[samplename][channel] = {}
                if not samples[samplename][channel].get(shift,False): samples[samplename][channel][shift] = {"files":[]}

                samples[samplename][channel][shift]["files"].append( "/".join([sample,root_file]))

        return samples

if __name__ == '__main__':
    main()
