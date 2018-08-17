import ROOT as R
import json
import os
import glob
import shutil
import subprocess as sp
import shlex
from runUtils import checkProxy, checkTokens, getSystem, getHeplxPublicFolder

def main():

    if not checkTokens(): sys.exit()

    M = Merger()
    M.mergeSamples()

class Merger():

    def __init__(self, force=False, version="v1"):

        self.force = force
        self.version = version

        if "cern" in getSystem():
            print "Merging from lxplus - You should switch to heplx to run faster."
  

        self.logpath = "/".join([ getHeplxPublicFolder(),"submit_log.log" ])
        with open(self.logpath,"r") as FSO:
            self.log = json.load(FSO)

        self.outdir = "/afs/hephy.at/data/higgs01"

        self.samples = self.collectFiles()
        self.mergekeys = self.mapCompletedJobs()

    def __del__(self):

        with open(self.logpath,"w") as FSO:
            json.dump(self.log, FSO, indent=2)

    def mergeSamples(self):
        for m in self.mergekeys:
            mergedir = "/".join([self.outdir,m[0],self.version])
            if os.path.exists(mergedir):
                shutil.rmtree(mergedir, ignore_errors=True)
            os.makedirs(mergedir)

            outfile = "/".join([mergedir, "{0}-{1}.root".format(m[1], m[2])])
            print "\033[1m{0}\033[0m".format(m[0])
            # Fallback
            # os.system("hadd -f {0} {1}".format(outfile, " ".join(self.samples[m[0]][m[1]][m[2]]["files"]) ) ) 

            self.mergeSample(outfile, mergekey=m)

    def mergeSample(self, outfile, mergekey):
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

    def stitchSamples(self):

        with open("stitchConfig.json","r") as FSO:
            stitch_config = json.load(FSO) 

    def mapCompletedJobs(self):

        mergekeys = []
        samples = self.log.keys()
        samples.sort()
        for sample in samples:
            for channel in self.log[sample]:
                for shift in self.log[sample][channel]:
                    if self.log[sample][channel][shift]["status"] == "MERGE" or (self.force and self.log[sample][channel][shift]["status"] == "DONE"):
                        mergekeys.append((sample,channel,shift))
        return mergekeys

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
