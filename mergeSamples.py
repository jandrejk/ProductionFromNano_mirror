import ROOT as R
import json
import os
import glob
import shutil
import subprocess as sp
import shlex

def main():
    M = Merger()
    M.mergeSamples()

class Merger():

    def __init__(self, force=False, version="v1"):

        self.force = force
        self.version = version

        host = os.environ["HOSTNAME"]
        if "heplx" in host: 
          system = "hephybatch"
        elif "lxplus" in host: 
          system = "lxbatch"

        with open("submit_log.log","r") as FSO:
            log = json.load(FSO)

        self.log = log[system]
        

        self.outdir = self.log.pop("outdir",None) + "/"
        self.outdir = self.outdir.replace("//","/")

        self.samples = self.collectFiles()
        self.mergekeys = self.mapCompletedJobs()
        

    def mergeSamples(self):
        for m in self.mergekeys:
            mergedir = "/".join([self.outdir,m[0],self.version])
            if os.path.exists(mergedir):
                shutil.rmtree(mergedir, ignore_errors=True)
            os.mkdir(mergedir)

            outfile = "/".join([mergedir, "{0}-{1}.root".format(m[1], m[2])])
            os.system("hadd {0} {1}".format(outfile, " ".join(self.samples[m[0]][m[1]][m[2]]["files"]) ))

    def mapCompletedJobs(self):

        mergekeys = []

        for sample in self.log:
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
