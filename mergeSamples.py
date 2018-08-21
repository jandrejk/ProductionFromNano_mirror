import ROOT as R
import json
import os
import glob
import shutil
import subprocess as sp
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
    M.collectSamples()
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

    def collectSamples(self):

        for m in self.mergekeys:
            mergedir = "/".join([self.outdir,m[0],self.version])
            if not os.path.exists(mergedir):
                os.makedirs(mergedir)

            outfile = "/".join([mergedir, "{1}-{2}_{0}.root".format( *m ) ])
            if os.path.exists(outfile):
                shutil.rmtree(outfile, ignore_errors=True)

            print "\033[1m{0}\033[0m".format(m[0])
            # Fallback
            # os.system("hadd -f {0} {1}".format(outfile, " ".join(self.samples[m[0]][m[1]][m[2]]["files"]) ) ) 

            self.collectSample(outfile, mergekey=m)

    def collectSample(self, outfile, mergekey):
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

        with open("tagMapping.json","r") as FSO:
            lumi_infos = json.load(FSO)

        mergedir = "/".join([self.outdir,self.version])
        if not os.path.exists(mergedir):
            os.makedirs(mergedir)

        for mergename in stitch_config:
            # if mergename == "BASIS_ntuple_VBF":
                to_merge = {}
                for sample in stitch_config[mergename]["samples"]:

                    to_merge[sample] = {"files":[],"xsec":0,"nevents":[] }
                    for ext in glob.glob( "/".join([ self.outdir, sample + "*" ]) ):
                        tag = ext.replace( self.outdir, "" )
                        to_merge[sample]["files"].append(ext)

                        if stitch_config[mergename]["data"]:
                            to_merge[sample]["xsec"] = 1
                            to_merge[sample]["nevents"] = [1]
                        else:
                            to_merge[sample]["xsec"] = lumi_infos[sample][1]
                            to_merge[sample]["nevents"].append(lumi_infos[tag][2])                        

                    to_merge[sample]["nevents"] = sum(to_merge[sample]["nevents"]) 
                    if to_merge[sample]["nevents"] == 0:
                        print "Samples for {0} missing...".format( mergename )
                        to_merge = {}
                        break

                if to_merge:
                    if stitch_config[mergename]["stitching"]:
                        self.mapJetsStitchingWeights( to_merge, stitch_config[mergename]["NNLO_xsec"] )
                    else:
                        self.mapLumiWeight(to_merge)

                    self.mergeSample(mergename, to_merge)


    def mergeSample(self, name, parts):

        mergedir = "/".join([self.outdir,self.version])
        shifts = ["NOMINAL"]
        # for es in ["TES","MES","EES"]:        
        #     for dm in ["1p0p0","1p1p0","3p0p0"]:
        #         for sh in ["Up","Down"]:
        #             shifts.append( es + dm + sh )

        for shift in shifts:


            filename = "_".join([ name.replace("BASIS",shift), self.channel ]) + ".root"
            outfile =  "/".join([mergedir, filename ])
            addfiles = []
            print outfile        
            complete = True

            for i,part in enumerate(parts):
                for j,file in enumerate(parts[part]["files"]):
                    samples = glob.glob( "/".join([ file, self.version, "{0}-{1}*".format(self.channel, shift) ]) )
                    if not samples:
                        complete = False
                    for sample in samples:
                        addfiles.append(sample)
            if complete:
                os.system("hadd -f {0} {1}".format(outfile, " ".join( addfiles ) ) ) 



    def mapLumiWeight(self, parts):
        for part in parts:
            parts[part]["lumiWeight"] = parts[part]["xsec"] / parts[part]["nevents"]

    def mapJetsStitchingWeights(self, parts, NNLO_xsec):

        lumis = [0,0,0,0,0]
        corr = NNLO_xsec
        for part in parts:
            use = -1
            for i,jetM in enumerate(["","1","2","3","4"]):
                if jetM + "Jets" in part:
                    use = i
            if use >= 0:
                lumis[use] = parts[part]["nevents"] / parts[part]["xsec"]
            if use == 0:
                corr /= parts[part]["xsec"]

        stitchWeights = []
        for i in xrange(5):
            if i == 0:
                stitchWeights.append( corr / lumis[i] )
            else:
                stitchWeights.append( corr / ( lumis[0] + lumis[i] ) )

        for part in parts:
            if "Jets" in part: parts[part]["lumiWeight"] = stitchWeights
            else: parts[part]["lumiWeight"] = parts[part]["xsec"] / parts[part]["nevents"]
        


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
