import subprocess as sp
import shlex
import json
import os 

def main():

    if not os.path.exists("samples"):
        os.makedirs("samples")

    print "Searching on DPM"
    samples = buildPathToSamples()

    for sample in samples:

        if not os.path.exists( "/".join( ["samples"] + sample.split("/")[:-1] ) ):
            os.makedirs( "/".join( ["samples"] + sample.split("/")[:-1] ) )

        if not os.path.exists( "samples/{0}.txt".format(sample) ):
            with open( "samples/{0}.txt".format(sample), "w" ) as FSO:
                FSO.write( getFilelist(samples[sample]) )
                FSO.write("\n")

def buildPathToSamples():
    with open("samples_config.json","r") as FSO:
        config = json.load(FSO)

    sampleFormat = {"mc":"NANOAODSIM","data":"NANOAOD"}
    dpm_path = "/dpm/oeaw.ac.at/home/cms/store/"

    sample_collection= {}

    for frmt in config:
        f = config[frmt]

        for run in f["run"]:
            r = f["run"][run]

            for sample in r["samples"]:
                parts = r["samples"][sample]

                for part in parts:
                    path = "/".join([ dpm_path, frmt, run, part, sampleFormat[frmt] ])
                    out = lsDPM(path)
                    for tag in out:

                        if r["creation"] in tag:
                            sample_collection[ "/".join([frmt, sample, "_".join([part, run, r["creation"] ]) ])  ] = "/".join([ path, tag ])

    return sample_collection

def getFilelist(sample):

    filelist = []
    for folder in lsDPM(sample):
        path = "/".join([ sample, folder ])
        files = lsDPM( path )

        for file in files:
            filelist.append( "/".join([ "root://hephyse.oeaw.ac.at/", path, file ]) )

    return "\n".join(filelist)



def lsDPM(path):
    proc = sp.Popen(shlex.split("dpns-ls {0}".format(path) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()

    return out.splitlines()



if __name__ == '__main__':
    main()

# /dpm/oeaw.ac.at/home/cms/store/mc/RunIIFall17NanoAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/
# /dpm/oeaw.ac.at/home/cms/store/data/Run2017B/SingleMuon/NANOAOD/31Mar2018-v1
