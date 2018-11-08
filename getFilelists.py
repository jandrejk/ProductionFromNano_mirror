import subprocess as sp
import shlex
import json
import shutil
import os 
from runUtils import checkProxy
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', dest='location', help='Get samples from local or DAS', choices=["das","private"], default = "das")
    args = parser.parse_args()

    if args.location == "private": 
        getfromPrivate("/dpm/oeaw.ac.at/home/cms/store/user/jaandrej/")
    else:
        getfromDAS()

def getfromPrivate(location):

    if not os.path.exists( "samples" ):
        os.mkdir( "samples" )

    for tier in ["samples/data","samples/mc"]:
        if not os.path.exists( tier ):
            os.mkdir(tier)
        folders = getPrivateFiles(location + tier )
        k = folders.keys()
        k.sort()
        for folder in k:
            if folders[folder]:
                sub = folder.replace(location+tier,"").split("/")[1]
                if not os.path.exists( "/".join([tier, sub]) ):
                    os.mkdir( "/".join([ tier, sub]) )

                with open( "/".join([tier,sub, folder.replace(location+tier,"").split("/")[2]]) + ".txt","w") as FSO:
                    FSO.write( "\n".join( folders[folder] ) )

def getPrivateFiles(source):

    proc = sp.Popen(shlex.split("dpns-ls -R {0}".format(source) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()
    folder = {}
    current = ""
    for line in out.splitlines():
        if ":" in line:
            folder[line.replace(":","")] = []
            current = line.replace(":","")
        if current and ".root" in line:
            folder[current].append("/".join(["root://hephyse.oeaw.ac.at",current, line]) )

    return folder

def getfromDAS():

    if not checkProxy(): sys.exit()

    if os.path.exists("samples"):
        shutil.rmtree("samples")

    with open("sample_collection.json","r") as FSO:
        config = json.load(FSO)

    for frmt in config:
        f = config[frmt]

        for run in f["run"]:
            r = f["run"][run]

            for sample in r["samples"]:
                links = r["samples"][sample]

                for link in links:
                    folder = "/".join(["samples",frmt, sample])
                    filename = buildFileName( link, run, r["creation"] )
                    content = getDASquery(link)

                    content = content.replace("/store/","root://cms-xrd-global.cern.ch//store/")


                    writeFile(content, filename, folder)
                    # print getDASquery(link,"file")


def getDASquery(link, info="file"):

    proc = sp.Popen(shlex.split('dasgoclient --query="{0} dataset={1}" '.format(info,link) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()

    return out

def writeFile( content, filename, folder ):
    if not os.path.exists(folder):
        os.makedirs(folder)

    with open("/".join([folder,filename]), "w" ) as FSO:
        FSO.write(content)

def buildFileName(link, run, creation):

    parts = link.split("/")

    ext = ""
    if "_ext" in parts[2]:
        ext = "_ext" + parts[2].split("_ext")[1].split("-")[0]

    return "_".join([parts[1],run, creation]) + ext + ".txt"



if __name__ == '__main__':
    main()

# /dpm/oeaw.ac.at/home/cms/store/mc/RunIIFall17NanoAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/
# /dpm/oeaw.ac.at/home/cms/store/data/Run2017B/SingleMuon/NANOAOD/31Mar2018-v1
