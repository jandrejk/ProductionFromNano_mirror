import subprocess as sp
import shlex
import json
import os 


def main():
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

                    if "Vienna" in getDASquery(link,"site"):
                        content = content.replace("/store/","root://hephyse.oeaw.ac.at//store/")
                    else:
                        content = content.replace("/store/","root://xrootd-cms.infn.it//store/")

                    writeFile(content, filename, folder)
                    # print getDASquery(link,"file")


def getDASquery(link, info="file"):

    proc = sp.Popen(shlex.split('dasgoclient --query="{0} dataset={1}" '.format(info,link) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()

    return out

def writeFile( content, filename, folder ):
    if not os.path.exists(folder):
        os.mkdirs(folder)

    with open("/".join([folder,filename]), "w" ) as FSO:
        FSO.write(content)





if __name__ == '__main__':
    main()

# /dpm/oeaw.ac.at/home/cms/store/mc/RunIIFall17NanoAOD/VBFHToTauTau_M125_13TeV_powheg_pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/
# /dpm/oeaw.ac.at/home/cms/store/data/Run2017B/SingleMuon/NANOAOD/31Mar2018-v1
