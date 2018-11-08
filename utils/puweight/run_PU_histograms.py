import sys, os
from glob import glob
import shlex
import argparse
import string
import json
import copy
import ROOT as R

# if the output directory does not exist, create it

def getSamples():
    samples = set()
    for file in glob("{CMSSW_BASE}/src/WawTools/NanoAODTools/samples/mc/*/*".format(**os.environ) ):
        samples.add( file.replace(".txt","").split("/")[-1] )
    return samples

def mergeFragments():
    if not os.path.exists("./MCpuHistograms/merged") :
        os.makedirs("./MCpuHistograms/merged")
    if not os.path.exists("./MCpuHistograms/usedHistos") :
        os.makedirs("./MCpuHistograms/usedHistos")

    for sample in getSamples() :
        #print key  
        save_name = './MCpuHistograms/merged/'+sample+'.root'
        if glob( 'MCpuHistograms/*'+ sample +'*'):
            cmd = 'hadd '+save_name+' MCpuHistograms/*'+ sample +'*'
            mvcmd = 'mv ./MCpuHistograms/*'+sample+'*root ./MCpuHistograms/usedHistos/'
            os.system(cmd)
            os.system(mvcmd)

def  createPUHistos():
    if not os.path.exists("./MCpuHistograms") :
        os.makedirs("./MCpuHistograms")
    if not os.path.exists("./rundir") :
        os.makedirs("./rundir")

    base_dir = "{CMSSW_BASE}/src/WawTools/NanoAODTools/".format(**os.environ)

    file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]

    job_counter = 0
    for file_name in file_list :
        for sample in open(file_name,'r') :
            os.chdir( base_dir + "/utils/puweight")
            with open("template.sh".format('batch') ) as FSO:
                templ = string.Template( FSO.read() )
            runscript = templ.substitute(sample_name = sample)
            os.chdir('./rundir')
            with open("submit_"+str(job_counter)+".sh","w") as FSO:
                FSO.write( runscript )
                
            os.chmod("submit_"+str(job_counter)+".sh", 0777)

            #print("sbatch submit_"+str(job_counter)+".sh") 
            os.system( "sbatch submit_"+str(job_counter)+".sh" )
            job_counter += 1

def makePuweights():
    os.chdir("{CMSSW_BASE}/src/WawTools/NanoAODTools/utils/puweight/".format(**os.environ))
    input_file = R.TFile("pileup_distribution_data2017.root")

    data_hist = input_file.Get("pileup")
    data_hist.Scale( 1/ data_hist.Integral() )

    histos = {}

    for histo in  glob('{CMSSW_BASE}/src/WawTools/NanoAODTools/utils/puweight/MCpuHistograms/merged/*'.format(**os.environ) ) :

        ifile = R.TFile(histo)
        tmp = ifile.Get("h1")
        tmp.Scale( 1/tmp.Integral() )
        histos[ histo.replace(".root","").split("/")[-1] ] =  copy.deepcopy(tmp)
        ifile.Close()
        # names.append(i)

    histos["pileup"] = copy.deepcopy(data_hist)

    if not os.path.exists("MCpuHistograms/result") :
        os.makedirs("MCpuHistograms/result")

    output_file = R.TFile("./MCpuHistograms/result/puweights.root","RECREATE")
    output_file.cd()

    for name, histo in histos.items():
        tmp_mc = copy.deepcopy(histo)
        tmp_data = copy.deepcopy(data_hist)

        print name
        tmp_data.SetName( name )
        maxbins = min( tmp_data.GetNbinsX(), tmp_mc.GetNbinsX() )
        #break
        for i in xrange(maxbins):

            db = tmp_data.GetBinContent(i)
            dm = tmp_mc.GetBinContent(i)

            if dm > 0:
                tmp_data.SetBinContent(i,db/dm )
            else:
                tmp_data.SetBinContent(i,0.)

        tmp_data.Write()

    output_file.Close()
    input_file.Close()

    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', dest='step', choices=["1","2"], default = "" )
    args = parser.parse_args()

    if args.step == "1":
        createPUHistos()
    if args.step == "2":
        mergeFragments()
        makePuweights()

