import ROOT as R
import copy
import os
from glob import glob



def makePuweights():

    input_file = R.TFile("pileup_distribution_data2017.root")

    data_hist = input_file.Get("pileup")
    data_hist.Scale( 1/ data_hist.Integral() )

    histos = {}



    for histo in  glob('{CMSSW_BASE}/src/WawTools/NanoAODTools/MCpuHistograms/merged/*'.format(**os.environ) ) :

        ifile = R.TFile(histo)
        tmp = ifile.Get("h1")
        tmp.Scale( 1/tmp.Integral() )
        histos[ histo.replace(".root","").split("/")[-1] ] =  copy.deepcopy(tmp)
        ifile.Close()
        # names.append(i)

    histos["data_pileup"] = copy.deepcopy(data_hist)

    if not os.path.exists("./MCpuHistogram/result") :
        os.makedirs("./MCpuHistograms/result")

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
    makePuweights()
