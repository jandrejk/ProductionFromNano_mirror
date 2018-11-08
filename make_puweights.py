import ROOT as R
import copy
import os



def getName (word) :
    partition = word.split('_')
    new_word = ['#']
    for i, p in enumerate(partition) :
	if 'RunIIFall17NanoAOD' in p or 'NANOAODSIM' in p :
	    new_word.append('#')
	else :
            if i > 1 :
 	    	new_word.append('_')
	new_word.append(p)

    #print new_word
    return ''.join(new_word)[:-5] # [:-5] to cut off the .root at the end of the name


input_file = R.TFile("./MCpuHistograms/merged/pileup_distribution_data2017.root")

data_hist = input_file.Get("pileup")
data_hist.Scale( 1/ data_hist.Integral() )

histos = []
name_list = os.listdir('/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/MCpuHistograms/merged/')
names = []
#print name_list
#chain = R.TChain("collection")


for i in name_list :
	if "NANOAODSIM" in i:	
		ifile = R.TFile('./MCpuHistograms/merged/'+i)
		tmp = ifile.Get("h1")
		tmp.Scale( 1/tmp.Integral() )
		histos.append(copy.deepcopy(tmp))
		ifile.Close()
		names.append(i)

histos.append(copy.deepcopy(data_hist))
names.append("pileup_distribution_data2017.root")

output_file = R.TFile("./MCpuHistograms/merged/puweights.root","RECREATE")
output_file.cd()

for j, h in enumerate(histos):
	tmp_mc = copy.deepcopy(h)
	tmp_data = copy.deepcopy(data_hist)
	hist_setname = getName(word = names[j])
	print hist_setname
	tmp_data.SetName( hist_setname )
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



