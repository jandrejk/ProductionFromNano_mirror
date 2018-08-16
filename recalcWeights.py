import ROOT as R
import root_pandas as rp # From CMSSW >= 9_4_0
import root_numpy as rn
import argparse
from rootpy.io import root_open

parser = argparse.ArgumentParser()

puweightFile = root_open("utils/puweight/puweights.root")
puweights_histo = puweightFile.Get("pileup")

zptweightFile = root_open("utils/zptweight/zpt_weights_2016_BtoH.root")
zpt_histo = zptweightFile.Get("zptmass_histo")

wsp = root_open("utils/CorrectionWorkspaces/htt_scalefactors_v17_1.root")
w = wsp.Get("w")

def main():

	parser.add_argument('-i', dest='inputfile', help='Sample to run over', type=str, metavar = 'INPUTFILE', default = "")
	parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt'], default = 'mt')
	parser.add_argument('-o', dest='outputfile', help='File the changes get written to', type=str, metavar = 'OUTPUTFILE', default = "output.root")
	parser.add_argument('-w', dest='weights', help='Weight that is to be recalculated', choices = ["all", "antilep_tauscaling", "puweight", "zweight", "trk_sf", "reco_sf", "top_weight"], default = "all" )
	
	global args 
	args = parser.parse_args()
	
	ignore_list = ['addlepton_p4'] # root_pandas can't handle vector<TLorentzVector> 
	tmp = rp.read_root( paths = args.inputfile, ignore=ignore_list) 
	
	if args.weights = "all" or args.weights = "puweight": tmp['puweight']	= tmp.apply( recalcPileupWeight, axis=1 )
	#if args.weights = "all" or args.weights = "zweight":  tmp[''] 			= tmp.apply( recalcZWeight, axis=1 )
	if args.weights = "all" or args.weights = "trk_sf":	  tmp['trk_sf']   	= tmp.apply( recalcTrkSF, axis=1 )
	if args.weights = "all" or args.weights = "reco_sf":  tmp['reco_sf'] 	= tmp.apply( recalcRecoSF, axis=1 )
	if args.weights = "all" or args.weights = "antilep_tauscaling":  
		tmp['eleTauFakeRateWeight'] = tmp.apply( recalcEleTauFakeRateWeight, axis=1 )
		tmp['muTauFakeRateWeight'] 	= tmp.apply( recalcMuTauFakeRateWeight, axis=1 )
		tmp['antilep_tauscaling']   = tmp.apply( recalcAntilepTauscaling, axis=1 ) # Always run eleTauFakeRateWeight and muTauFakeRateWeight first!

	#write new file:
	tmp.to_root(args.outputfile, key="TauCheck")

def recalcEleTauFakeRateWeight(row):
	eleTauFakeRateWeight = 1.0
	if args.channel == "mt":
		if( (row['gen_match_2'] == 1 or row['gen_match_2'] == 3) and row['againstElectronVLooseMVA6_2'] > 0.5):
			if(   rn.abs(row['eta_2']) < 1.460 ): eleTauFakeRateWeight *= 1.09
			elif( rn.abs(row['eta_2']) < 1.558 ): eleTauFakeRateWeight *= 1.19
	elif args.channel == "et":
		if( (row['gen_match_2'] == 1 or row['gen_match_2'] == 3) and row['againstElectronTightMVA6_2'] > 0.5):
			if(   rn.abs(row['eta_2']) < 1.460 ): eleTauFakeRateWeight *= 1.80
			elif( rn.abs(row['eta_2']) < 1.558 ): eleTauFakeRateWeight *= 1.53
	elif args.channel == "tt":
		if( (row['gen_match_1'] == 1 or row['gen_match_1'] == 3) and row['againstElectronVLooseMVA6_1'] > 0.5):
			if(   rn.abs(row['eta_1']) < 1.460 ): eleTauFakeRateWeight *= 1.09
			elif( rn.abs(row['eta_1']) < 1.558 ): eleTauFakeRateWeight *= 1.19
		if( (row['gen_match_2'] == 1 or row['gen_match_2'] == 3) and row['againstElectronVLooseMVA6_2'] > 0.5):
			if(   rn.abs(row['eta_2']) < 1.460 ): eleTauFakeRateWeight *= 1.09
			elif( rn.abs(row['eta_2']) < 1.558 ): eleTauFakeRateWeight *= 1.19
	return eleTauFakeRateWeight
	
def recalcMuTauFakeRateWeight(row):
	muTauFakeRateWeight = 1.0
	if args.channel == "mt":
		if( (row['gen_match_2'] == 2 or row['gen_match_2'] == 4) and row['againstMuonTight3_2'] > 0.5):
			if  ( rn.abs(row['eta_2']) < 0.4 ): muTauFakeRateWeight *= 1.17
			elif( rn.abs(row['eta_2']) < 0.8 ): muTauFakeRateWeight *= 1.29
			elif( rn.abs(row['eta_2']) < 1.2 ): muTauFakeRateWeight *= 1.14
			elif( rn.abs(row['eta_2']) < 1.7 ): muTauFakeRateWeight *= 0.93
			elif( rn.abs(row['eta_2']) < 2.3 ): muTauFakeRateWeight *= 1.61
	elif args.channel == "et":
		if( (row['gen_match_2'] == 2 or row['gen_match_2'] == 4) and row['againstMuonLoose3_2'] > 0.5):
			if  ( rn.abs(row['eta_2']) < 0.4 ): muTauFakeRateWeight *= 1.06
			elif( rn.abs(row['eta_2']) < 0.8 ): muTauFakeRateWeight *= 1.02
			elif( rn.abs(row['eta_2']) < 1.2 ): muTauFakeRateWeight *= 1.10
			elif( rn.abs(row['eta_2']) < 1.7 ): muTauFakeRateWeight *= 1.03
			elif( rn.abs(row['eta_2']) < 2.3 ): muTauFakeRateWeight *= 1.94
	elif args.channel == "tt":
		if( (row['gen_match_1'] == 2 or row['gen_match_1'] == 4) and row['againstMuonLoose3_1'] > 0.5):
			if  ( rn.abs(row['eta_1']) < 0.4 ): muTauFakeRateWeight *= 1.06
			elif( rn.abs(row['eta_1']) < 0.8 ): muTauFakeRateWeight *= 1.02
			elif( rn.abs(row['eta_1']) < 1.2 ): muTauFakeRateWeight *= 1.10
			elif( rn.abs(row['eta_1']) < 1.7 ): muTauFakeRateWeight *= 1.03
			elif( rn.abs(row['eta_1']) < 2.3 ): muTauFakeRateWeight *= 1.94
		if( (row['gen_match_2'] == 2 or row['gen_match_2'] == 4) and row['againstMuonLoose3_2'] > 0.5):
			if  ( rn.abs(row['eta_2']) < 0.4 ): muTauFakeRateWeight *= 1.06
			elif( rn.abs(row['eta_2']) < 0.8 ): muTauFakeRateWeight *= 1.02
			elif( rn.abs(row['eta_2']) < 1.2 ): muTauFakeRateWeight *= 1.10
			elif( rn.abs(row['eta_2']) < 1.7 ): muTauFakeRateWeight *= 1.03
			elif( rn.abs(row['eta_2']) < 2.3 ): muTauFakeRateWeight *= 1.94
	return muTauFakeRateWeight
	
def recalcAntilepTauscaling(row):
	return row['eleTauFakeRateWeight'] * row['muTauFakeRateWeight']

def recalcPileupWeight(row):
	npu = row['npu']
	return puweights_histo.GetBinContent( puweights_histo.GetXaxis().FindBin(npu) );

def recalcZWeight(row): 
	gen_pt = rn.sqrt( rn.power(row['gen_ll_px'],2) + rn.power(row['gen_ll_py'],2) )
	return zpt_histo.GetBinContent(zpt_histo.GetXaxis().FindBin(row['gen_Mll']),zpt_histo.GetYaxis().FindBin(gen_pt))

def recalcTopWeight(row): 
	top_pt = []
	for i in range(row['ngenSum']):
		if ( rn.abs( row['genSum_pdgId'][i]) == 6 and row['genSum_fromHardProcess'][i] and row['genSum_isLastCopy'][i] ) :
			top_pt.append(row['genSum_pt'])
	if len(top_pt) != 2:
		return 1.
	return rn.sqrt( rn.exp(0.0615-0.0005*top_pt[0]) * rn.exp(0.0615-0.0005*top_pt[1]) )

def recalcTrkSF(row):
	return w.function("m_trk_ratio").getVal()
	
def recalcRecoSF(row):
	return w.function("e_reco_ratio").getVal()
	
if __name__ == '__main__':
	main()

# https://github.com/mflechl/HephyHiggs/blob/ML_SM2016/production2016_root6/src/syncDATA.cc#L2626-L2672
# https://github.com/mflechl/HephyHiggs/blob/ML_SM2016/production2016_root6/src/TNtupleAnalyzer.cc#L176-L180
# https://github.com/mflechl/HephyHiggs/blob/ML_SM2016/production2016_root6/src/TNtupleAnalyzer.cc#L191-L195
# https://github.com/mflechl/HephyHiggs/blob/ML_SM2016/production2016_root6/src/syncDATA.cc#L1757-L1857
