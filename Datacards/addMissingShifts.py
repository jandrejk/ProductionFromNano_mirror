import root_pandas as rp
import sys
import numpy as np
from ROOT import TFile, gDirectory, TH1F

def GetKeyNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetName() for key in gDirectory.GetListOfKeys()]

channel = 'tt'

def MissingShift(name) :
	missing_shifts = [
		'CMS_htt_boson_reso_met_13TeV',
		'CMS_htt_boson_scale_met_13TeV',
		'CMS_htt_eff_b_13TeV',
		'CMS_htt_mistag_b_13TeV',
		'CMS_scale_met_unclustered_13TeV']
	for word in missing_shifts :
		if word in name :
			return True

	return False


def WriteHisto (directory, Histkey) :
	out_histo = directory.Get(Histkey.GetName())
	out_histo.SetDirectory(0)
	out_histo.Write()

path_Vienna = "/afs/hephy.at/user/m/mspanring/public/forJanik/datacards/vienna/htt_{0}.inputs-sm-13TeV-ML.root".format(channel)
#df_Vienna = rp.read_root(path_Vienna, key= "TauCheck", ignore=["addlepton*"] )
#df.to_root("{0}_{1}.root".format(name,channel), key="TauCheck" )

path_KIT = "/afs/hephy.at/user/m/mspanring/public/forJanik/datacards/kit/htt_{0}.inputs-sm-13TeV-ML.root".format(channel)
#df_KIT = rp.read_root(path_KIT)

#print df_KIT.columns.T
Vienna_file = TFile.Open(path_Vienna,'read')
Tdirectories_Vienna = [key for key in Vienna_file.GetListOfKeys()]


KIT_file = TFile.Open(path_KIT,'read')
Tdirectory_names_KIT = [key.GetName() for key in KIT_file.GetListOfKeys()]


output = TFile.Open("/afs/hephy.at/user/j/jandrejkovic/public/forMarkus/datacards/vienna/htt_{0}.inputs-sm-13TeV-ML_copyMissingDatacards.root".format(channel),'recreate')

#print Tdirectories_KIT
#print Tdirectories_Vienna

counter = 0
loop = 0
#print Tdirectories_Vienna 
for Tdir in Tdirectories_Vienna :
	
	output.mkdir(Tdir.GetName())
	output.cd(Tdir.GetName())
	#print Tdir

	d = Vienna_file.Get(Tdir.GetName())
	#print d
	#print d.GetListOfKeys()
	histo_name_keys = [key for key in d.GetListOfKeys()]


	# is Tdirectory name from Vienna also in KIT ?
	if Tdir.GetName() in Tdirectory_names_KIT :
		#print 'copying missing datacards'
		for h_key in histo_name_keys :

			# 
			if MissingShift(name=h_key.GetName()) :

				d_KIT = KIT_file.Get(Tdir.GetName())
				KIT_histo_names = [key.GetName() for key in d_KIT.GetListOfKeys()]

				# check if histogram exists in KIT file
				if (h_key.GetName() in KIT_histo_names) :
					#print Tdir.GetName()
					#print h_key.GetName()
					#print 'copy from KIT'
					loop += 1
					WriteHisto(directory=d_KIT, Histkey=h_key)
				else :
					#print 'copy from Vienna'
					WriteHisto(directory=d, Histkey=h_key)
					
				
			else :
				#print 'copy from Vienna'
				WriteHisto(directory=d, Histkey=h_key)
				#out_histo = d.Get(h_key.GetName())
				#out_histo.SetDirectory(0)
				#ut_histo.Write()

			counter +=1
	# copy all histograms from Vienna
	else :
		for h_key in histo_name_keys[:2] :
			WriteHisto(directory=d, Histkey=h_key)
			counter +=1
		
print 'No. of Datacards: ', counter
print "No. of Datacards taken from KTI: ", loop
print "---END--- addMissingShifts.py"