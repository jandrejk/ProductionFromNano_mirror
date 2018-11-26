from ROOT import TFile, gDirectory, TH1F
import argparse

def GetKeyNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetName() for key in gDirectory.GetListOfKeys()]

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

def main () :
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt'], default = 'mt')
	args = parser.parse_args()
	channel = args.channel
	print "=== Decay channel: {0} ===".format(channel)

	path_Vienna = "/afs/hephy.at/user/m/mspanring/public/forJanik/datacards/vienna/htt_{0}.inputs-sm-13TeV-ML.root".format(channel)
	path_KIT = "/afs/hephy.at/user/m/mspanring/public/forJanik/datacards/kit/htt_{0}.inputs-sm-13TeV-ML.root".format(channel)

	Vienna_file = TFile.Open(path_Vienna,'read')
	Tdirectories_Vienna = [key for key in Vienna_file.GetListOfKeys()]

	KIT_file = TFile.Open(path_KIT,'read')
	Tdirectory_names_KIT = [key.GetName() for key in KIT_file.GetListOfKeys()]


	output = TFile.Open("/afs/hephy.at/user/j/jandrejkovic/public/forMarkus/datacards/vienna/htt_{0}.inputs-sm-13TeV-ML_copyMissingDatacards.root".format(channel),'recreate')

	counter = 0
	loop = 0
	# Loop over all TDirectories in Vienna file
	for Tdir in Tdirectories_Vienna :
		
		output.mkdir(Tdir.GetName())
		output.cd(Tdir.GetName())
		
		d = Vienna_file.Get(Tdir.GetName())
		histo_name_keys = [key for key in d.GetListOfKeys()]


		# is Tdirectory from Vienna also in KIT ?
		if Tdir.GetName() in Tdirectory_names_KIT :
			
			# loop over all Datacards in TDirectory
			for h_key in histo_name_keys :

				# Check if Datacard name belongs to the list of missing names
				if MissingShift(name=h_key.GetName()) :

					d_KIT = KIT_file.Get(Tdir.GetName())
					KIT_histo_names = [key.GetName() for key in d_KIT.GetListOfKeys()]

					# check if Datacard exists in KIT file
					if (h_key.GetName() in KIT_histo_names) : # if yes copy Datacard from KIT
						loop += 1
						WriteHisto(directory=d_KIT, Histkey=h_key)
					else : # copy Datacard from Vienna
						WriteHisto(directory=d, Histkey=h_key)
						
					
				else : # Datacard is not missing, therefore take it from Vienna
					WriteHisto(directory=d, Histkey=h_key)
				
				counter +=1
	
		else : # Tdirectory only in Vienna, therefore copy all Datacards from Vienna that are inside this TDirectory
			for h_key in histo_name_keys[:2] :
				WriteHisto(directory=d, Histkey=h_key)
				counter +=1

	print 'No. of Datacards: ', counter
	print "No. of Datacards taken from KIT: ", loop
if __name__ == '__main__':
  	main()		
	print "---END--- addMissingShifts.py"