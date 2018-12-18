from ROOT import TFile, gDirectory, TH1F
import argparse

def GetKeyNames( self, dir = "" ):
        self.cd(dir)
        return [key.GetName() for key in gDirectory.GetListOfKeys()]

def MissingShift(name) :
	missing_shifts = [
		'CMS_htt_boson_reso_met',
		'CMS_htt_boson_scale_met',
		'CMS_htt_eff_b',
		'CMS_htt_mistag_b',
		'CMS_scale_met_unclustered']
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
	parser.add_argument('-y', dest='year', help='data taking year',choices = ['2016','2017'], default = '2016')
	args = parser.parse_args()
	channel = args.channel
	year = args.year
	print "=== Decay channel: {0} ===".format(channel)
	print "=== year: {0} ===".format(year)

	path_Vienna = "/afs/hephy.at/data/higgs02/shapes/Vienna/htt_{0}.inputs-sm-Run{1}-ML.root".format(channel,year)
	path_KIT = "/afs/hephy.at/data/higgs02/shapes/KIT/htt_{0}.inputs-sm-Run{1}-ML.root".format(channel,year)

	Vienna_file = TFile.Open(path_Vienna,'read')
	Tdirectories_Vienna = [key for key in Vienna_file.GetListOfKeys()]

	KIT_file = TFile.Open(path_KIT,'read')
	Tdirectory_names_KIT = [key.GetName() for key in KIT_file.GetListOfKeys()]


	output = TFile.Open("/afs/hephy.at/user/j/jandrejkovic/public/forMarkus/datacards_new/vienna/htt_{0}.inputs-sm-Run{1}-ML.root".format(channel,year),'recreate')

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