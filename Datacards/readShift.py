from ROOT import TFile, gDirectory, TH1F
import argparse
import os
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
	out_histo.Delete()

def CheckDir (directory) :
	if not os.path.exists(directory):
    		os.makedirs(directory)

def PrintInputSettings (args) :
	print "==================== INPUT SETTINGS ===================="
	print "%-15s : %-20s" % ("decay channel", args.channel)
	print "%-15s : %-20s" % ("year", args.year)
	print "%-15s : %-20s" % ("1st input path", args.path1)
	print "%-15s : %-20s" % ("2nd input path", args.path2)
	print "%-15s : %-20s" % ("output path", args.outdirectory)
	print "%-15s : %-20s" % ("simple copy", args.copy)
	print "========================================================"
	

def ApplyShifts (h_key,d,d_KIT,KIT_histo_names,nominal_vienna,nominal_KIT,loop) :	
	if (h_key.GetName() in KIT_histo_names) : # Check if shape exists in KIT - if True copy Datacard from KIT
		datacard_shift_KIT = d_KIT.Get(h_key.GetName())
		datacard_shift_Vienna = nominal_vienna.Clone() 
		datacard_shift_Vienna.SetName(h_key.GetName())
		datacard_shift_Vienna.SetTitle(h_key.GetName())

		datacard_shift_Vienna.Multiply(datacard_shift_KIT)
		datacard_shift_Vienna.Divide(nominal_KIT)
		# Write shifted histogram
		datacard_shift_Vienna.SetDirectory(0)
		datacard_shift_Vienna.Write()
		# release histos from memory
		datacard_shift_KIT.Delete()
		datacard_shift_Vienna.Delete()	
		loop += 1
	else : # copy Datacard from Vienna
		WriteHisto(directory=d, Histkey=h_key)
	return loop

def CopyShifts (h_key,KIT_histo_names,loop,d_KIT,d) :
	if MissingShift(name=h_key.GetName()) :				
		# check if Datacard exists in KIT file
		if (h_key.GetName() in KIT_histo_names) : # if yes copy Datacard from KIT
			loop += 1
			WriteHisto(directory=d_KIT, Histkey=h_key)
		else : # copy Datacard from Vienna
			WriteHisto(directory=d, Histkey=h_key)
			
					
	else : # Datacard is not missing, therefore take it from Vienna
		WriteHisto(directory=d, Histkey=h_key)
	return loop

def main () :
	parser = argparse.ArgumentParser()
	parser.add_argument('-c', dest='channel', help='Dataset channel',choices = ['mt','et','tt'], default = 'mt')
	parser.add_argument('-y', dest='year', help='data taking year',choices = ['2016','2017'], default = '2016')
	parser.add_argument('-o', dest='outdirectory', help='output directory', default = '/afs/hephy.at/user/j/jandrejkovic/public/forMarkus/datacards_2018_12_19/vienna/')
	parser.add_argument('-in1', dest='path1', help='directory of Vienna shape datacards', default = '/afs/hephy.at/data/higgs02/shapes/Vienna')
	parser.add_argument('-in2', dest='path2', help='directory of KIT shape datacards', default = '/afs/hephy.at/data/higgs02/shapes/KIT')
	parser.add_argument('--copy', dest='copy', help='simply copy datacards from KIT', action = "store_true")
	
	args = parser.parse_args()
	
	channel = args.channel
	year = args.year
	directory = args.outdirectory
	CheckDir(directory=directory)
	
	PrintInputSettings(args=args)
	

	path_Vienna = "{2}/htt_{0}.inputs-sm-Run{1}-ML.root".format(channel,year,args.path1)
	path_KIT = "{2}/htt_{0}.inputs-sm-Run{1}-ML.root".format(channel,year,args.path2)

	Vienna_file = TFile.Open(path_Vienna,'read')
	Tdirectories_Vienna = [key for key in Vienna_file.GetListOfKeys()]

	KIT_file = TFile.Open(path_KIT,'read')
	Tdirectory_names_KIT = [key.GetName() for key in KIT_file.GetListOfKeys()]

	output = TFile.Open("{2}/htt_{0}.inputs-sm-Run{1}-ML.root".format(channel,year,directory),'recreate')
	
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
			
			d_KIT = KIT_file.Get(Tdir.GetName())
			KIT_histo_names = [key.GetName() for key in d_KIT.GetListOfKeys()]

			hist_dict = dict(zip(histo_name_keys,[k.GetName() for k in histo_name_keys]))

			# loop over all shapes in the given Tdirectory
			for h_key in sorted(hist_dict,key=hist_dict.get) :
				
				if args.copy : # just copy the missing datacards from KIT
					loop = CopyShifts(h_key=h_key, KIT_histo_names=KIT_histo_names, loop=loop, d_KIT=d_KIT, d=d)
				else : # apply shifts from KIT datacards to nominal Vienna datacards to get the shifted Vienna datacards
					if ('CMS' in h_key.GetName()) == False : # check if it is the nominal sample
						WriteHisto(directory=d, Histkey=h_key)
						nominal_vienna = d.Get(h_key.GetName())
						nominal_KIT    = d_KIT.Get(h_key.GetName())
					else :
						loop = ApplyShifts(h_key=h_key, d=d, d_KIT=d_KIT, KIT_histo_names=KIT_histo_names,nominal_vienna=nominal_vienna,nominal_KIT=nominal_KIT,loop=loop)

				counter += 1
				
		else : # Tdirectory only in Vienna, therefore copy all Datacards from Vienna that are inside this TDirectory
			for h_key in histo_name_keys[:2] :
				WriteHisto(directory=d, Histkey=h_key)
				counter +=1
				
					
	print 'No. of Datacards: ', counter
	print "No. of Datacards taken/shifted from KIT: ", loop
if __name__ == '__main__':
  	main()		
	print "---END--- readShift.py"
