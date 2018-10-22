import shutil
import glob
import sys
import os
import ROOT as R

channel = str(sys.argv[1])
outdir =  str(sys.argv[2])
rundir =  str(sys.argv[3])



files = glob.glob("{channel}-*.root".format(channel = channel) )

summary = ""
print "#####   VALIDATING ######"
glob_error = 0
for filename in files:
	error = 0
	tmpfile = R.TFile(filename)
	if tmpfile.IsOpen() and not tmpfile.IsZombie() and not tmpfile.TestBit(R.TFile.kRecovered) and len(tmpfile.GetListOfKeys()) > 0:
		shutil.copyfile(filename, "/".join([outdir,filename]) )
		remotefile = R.TFile("/".join([outdir,filename]))
		if remotefile.IsOpen() and not remotefile.IsZombie():
			os.chmod( "/".join([outdir,filename]), 0777)
			summary += "file:{0}#DONE#\n".format(filename)
		else: error = 1
	else: error = 2

	if error:
		glob_error = error
		summary += "file:{0}#ERROR{1}#\n".format(filename,error)

	tmpfile.Close()

with open("copy_status.log","w") as FSO:
	FSO.write(summary)

if os.getcwd() != rundir:
	shutil.copy("copy_status.log", "/".join([rundir,"copy_status.log"]) )

sys.exit(glob_error)