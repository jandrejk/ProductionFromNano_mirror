import shutil
import glob
import sys
import ROOT as R

channel = str(sys.argv[1])
outdir =  str(sys.argv[2])
rundir =  str(sys.argv[3])



files = glob.glob("{channel}-*.root".format(channel = channel) )

summary = ""

for filename in files:
	error = 0
	tmpfile = R.TFile(filename)
	if tmpfile.IsOpen() and not tmpfile.IsZombie():
		shutil.copyfile(filename, "/".join([outdir,filename]) )
		remotefile = R.TFile("/".join([outdir,filename]))
		if remotefile.IsOpen() and not remotefile.IsZombie():
			summary += "file:{0}#DONE#\n".format(filename)
		else: error = 1
	else: error = 2

	if error: summary += "file:{0}#ERROR{1}#\n".format(filename,error)

	tmpfile.Close()

with open("copy_status.log","w") as FSO:
	FSO.write(summary)

if os.getcwd() != rundir:
	shutil.copy("copy_status.log", "/".join([rundir,"copy_status.log"]) )