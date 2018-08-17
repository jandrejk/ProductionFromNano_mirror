import shlex
import os
import json
import sys
import glob
import shutil
import subprocess as sp
import argparse
from runUtils import checkProxy, checkTokens, getSystem, getHeplxPublicFolder

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('-r', dest='resubmit', help='Resubmit failed jobs', action="store_true")
  parser.add_argument('-v', dest='verbose',  help='Print rundirs of failed jobs', action="store_true")

  args = parser.parse_args()

  if not checkTokens(): sys.exit()

  B = Bookkeeping(args.verbose)
  B.printStatus()

  if args.resubmit:
    B.resubmitFailedjobs()


class Bookkeeping():

  def __init__(self, verbose = False):

    self.verbose = verbose

    host = getSystem()
    if "hephy" in host: 
      self.system = "hephybatch"
    elif "cern" in host: 
      self.system = "lxbatch"
    else:
      print "Dont know host you are on. Only heplx and lxplus!"   
      sys.exit()
      
    self.cwd = os.getcwd()

    self.log = {}
    self.logpath = "/".join([ getHeplxPublicFolder(),"submit_log.log" ])
    if os.path.exists(self.logpath):
      with open(self.logpath,"r") as FSO:
        self.log = json.load(FSO)

     
    self.runningJobs = self.getRunningJobs()
    self.summary = {}
    self.failed_paths = []
    self.outdir = "/afs/hephy.at/data/higgs01"

    for sample in self.log:
      with open( glob.glob( "samples/*/*/{0}.txt".format(sample) )[0], "r" ) as FSO:
        ntotal = len(FSO.read().splitlines() )
      for channel in self.log[sample]:
        for shift in self.log[sample][channel]:


          if not self.summary.get(sample,False): self.summary[sample] = {}
          if not self.summary[sample].get(shift,False): self.summary[sample][shift] = {}
          if not self.summary[sample][shift].get(channel,False):
            self.summary[sample][shift][channel] = {"total":0,
                                                    "finished":0,
                                                    "finished_files":[],
                                                    "running":0,
                                                    "running_files":[],
                                                    "pending":0,
                                                    "unknown":0,
                                                    "failed":0,
                                                    "failed_files":[],
                                                    "jobids":[] }

          for file in glob.glob(  "{0}/{1}/{2}-{3}*root".format(self.outdir, sample, channel, shift) ):
            if self.log[sample][channel][shift]["submit_time"] < os.path.getmtime(file):
              self.summary[sample][shift][channel]["finished_files"].append(file)

          self.summary[sample][shift][channel]["total"] = ntotal
          self.summary[sample][shift][channel]["finished"] = len(self.summary[sample][shift][channel]["finished_files"])

    self.matchRunInfo()
    self.getFullStatus()

  def __del__(self):
    with open(self.logpath,"w") as FSO:
      json.dump(self.log, FSO, indent = 2)

  def getRunningJobs(self):
    runningJobs = {}
    if self.system == "lxbatch":


      proc = sp.Popen( shlex.split('bjobs -UF'), stdout=sp.PIPE )
      (out, err) = proc.communicate()

      for entry in out.splitlines():
        if 'Job <' in entry:
          jobinfo = entry.replace('>','').replace(" ","").split(',')
          for i,part in enumerate(jobinfo):
            if "<" in part:
              jobinfo[i] = tuple( part.split("<") )

          jobinfo = dict(jobinfo)
          jobid = jobinfo.pop("Job") 
          runningJobs[ jobid ] = jobinfo

    elif self.system == "hephybatch":
      runningJobs = {}
      statusmap = {
        "RUNN":"RUN",
        "PEND":"PEND"
      }
      proc = sp.Popen( shlex.split('squeue -u {0} -l --format="%.18i#%.150j#%.4T" -h '.format(os.environ["USER"]) ), stdout=sp.PIPE )
      (out, err) = proc.communicate()
      for entry in out.splitlines():
        jobinfo = entry.replace(" ","").split("#")
        if len(jobinfo) != 3: continue

        runningJobs[jobinfo[0]] = {}
        runningJobs[jobinfo[0]]["JobName"] = jobinfo[1]
        runningJobs[jobinfo[0]]["Status"] =  statusmap.get(jobinfo[2], "UNK")



    return runningJobs

  def matchRunInfo(self):

    for sample in self.summary:
      for shift in self.summary[sample]:
        for channel in self.summary[sample][shift]:

          fjob = "{0}+{1}-{2}".format( sample,channel,shift )

          for ojob in self.runningJobs:
            if fjob in self.runningJobs[ojob]["JobName"]:
              self.summary[sample][shift][channel]["jobids"].append( ojob )
              self.summary[sample][shift][channel]["running_files"].append( self.runningJobs[ojob]["JobName"] )
              if self.runningJobs[ojob]["Status"] == "RUN": self.summary[sample][shift][channel]["running"] += 1
              elif self.runningJobs[ojob]["Status"] == "PEND": self.summary[sample][shift][channel]["pending"] += 1
              else: self.summary[sample][shift][channel]["unknown"] += 1

  def getFullStatus(self):

    for sample in self.summary:
      for shift in self.summary[sample]:
        for channel in self.summary[sample][shift]:

          t = self.summary[sample][shift][channel]["total"]
          f = self.summary[sample][shift][channel]["finished"]
          r = self.summary[sample][shift][channel]["running"]
          p = self.summary[sample][shift][channel]["pending"]

          if (f + r + p) < t and self.log[sample][channel][shift]["site"] == self.system:
            self.matchRundirs(sample, shift, channel)
          if f == t and self.log[sample][channel][shift]["status"] != "DONE":
            self.log[sample][channel][shift]["status"] = "MERGE"

  def matchRundirs(self, sample, shift, channel):

    rundirs = {}
    for configBall in glob.glob("out/{sample}/rundir_{channel}_{shift}*/configBall.json".format(sample=sample,channel=channel,shift=shift)):
      with open(configBall,"r") as FSO:
        rundir = "/".join(configBall.split("/")[:-1])
        run_file = json.load(FSO)["file"].split("/")[-1]
        jobname = "{0}+{1}-{2}_{3}".format( sample, channel, shift, run_file )

      match = False
      for file in self.summary[sample][shift][channel]["finished_files"]:
        if run_file in file: match = True

      for file in self.summary[sample][shift][channel]["running_files"]:
        if run_file in file: match = True        

      if not match:
        self.summary[sample][shift][channel]["failed_files"].append( (rundir,jobname) )

    self.summary[sample][shift][channel]["failed"] = len(self.summary[sample][shift][channel]["failed_files"])
    self.failed_paths += self.summary[sample][shift][channel]["failed_files"]

  def resubmitFailedjobs(self):

    print "_"*105
    print "Resubmitting failed jobs"
    print "_"*105
    if not checkProxy(): return

    for failed in self.failed_paths:
      os.chdir( "/".join([self.cwd, failed[0]]) )
      if os.path.exists("proxy"):
          shutil.rmtree("proxy")
      shutil.copytree("/".join([self.cwd,"proxy" ]), "proxy")

      if self.system == "lxbatch":
        os.system( " bsub -q 1nd -J {0} submit.sh".format( failed[1] ) )

      if self.system == "hephybatch":
        os.system( "sbatch submit.sh" )        

    print "_"*105



  def printStatus(self):

    samples = self.summary.keys()
    samples.sort()

    for sample in samples:
      should_display = False
      print_summary = ""
      print_summary += "\n{0}\033[1m{1}\033[0m{0}\n\n".format( " "*((110 - len(sample))/2), sample )
      print_summary += "{0}_{1} ET {1}_{1} MT {1}_{1} TT {1}_\n".format("_"*16, "_"*13)
      print_summary += "{0}|{1}|{1}|{1}|\n".format(" "*16, " "*30)
      for shift in self.summary[sample]:
        line = {"et":" "*28,"mt":" "*28,"tt":" "*28}
        for channel in self.summary[sample][shift]:
          
          if self.log[sample][channel][shift]["status"] == "DONE":
            line[channel] = "         \033[1;32mFinished\033[0m           "

          elif self.log[sample][channel][shift]["status"] == "MERGE":
            should_display = True
            line[channel] = "       \033[1;33mReady to merge\033[0m       "

          else:
            should_display = True
            finished = self.summary[sample][shift][channel]["finished"]
            total = self.summary[sample][shift][channel]["total"]
            r = self.summary[sample][shift][channel]["running"]
            p = self.summary[sample][shift][channel]["pending"]
            u = self.summary[sample][shift][channel]["failed"]

            if self.log[sample][channel][shift]["site"] != self.system:
              r = p = u = "?"

            if self.verbose:
              for ff in self.summary[sample][shift][channel]["failed_files"]:
                print_summary += ff[0]+"\n"

            ft = "{0}{1}/{2}/{3}/{4}/{5} {0}".format( " "*4,cS(u,"r") ,cS(p,"b"),cS(r,"y"), cS(finished,"g"),cS(total,"") )
            line[channel] = ft

        print_summary += "{0} | {1} | {2} | {3} |\n".format(" "*(15 - len(shift)) + shift, line["et"], line["mt"], line["tt"])

      print_summary += "{0}|{1}|{1}|{1}|\n".format(" "*16, " "*30)
      print_summary += "="*110 + "\n"

      if should_display or self.verbose: print print_summary



def cS(string, color):
  
  if string == 0: s = " "
  else: s = str(string)

  if not color: c = "\033[1;30m"
  if color == "r": c = "\033[1;31m"
  if color == "g": c = "\033[1;32m"
  if color == "y": c = "\033[1;33m"
  if color == "b": c = "\033[1;34m"

  return c+" "*(3 - len(s)) + s + "\033[0m"
    

if __name__ == '__main__':
  main()
