import shlex
import os
import json
import glob
import subprocess as sp
import argparse

def main():

  parser = argparse.ArgumentParser()
  parser.add_argument('-r', dest='resubmit', help='Resubmit failed jobs', action="store_true")

  args = parser.parse_args()

  B = Bookkeeping()
  B.printStatus()

  if args.resubmit:
    B.resubmitFailedjobs()


class Bookkeeping():

  def __init__(self):
    self.system = "lxbatch"
    self.cwd = os.getcwd()

    self.log = {}
    if os.path.exists("submit_log.log"):
      with open("submit_log.log","r") as FSO:
        self.log = json.load(FSO)

     
    self.runningJobs = self.getRunningJobs()
    self.summary = {}
    self.failed_paths = []
    self.outdir = self.log[self.system].pop("outdir", "")

    for sample in self.log[self.system]:
      with open( glob.glob( "samples/*/*/{0}*".format(sample) )[0], "r" ) as FSO:
        ntotal = len(FSO.read().splitlines() )
      for channel in self.log[self.system][sample]:
        for shift in self.log[self.system][sample][channel]:


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

          for file in glob.glob(  "{0}/v1/{1}/{2}-{3}*root".format(self.outdir, sample, channel, shift) ):
            if self.log[self.system][sample][channel][shift]["submit_time"] < os.path.getmtime(file):
              self.summary[sample][shift][channel]["finished_files"].append(file)

          self.summary[sample][shift][channel]["total"] = ntotal
          self.summary[sample][shift][channel]["finished"] = len(self.summary[sample][shift][channel]["finished_files"])

    self.matchRunInfo()
    self.getFullStatus()

  def getRunningJobs(self):
    if self.system == "lxbatch":

      runningJobs = {}
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

      return runningJobs

  def matchRunInfo(self):

    for sample in self.summary:
      for shift in self.summary[sample]:
        for channel in self.summary[sample][shift]:

          fjob = "{0}#{1}-{2}".format( sample,channel,shift )

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

          if (f + r + p) < t:
            self.matchRundirs(sample, shift, channel)
          if f == t:
            self.log[self.system][sample][channel][shift]["status"] = "DONE"

  def matchRundirs(self, sample, shift, channel):

    rundirs = {}
    for configBall in glob.glob("out/v1/{sample}/rundir_{channel}_{shift}*/configBall.json".format(sample=sample,channel=channel,shift=shift)):
      with open(configBall,"r") as FSO:
        rundir = "/".join(configBall.split("/")[:-1])
        run_file = json.load(FSO)["file"].split("/")[-1]
        jobname = "{0}#{1}-{2}_{3}".format( sample, channel, shift, run_file )

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

    for failed in self.failed_paths:
      if self.system == "lxbatch":
        os.chdir( "/".join([self.cwd, failed[0]]) )
        os.system( " bsub -q 1nd -J {0} submit.sh".format( failed[1] ) )

    print "_"*105



  def printStatus(self):

    samples = self.summary.keys()
    samples.sort()

    print "_"*110
    for sample in samples:
      print "\033[1m" + sample + "\033[0m","\n"
      for shift in self.summary[sample]:
        line = {"et":" "*30,"mt":" "*30,"tt":" "*30}
        for channel in self.summary[sample][shift]:



          if self.log[self.system][sample][channel][shift]["status"] == "DONE":
            line[channel] = " "*10 +"{0}: (      \033[1;32mFinished\033[0m     )".format(channel)
          else:

            finished = self.summary[sample][shift][channel]["finished"]
            total = self.summary[sample][shift][channel]["total"]
            r = self.summary[sample][shift][channel]["running"]
            p = self.summary[sample][shift][channel]["pending"]
            u = self.summary[sample][shift][channel]["failed"]

            ft = "{0}({1}/{2}/{3}/{4}/{5})".format( channel+": ",cS(u,"r") ,cS(p,"y"),cS(r,"b"), cS(finished,"g"),cS(total,"") )
            line[channel] = " "*(85-len(ft)) + ft

        print "{0}  {1} {2} {3}".format(" "*(15 - len(shift)) + shift, line["et"], line["mt"], line["tt"])
      print "_"*110


def cS(string, color):
  s = str(string)

  if not color: c = "\033[1;30m"
  if color == "r": c = "\033[1;31m"
  if color == "g": c = "\033[1;32m"
  if color == "y": c = "\033[1;33m"
  if color == "b": c = "\033[1;34m"

  return c+" "*(3 - len(s)) + s + "\033[0m"
    


if __name__ == '__main__':
  main()