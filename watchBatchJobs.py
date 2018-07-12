import subprocess as sp
import os
import shlex
import argparse


def listify(string):

	l = []
	for s in string.split(" "):
		if s: l.append(s)
	return l

def usefullInfo(string):

	infolist = listify(string)

	return [ infolist[0], infolist[4] ]


def makeSummary( output ):

	summary = {}
	for line in output.splitlines():
		info = usefullInfo(line)
		
		if not summary.get( info[1], False ):
			summary[info[1]] = [ info[0] ]
		else:
			summary[info[1]].append(info[0])

	return summary


parser = argparse.ArgumentParser()
parser.add_argument('-v', dest='verbose', help='Verbose output', action="store_true")
args = parser.parse_args()
# user = "fspreitzer"
user = os.environ["USER"]

proc = sp.Popen(shlex.split('squeue -h -l -u {0}'.format(user) ) , stdout=sp.PIPE, stderr=sp.PIPE)
(out, err) = proc.communicate()

print user
for status, jobs in makeSummary(out).items():
	print "{0}: {1} jobs".format( status, len(jobs)  )
	if args.verbose:
		for job in jobs:
			print "\t", job

