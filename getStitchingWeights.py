import json

def main():

	with open("stitchConfig.json","r") as FSO:
		stitch_info = json.load(FSO)

	with open("tagMapping.json","r") as FSO:
		tag_map = json.load(FSO)

	for mergename in stitch_info:
		if  stitch_info[mergename]["stitching"]:
			xsecs =   [0]*5
			nevents = [0]*5
			for sample in stitch_info[mergename]["samples"]:
				nj = getJets(sample)
				xsecs[nj] = tag_map[sample]["xsec"]
				nevents[nj] = tag_map[sample]["nevents"]

			weights = calcStitchingWeights( stitch_info[mergename]["NNLO_xsec"], xsecs, nevents )
			print mergename
			print showWeightstring(weights)

def showWeightstring(weights):
	wstr = "("
	for i in xrange(5):
		if i == 0:
			wstr += "(NUP == {0})*{1}".format(i,weights[i])
		else:
			wstr += " + (NUP == {0})*{1}".format(i,weights[i])

	return wstr + ")"

def getJets(sample):
	idx = 0
	for i,j in enumerate(["","1","2","3","4"]):
		if j+"Jets" in sample:
			idx = i
	return idx

def calcStitchingWeights(nnlo_xsec, xsecs, nevents):

	corr = nnlo_xsec / xsecs[0]
	weights = [0]*5
	for i in xrange(len(xsecs)):
		if i == 0:
			weights[i] = corr / ( nevents[i] / xsecs[i] )
		else:
			weights[i] = corr / ( ( nevents[0] / xsecs[0] ) + ( nevents[i] / xsecs[i] )  )

	return weights

if __name__ == '__main__':
	main()