import sys, os
from glob import glob


def getSamples():
    samples = set()
    for file in glob("samples/mc/*/*"):
        samples.add( file.replace(".txt","").split("/")[-1] )
    return samples

def mergeFragments():    
    if not os.path.exists("./MCpuHistograms/merged") :
        os.makedirs("./MCpuHistograms/merged")
    if not os.path.exists("./MCpuHistograms/usedHistos") :
        os.makedirs("./MCpuHistograms/usedHistos")

    for sample in getSamples() :
        #print key	
        save_name = './MCpuHistograms/merged/'+sample+'.root'
        if glob( 'MCpuHistograms/*'+ sample +'*'):
            cmd = 'hadd '+save_name+' MCpuHistograms/*'+ sample +'*'
            mvcmd = 'mv ./MCpuHistograms/*'+sample+'*root ./MCpuHistograms/usedHistos/'
            os.system(cmd)
            os.system(mvcmd)

if __name__ == '__main__':
    mergeFragments()






