from ROOT import TFile, gDirectory, TH1F
import root_pandas as rp

path_Vienna_emb = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/out/emb_test/mt-NOMINAL_all.root"
path_emb = "/afs/hephy.at/data/higgs01/v9/mt-NOMINAL_ntuple_EMB.root"

df_new = rp.read_root(path_Vienna_emb, key= "TauCheck",ignore='addlepton_p4')

print df_new['evt']

# new = TFile.Open(path_Vienna_emb,'read')
# emb = TFile.Open(path_emb,'read')

# tree = new.Get("TauCheck")
# events = tree.GetBranch('evt')
# l = tree.GetLeaf("evt")


# event_ids = []#[int(l.GetValue()) ]
# for it in xrange(l.GetBranch().GetEntries()) :
# 	l.GetBranch().GetEntry(it)
# 	event_ids.append(int(l.GetValue()))

#df.to_root("{0}_{1}.root".format(name,channel), key="TauCheck" )
