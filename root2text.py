import sys
import ROOT
fi = ROOT.TFile.Open(sys.argv[1])
oname = sys.argv[2]

# ToDo: add some options to allow relabling of the parameters
veto_brs = ["limit", "limitErr", "mh", "syst", "iToy", "iSeed", "iChannel", "t_cpu", "t_real", "quantileExpected","deltaNLL"]

def printPoint(tr, keys,e):
  tr.GetEntry(e)
  ret_str = ["%g"%(getattr(tr,k)) for k in keys]
  ret_str.append("%g"%(2*tr.deltaNLL))
  return " ".join(ret_str)

tr = fi.Get("limit")
keys  = tr.GetListOfBranches()

params = []
for k in keys: 
 nm = k.GetName()
 if nm in veto_brs: continue
 params.append(nm)

out = open(oname,"w")
out.write(" ".join(params))
out.write(" chi2 ")
out.write("\n")

entries = tr.GetEntries()
for e in range(entries):
  out.write(printPoint(tr,params,e))
  out.write("\n")
out.close()
 

