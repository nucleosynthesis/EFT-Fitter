import sys
import ROOT
fi = ROOT.TFile.Open(sys.argv[1])
oname = sys.argv[2]
MAXc2 = 20

pick_branches = []
if len(sys.argv)>3:
  pick_branches = sys.argv[3:]
  print("selecting only branches",pick_branches)

# ToDo: add some options to allow relabling of the parameters
veto_brs = ["limit", "limitErr", "mh", "syst", "iToy", "iSeed", "iChannel", "t_cpu", "t_real", "quantileExpected","deltaNLL"]

def checkc2(tr,e):
  tr.GetEntry(e)
  if (2*tr.deltaNLL) > MAXc2: return False
  return True

def printPoint(tr, keys,e):
  tr.GetEntry(e)
  ret_str = ["%g"%(getattr(tr,k)) for k in keys]
  ret_str.append("%g"%(2*tr.deltaNLL))
  return " ".join(ret_str)

def printSPoint(pt):
  return " ".join(pt)

def getPoint(tr,keys,e):
  ret = ["%g"%(getattr(tr,k)) for k in keys]
  ret.append("%g"%(2*tr.deltaNLL))
  return ret
  

tr = fi.Get("limit")
keys  = tr.GetListOfBranches()

params = []
for k in keys: 
 nm = k.GetName()
 if len(pick_branches): 
   if not nm in pick_branches: continue
 else:
   if nm in veto_brs: continue
 params.append(nm)

out = open(oname,"w")
out.write(" ".join(params))
out.write(" chi2 ")
out.write("\n")

entries = tr.GetEntries()
points = []
for e in range(entries):
  if not checkc2(tr,e): continue
  pt = getPoint(tr,params,e)
  if pt in points: continue 
  points.append(pt)
  #print(points)
for pt in points:
  out.write(printSPoint(pt))
  out.write("\n")
out.close()
 

