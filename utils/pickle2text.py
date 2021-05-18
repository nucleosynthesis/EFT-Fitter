# script to covert picked results into a simple txt file 
# python pickle2text.py results.pkl parameter outname
import sys,os
import pickle
import numpy as np
from collections import OrderedDict as od
from importlib import import_module
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import sys 

inputPkl = sys.argv[1]
which    = sys.argv[2]
oname    = sys.argv[3]

with open(inputPkl,"rb") as fpkl: results = pickle.load(fpkl)
mode = "profiled" #  choose profiled or fixedd

ext=False

params = []
if ".py" in which:
  mod = import_module(which.replace(".py",""),package=None)
  for k in mod.pois.keys(): 
    if "freeze" in mod.pois[k].keys():
      if mod.pois[k]['freeze']: continue
    params.append(k)
  ext=True
   
else : params = [which]
print(params)
for par in params: 
    dchi2 = results[par][mode]['chi2']
    p     = results[par][mode]['pvals']

    if ext : onamev = par+"_"+oname
    else: onamev = oname
    out = open(onamev,"w")
    out.write(par)
    out.write(" chi2 ")
    out.write("\n")

    for e in range(len(p)):
      out.write(" %g  %g "%(p[e],dchi2[e]))
      out.write("\n")
    out.close()
    print("1d scan for ",par,"written to ",onamev)

