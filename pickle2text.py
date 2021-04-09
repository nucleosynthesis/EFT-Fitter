import pickle
import numpy as np
from collections import OrderedDict as od

import sys 

inputPkl = sys.argv[1]
which    = sys.argv[2]
oname    = sys.argv[3]

with open(inputPkl,"rb") as fpkl: results = pickle.load(fpkl)
mode = "profiled" #  profiled
#mode = "fixed" #  profiled

dchi2 = results[which][mode]['dchi2']
p     = results[which][mode]['pvals']

out = open(oname,"w")
out.write(which)
out.write(" chi2 ")
out.write("\n")

for e in range(len(p)):
  out.write(" %g  %g "%(p[e],dchi2[e]))
  out.write("\n")
out.close()

