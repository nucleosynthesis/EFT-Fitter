import numpy
import matplotlib.pyplot as plt

import sys 
outp   = sys.argv[1]
inputs = sys.argv[2:]

MAXC2 = 100

for in_p in inputs:
 fi = open(in_p,"r")
 points  = []
 for i,line in enumerate(fi.readlines()): 
  if i==0: par  = line.split()[0]
  else:
   if float(line.split()[1]) > MAXC2: continue
   points.append([float(line.split()[0]),float(line.split()[1])])

 points.sort()
 p  = [pp[0] for pp in points]
 c2 = [pp[1] for pp in points]
 plt.plot(p,c2,label=in_p.replace(".txt",""),marker="o",markersize=2.0)

plt.xlabel(par)
plt.ylabel("$\Delta \chi^{2}$")
plt.legend()
if outp=="show":
    plt.show()
else:
    print("Saved as ",outp+".pdf",outp+".png")
    plt.savefig(outp+".pdf")
    plt.savefig(outp+".png")
