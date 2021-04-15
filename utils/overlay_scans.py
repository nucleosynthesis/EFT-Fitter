# very quick tool to plot likelihood scans from txt files (1D only)
# python overlay_scans.py outname input1.txt [input2.txt input3.txt....] 
# set outname="show" to just show the results rather than save as pdf/png

import numpy
import matplotlib.pyplot as plt

import sys 
outp   = sys.argv[1]
inputs = sys.argv[2:]

MAXC2 = 10

for in_p in inputs:
 fi = open(in_p,"r")
 points  = []
 for i,line in enumerate(fi.readlines()): 
  if i==0: par  = line.split()[0]
  else:
   #if float(line.split()[1]) > MAXC2-0.001: continue
   points.append([float(line.split()[0]),float(line.split()[1])])

 points.sort()
 p  = [pp[0] for pp in points]
 c2 = [pp[1] for pp in points]
 plt.plot(p,c2,label=in_p.replace(".txt",""),marker="o",markersize=2.0)

plt.xlabel(par)
plt.ylim(0,MAXC2)
plt.ylabel("$\Delta \chi^{2}$")
plt.legend()
if outp=="show":
    plt.show()
else:
    print("Saved as ",outp+".pdf",outp+".png")
    plt.savefig(outp+".pdf")
    plt.savefig(outp+".png")
