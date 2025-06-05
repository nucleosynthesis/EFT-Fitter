# very quick tool to plot likelihood scans from txt files (1D only)
# python overlay_scans.py outname input1.txt [input2.txt input3.txt....] 
# set outname="show" to just show the results rather than save as pdf/png

import numpy
import matplotlib.pyplot as plt
import pickle
import sys 


from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--poi', dest='poi', help="Parameter of interest")
  parser.add_option('--out', dest='out', default='out',help="Output file name (without extension)")
  return parser.parse_args()
(opt,args) = get_options()

outp   = opt.out
inputs = args

from intervals import * 

for in_p in inputs:
 if in_p.endswith(".txt"):
  fi = open(in_p,"r")
  points  = []
  for i,line in enumerate(fi.readlines()): 
    if i==0: 
      par  = line.split()[0]
      opt.poi = par

    else:
   #if float(line.split()[1]) > MAXC2-0.001: continue
   # we assume they should be from 0 so reset why not
      points.append([float(line.split()[0]),float(line.split()[1])])
 else:
  with open(in_p,"rb") as fpkl: results = pickle.load(fpkl)
  mode = 'profiled' #  choose profiled or fixedd
  which = opt.poi
  dchi2 = results[which][mode]['dchi2']
  p     = results[which][mode]['pvals']
  points = [[pp,dd] for pp,dd in zip(p,dchi2)]

 #print("Off-setting to 0 and minimum")
 #min_c2 = min([pp[1] for pp in points])
 #points = [[pp[0],pp[1]-min_c2] for pp in points]
 points.sort()
 p  = [pp[0] for pp in points]
 c2 = [pp[1] for pp in points]
 lab = in_p.replace(".txt","")
 crossings1 = findCrossings(points,1)
 crossings4 = findCrossings(points,4)
 minx = findMin(points)
 print(in_p)
 print("%s ->\nMinimum at %s=%.4f"%(lab,opt.poi,minx))
 print("... crossings 1 at ",",".join(["%4.4f"%f for f in crossings1]))
 print("... crossings 4 at ",",".join(["%4.4f"%f for f in crossings4]))
 crossings1 = numpy.array(crossings1)
 
 pu  = mymin(abs(crossings1[crossings1>=minx]-minx))
 pd  = mymin(abs(crossings1[crossings1<=minx]-minx))
 print("uncert = +%4.4f -%4.4f"%(pu,pd))
 print(" ... +/- %4.4f"%(abs(pu+pd)/2.))
 lab+= " "+opt.poi+"=%.4f"%(minx)+"$^{+%.4f}$"%(pu)+"$_{-%4.4f}$"%(pd)
 plt.plot(p,c2,label=lab,marker="o",markersize=2.0)

 for x in crossings1: 
   plt.plot([x,x],[0,1],color='red',linestyle='--')
 for x in crossings4: 
   plt.plot([x,x],[0,4],color='red',linestyle='--')

plt.xlabel(opt.poi)
plt.ylim(0,MAXC2)
plt.ylabel("$\\Delta \\chi^{2}$")
plt.legend()

if outp=="show":
    plt.show()
else:
    print("Saved as ",outp+".pdf",outp+".png")
    plt.savefig(outp+".pdf")
    plt.savefig(outp+".png")
