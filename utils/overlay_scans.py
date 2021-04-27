# very quick tool to plot likelihood scans from txt files (1D only)
# python overlay_scans.py outname input1.txt [input2.txt input3.txt....] 
# set outname="show" to just show the results rather than save as pdf/png

import numpy
import matplotlib.pyplot as plt

import sys 
outp   = sys.argv[1]
inputs = sys.argv[2:]

MAXC2 = 10

def findCrossings(pts,c):
  # start from the left and find the point where we cross c
  # pts are poi,chi2
  xings =[]
  for i,pt in enumerate(pts):     
    if i==len(pts)-1: break
    pt2 = pts[i+1]
    #print(pt,pt2,c)
    #print(pt[1] < c and pt2[1] > c)
    #print(pt[1] > c and pt2[1] < c)
    if (pt[1] < c and pt2[1] > c) or (pt[1] > c and pt2[1] < c): 
      m = (pt2[1]-pt[1])/(pt2[0]-pt[0])
      b = pt[1] - m*pt[0]
      x = (c-b)/m 
      xings.append(x)
  return xings

def findMin(pts):
  pts2 = [[pt[1],pt[0]] for pt in pts]
  pts2.sort()
  return pts2[0][1]

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
 lab = in_p.replace(".txt","")
 crossings1 = findCrossings(points,1)
 crossings4 = findCrossings(points,4)
 print("%s, Minimum at %s=%.3f"%(lab,par,findMin(points)))
 print("... crossings 1 at ",",".join(["%4.3f"%f for f in crossings1]))
 print("... crossings 4 at ",",".join(["%4.3f"%f for f in crossings4]))
 plt.plot(p,c2,label=lab,marker="o",markersize=2.0)

 for x in crossings1: 
   plt.plot([x,x],[0,1],color='red',linestyle='--')
 for x in crossings4: 
   plt.plot([x,x],[0,4],color='red',linestyle='--')

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
