import matplotlib.pyplot as plt
from importlib import import_module
import sys 
import numpy 

mod = import_module(sys.argv[1])

lh = mod.X['likelihood']['likelihood']
npartotal = len(lh.getParameters())
xwid = int(npartotal**0.5)+1

total = xwid*xwid
ywid = xwid
if total-ywid > npartotal : ywid-=1

fig, axs = plt.subplots(xwid, ywid)
jindex=0
iindex=0
print("total pars = ",npartotal," xwid = ", xwid)
for p in range(lh.getNSplines()):
  if iindex>xwid-1:
    jindex+=1
    iindex = 0
  sp = lh.getSpline(p)
  par = sp.getParameters()
  if len(par)>1 : continue
  par = par[0]
  minp,maxp=sp.getMinMax(par)
  x = numpy.arange(minp,maxp,(maxp-minp)/30)
  y = [sp.evaluate({par:xx}) for xx in x]
  #print(iindex,jindex,p)
  x2 = sp.getPoints([par])
  y2 = sp.getF()
  axs[iindex,jindex].plot(x,y)
  axs[iindex,jindex].plot(x2,y2,marker="o",markersize=0.2,linestyle='none',color='red')
  axs[iindex,jindex].set_xlabel(par)
  iindex+=1

plt.show()

