# check your inputs with this script 
# python checkinputs.py module.input

# if the input is a data (eg likelhood scan or best fits) it will plot those
# if the input is a theory uncerts, it will plot the uncertainties and correlations 

import matplotlib.pyplot as plt
import matplotlib as mpl
from importlib import import_module
import sys,os
import numpy 
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
mod = import_module(sys.argv[1],package=None)

if hasattr(mod,"X"):
    if "likelihood" in mod.X.keys():
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
          axs[iindex,jindex].plot(x,y,label="Interp",color="gray",lw=0.8)
          axs[iindex,jindex].scatter(x2,y2,marker="o",s=1.0,color='red',label="Input")
          axs[iindex,jindex].set_xlabel(par)
          iindex+=1
        axs[0,0].legend()
        plt.show()
    else:  
      fig,axs = plt.subplots(1,2)
      objects = mod.X.keys()
      y_pos  = numpy.arange(len(objects))
      bestfit = [ mod.X[k]['bestfit'] for k in objects ]
      sigma   = [ 0.5*(mod.X[k]['Up01Sigma']+mod.X[k]['Down01Sigma']) for k in objects ]
      axs[0].errorbar(bestfit, y_pos,xerr=sigma,linestyle='none',marker='o',markersize=3,color='black')
      axs[0].set_xlabel('Measured Values.')
      axs[0].set_yticks(y_pos)
      axs[0].set_yticklabels(objects)
      data = []
      if hasattr(mod,"rho"):
        print("found correlations")
        for i,k in enumerate(objects):
         row = []
         for j,l in enumerate(objects):
           if i==j : row.append(1.)
           elif i>j: row.append(mod.rho[(l,k)])
           else: row.append(mod.rho[(k,l)])
         data.append(row)
      else: 
        for i,k in enumerate(objects):
         row = []
         for j,l in enumerate(objects):
           if i==j : row.append(1.)
           else: row.append(0.)
         data.append(row)
      pos = axs[1].imshow(data,interpolation='none',cmap=plt.get_cmap('seismic'), norm=mpl.colors.Normalize(vmin=-1, vmax=1))
      fig.colorbar(pos,ax=axs[1])
      axs[1].set_xticks(y_pos)
      axs[1].set_yticks(y_pos)
      axs[1].set_xticklabels(objects,rotation=90)
      axs[1].set_yticklabels(objects)
      plt.show()

elif hasattr(mod,"TH"):
  #import matplotlib.pyplot as plt; plt.rcdefaults()
  fig,axs = plt.subplots(1,2)
  objects = mod.TH.keys()
  y_pos  = numpy.arange(len(objects))
  uncert = [ mod.TH[k] for k in objects ]
  axs[0].barh(y_pos, uncert, align='center', alpha=0.5)
  axs[0].set_xlabel('Theory Uncert.')
  axs[0].set_yticks(y_pos)
  axs[0].set_yticklabels(objects)
  #plt.title('Category')
  
  data = []
  if hasattr(mod,"rhoTH"):
    print("found correlations")
    for i,k in enumerate(objects):
     row = []
     for j,l in enumerate(objects):
       if i==j : row.append(1.)
       elif i>j: row.append(mod.rhoTH[(l,k)])
       else: row.append(mod.rhoTH[(k,l)])
     data.append(row)
  else: 
    for i,k in enumerate(objects):
     row = []
     for j,l in enumerate(objects):
       if i==j : row.append(1.)
       else: row.append(0.)
     data.append(row)
  pos = axs[1].imshow(data,interpolation='none',cmap=plt.get_cmap('seismic'), norm=mpl.colors.Normalize(vmin=-1, vmax=1))
  fig.colorbar(pos,ax=axs[1])
  axs[1].set_xticks(y_pos)
  axs[1].set_yticks(y_pos)
  axs[1].set_xticklabels(objects,rotation=90)
  axs[1].set_yticklabels(objects)
  plt.show()

else: 
      print("Incorrect data format for",sys.argv[1])
      sys.exit()
