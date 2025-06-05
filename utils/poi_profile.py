#plot a poi from a scan in a results pickle and show the other profiled parameters too

# python poi_profile.py results.pkl parameter 

import sys,os
import numpy as np
import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 6})
from intervals import * 

from optparse import OptionParser

def get_options():
  parser = OptionParser()
  parser.add_option('--pois', dest='pois', default='all', help="Parameters of interest")
  parser.add_option('--save', dest='save', default=False,action='store_true', help="Save the plot")
  parser.add_option('--show', dest='show', default=False,action='store_true', help="Show the plot")
  parser.add_option('--which', dest='which', default='profiled', help="profile/fixed")
  parser.add_option('--chi2max', dest='chi2max', type='float', default=-1, help="Set a maximum chi2 value to plot, default is -1 (no limit)")
  return parser.parse_args()
(opt,args) = get_options()

inputPkl = sys.argv[1]

with open(inputPkl,"rb") as fpkl: results = pickle.load(fpkl)
mode = opt.which #  choose profiled or fixedd

if opt.pois=="all":
    pois = ','.join(results.keys())
    pois = pois.split(",")
else:
    pois    = opt.pois.split(",")

jet = plt.cm.jet
styles = ["-","--",":"]

print("POIs: ",pois)

print(results['kappa_Z'].keys())
for i,which in enumerate(pois):
    fig,axs = plt.subplots(1,3,figsize=(28, 7),constrained_layout=True)
    dchi2 = results[which][mode]['dchi2']
    p     = results[which][mode]['pvals']

    # print the crossings etc 
    points = [[pp,dd] for pp,dd in zip(p,dchi2)]
    points.sort()
    crossings1 = findCrossings(points,1)
    crossings4 = findCrossings(points,4)
    minx = findMin(points)
    par = which 
    print("\n%s=%.4f"%(par,minx))
    print(" ... crossings 1 at ",",".join(["%4.4f"%f for f in crossings1]))
    print(" ... crossings 4 at ",",".join(["%4.4f"%f for f in crossings4]))
    crossings1 = np.array(crossings1)
 
    pu  = mymin(abs(crossings1[crossings1>=minx]-minx))
    pd  = mymin(abs(crossings1[crossings1<=minx]-minx))
    print(" uncert = +%4.4f -%4.4f"%(pu,pd))
    print("     ... +/- %4.4f"%(abs(pu+pd)/2.))
    ######## 

    axs[0].plot(p,dchi2,color='black',marker="o")
    axs[0].set_xlabel(which)
    axs[0].set_ylabel("$\Delta \chi^{2}$")
    npoints = len(p)
    npoi = len(list(results[which][mode]['otherpoi']))
    npredictions = len(results[which][mode]['predictions'])


    colors   = jet(np.linspace(0, 1, npoi))
    colors2  = jet(np.linspace(0, 1, npredictions))

    if opt.chi2max > 0:
        axs[0].set_ylim(0,opt.chi2max)
    for i,poi in enumerate(results[which][mode]['otherpoi'].keys()):
        if poi == which : 
            i-=1
            continue
        axs[1].plot(p,results[which][mode]['otherpoi'][poi],label=poi,color=colors[i],linestyle=styles[i%3])

    for i,label in enumerate(results[which][mode]["predictions"].keys()):
        axs[2].plot(p,results[which][mode]['predictions'][label],label=label,color=colors2[i],linestyle=styles[i%3])

    axs[0].set_xlabel(which)
    axs[0].grid()
    axs[1].set_xlabel(which)
    axs[1].set_ylabel("param value")
    axs[1].grid()
    axs[2].set_xlabel(which)
    axs[2].set_ylabel("prediction value")
    axs[2].grid()

    axs[1].legend(bbox_to_anchor=(1.1, 1.0))
    axs[2].legend(bbox_to_anchor=(1.1, 1.0),ncol=2)
    plt.margins(x=0.1)


    if opt.save:
        plt.savefig("poi_profile_%s_%s.pdf"%(which,mode),bbox_inches='tight')
        plt.savefig("poi_profile_%s_%s.png"%(which,mode),bbox_inches='tight')
                    
    elif opt.show:
        plt.show()
