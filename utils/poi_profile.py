#plot a poi from a scan in a results pickle and show the other profiled parameters too

# python poi_profile.py results.pkl parameter 

import sys,os
import numpy as np
import pickle
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 6})
fig,axs = plt.subplots(1,3,figsize=(28, 7),constrained_layout=True)


inputPkl = sys.argv[1]
which    = sys.argv[2]

with open(inputPkl,"rb") as fpkl: results = pickle.load(fpkl)
mode = "profiled" #  choose profiled or fixedd

dchi2 = results[which][mode]['dchi2']
p     = results[which][mode]['pvals']

axs[0].plot(p,dchi2,color='black',marker="o")
axs[0].set_xlabel(which)
axs[0].set_ylabel("$\Delta \chi^{2}$")
npoints = len(p)
npoi = len(list(results.keys()))
npredictions = len(results[which][mode]['predictions'])

jet = plt.cm.jet
colors   = jet(np.linspace(0, 1, npoi))
colors2  = jet(np.linspace(0, 1, npredictions))
styles = ["-","--",":"]

for i,poi in enumerate(results[which][mode]['otherpoi'].keys()):
    if poi == which : continue
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

#fig.tight_layout()

plt.show()
