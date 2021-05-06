#plot a poi from a scan in a results pickle and show the other profiled parameters too

# python poi_profile.py results.pkl 

import sys,os
import numpy as np
import pickle
import matplotlib.pyplot as plt
fig,axs = plt.subplots(1,2,figsize=(12, 6))


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

jet = plt.cm.jet
colors = jet(np.linspace(0, 1, npoi))
styles = ["-","--",":"]
for i in range(npoi):
    if list(results.keys())[i] == which : continue 
    axs[1].plot(p,[results[which]["profiled"]['allpvals'][j][i] for j in range(npoints)],label=list(results.keys())[i],color=colors[i],linestyle=styles[i%3])

axs[0].set_xlabel(which)
axs[0].grid()
axs[1].set_ylabel("param value")
axs[1].grid()


plt.legend(bbox_to_anchor=(1.2, 1.05))
fig.tight_layout()
plt.show()
