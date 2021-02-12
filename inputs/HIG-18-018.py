from collections import OrderedDict as od

name = "HIG-18-018"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["ttH_hgg"] = {"bestfit":1.7, "Up01Sigma":0.6, "Down01Sigma":0.5, "Up01SigmaExp":0.6, "Down01SigmaExp":0.5, "merged":False}

# Correlations
rho = od()
