from collections import OrderedDict as od

# H->bb (2016+2017): ttH

name = "HIG-18-030"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["ttH_hbb"] = {"bestfit":1.15, "Up01Sigma":0.32, "Down01Sigma":0.29, "Up01SigmaExp":0.32, "Down01SigmaExp":0.29, "merged":False}

# Correlations
rho = od()


