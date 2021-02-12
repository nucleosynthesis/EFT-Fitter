from collections import OrderedDict as od

# H->multilepton (2016+2017): ttH

name = "HIG-18-019"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["ttH_hmultilepton"] = {"bestfit":0.96, "Up01Sigma":0.34, "Down01Sigma":0.31, "Up01SigmaExp":0.30, "Down01SigmaExp":0.27, "merged":False}

# Correlations
rho = od()


