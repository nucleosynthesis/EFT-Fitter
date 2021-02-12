from collections import OrderedDict as od

# H->bb (2016+2017): VH lep

name = "HIG-18-016"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["WH_lep_hbb"] = {"bestfit":1.24, "Up01Sigma":0.38, "Down01Sigma":0.38, "Up01SigmaExp":0.38, "Down01SigmaExp":0.38, "merged":False}
X["ZH_lep_hbb"] = {"bestfit":0.88, "Up01Sigma":0.29, "Down01Sigma":0.29, "Up01SigmaExp":0.29, "Down01SigmaExp":0.29, "merged":False}

# Correlations: assumed zero correlation between measured params
rho = od()
rho[("WH_lep_hbb","ZH_lep_hbb")] = 0.

