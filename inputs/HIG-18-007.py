from collections import OrderedDict as od

# H->tautau (2016): VH lep

name = "HIG-18-007"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["WH_lep_htt"] = {"bestfit":3.6, "Up01Sigma":1.8, "Down01Sigma":1.6, "Up01SigmaExp":1.6, "Down01SigmaExp":1.4, "merged":False}
X["ZH_lep_htt"] = {"bestfit":1.4, "Up01Sigma":1.6, "Down01Sigma":1.5, "Up01SigmaExp":1.5, "Down01SigmaExp":1.3, "merged":False}

# Correlations: assumed zero correlation between measured params
rho = od()
rho[("WH_lep_htt","ZH_lep_htt")] = 0.


