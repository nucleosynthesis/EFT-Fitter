from collections import OrderedDict as od

# FCCee projections
name = "FCCee"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
# Numbers below are in % of the SM value (i.e 0.5 = 0.5% = 0.005)
X = od()
X["ZH_lep"]        = {"Up01SigmaExp":0.52 , "Down01SigmaExp":0.52, "merged":False}
X["ZH_had"]        = {"Up01SigmaExp":0.38 , "Down01SigmaExp":0.38, "merged":False}
X["ZH_hbb"]        = {"Up01SigmaExp":0.21 , "Down01SigmaExp":0.21, "merged":False}
X["ZH_hcc"]        = {"Up01SigmaExp":1.60 , "Down01SigmaExp":1.60 , "merged":False}
X["ZH_hss"]        = {"Up01SigmaExp":120. , "Down01SigmaExp":120.  , "merged":False}
X["ZH_hgluglu"]    = {"Up01SigmaExp":0.80 , "Down01SigmaExp":0.80 , "merged":False}
X["ZH_tautau"]     = {"Up01SigmaExp":0.58 , "Down01SigmaExp":0.58, "merged":False}
X["ZH_mumu"]       = {"Up01SigmaExp":11.0 , "Down01SigmaExp":11.0 , "merged":False}
X["ZH_hWW"]        = {"Up01SigmaExp":0.80 , "Down01SigmaExp":0.80 , "merged":False}
X["ZH_hZZ"]        = {"Up01SigmaExp":2.50 , "Down01SigmaExp":2.50  , "merged":False}
X["ZH_hgamgam"]    = {"Up01SigmaExp":3.60 , "Down01SigmaExp":3.60 , "merged":False}
X["ZH_hZgam"]      = {"Up01SigmaExp":11.8 , "Down01SigmaExp":11.8, "merged":False}

# Correlations: Assume same correlations as FCC-ee numbers

rho = od()
rho[("ZH_hss","ZH_hbb")] = 0.0797
rho[("ZH_hss","ZH_hcc")] = 0.0623 
rho[("ZH_hss","ZH_hgluglu")] = 0.0190
rho[("ZH_hbb","ZH_hcc")] = 0.0752 
rho[("ZH_hss","ZH_hgluglu")] = 0.0238
rho[("ZH_hcc","ZH_hgluglu")] = 0.0385 
#rho[("ZH_hbb","ZH_had")] = 0.99











