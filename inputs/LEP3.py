from collections import OrderedDict as od

# LEP3 projections
name = "LEP3"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
# Numbers below are in % of the SM value (i.e 0.5 = 0.5% = 0.005)
X = od()
X["ZH_lep"]        = {"Up01SigmaExp":1.14 , "Down01SigmaExp":1.14, "merged":False}
X["ZH_had"]        = {"Up01SigmaExp":0.84 , "Down01SigmaExp":0.84, "merged":False}
X["ZH_hbb"]        = {"Up01SigmaExp":0.462, "Down01SigmaExp":0.462, "merged":False}
X["ZH_hcc"]        = {"Up01SigmaExp":3.52 , "Down01SigmaExp":3.52 , "merged":False}
X["ZH_hss"]        = {"Up01SigmaExp":264.  , "Down01SigmaExp":264  , "merged":False}
X["ZH_hgluglu"]    = {"Up01SigmaExp":1.76 , "Down01SigmaExp":1.76 , "merged":False}
X["ZH_tautau"]     = {"Up01SigmaExp":1.276, "Down01SigmaExp":1.276, "merged":False}
X["ZH_mumu"]       = {"Up01SigmaExp":24.2 , "Down01SigmaExp":24.2 , "merged":False}
X["ZH_hWW"]        = {"Up01SigmaExp":1.76 , "Down01SigmaExp":1.76 , "merged":False}
X["ZH_hZZ"]        = {"Up01SigmaExp":5.5  , "Down01SigmaExp":5.5  , "merged":False}
X["ZH_hgamgam"]    = {"Up01SigmaExp":7.92 , "Down01SigmaExp":7.92 , "merged":False}
X["ZH_hZgam"]      = {"Up01SigmaExp":25.96, "Down01SigmaExp":25.96, "merged":False}

# Correlations: Assume same correlations as FCC-ee numbers

rho = od()
rho[("ZH_hss","ZH_hbb")] = 0.0797
rho[("ZH_hss","ZH_hcc")] = 0.0623 
rho[("ZH_hss","ZH_hgluglu")] = 0.0190
rho[("ZH_hbb","ZH_hcc")] = 0.0752 
rho[("ZH_hss","ZH_hgluglu")] = 0.0238
rho[("ZH_hcc","ZH_hgluglu")] = 0.0385 




