from collections import OrderedDict as od

# 2016: H->WW analysis (Stage 0 only)
name = "HIG-16-042"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["ggH_hww"] = {"bestfit":1.24, "Up01Sigma":0.2, "Down01Sigma":0.26, "Up01SigmaExp":0.20, "Down01SigmaExp":0.26, "merged":False}
X["qqH_hww"] = {"bestfit":0.24, "Up01Sigma":0.74, "Down01Sigma":0.74, "Up01SigmaExp":0.74, "Down01SigmaExp":0.74, "merged":False}
X["qqH_VH2JET_hww"] = {"bestfit":12.88, "Up01Sigma":4.92, "Down01Sigma":4.89, "Up01SigmaExp":4.92, "Down01SigmaExp":4.89, "merged":False}
X["WH_lep_hww"] = {"bestfit":1.8, "Up01Sigma":1.96, "Down01Sigma":1.64, "Up01SigmaExp":1.96, "Down01SigmaExp":1.64, "merged":False}
X["ZH_lep_hww"] = {"bestfit":0.71, "Up01Sigma":1.68, "Down01Sigma":1.68, "Up01SigmaExp":1.68, "Down01SigmaExp":1.68, "merged":False}

# Correlations: no correlation info available, assume same as HIG-19-015
rho = od()
rho[("ggH_hww","qqH_hww")] = -0.32
rho[("ggH_hww","qqH_VH2JET_hww")] = -0.2
rho[("ggH_hww","WH_lep_hww")] = 0.01
rho[("ggH_hww","ZH_lep_hww")] = 0.01

rho[("qqH_hww","qqH_VH2JET_hww")] = -0.03
rho[("qqH_hww","WH_lep_hww")] = -0.01
rho[("qqH_hww","ZH_lep_hww")] = 0.0

rho[("qqH_VH2JET_hww","WH_lep_hww")] = 0.02
rho[("qqH_VH2JET_hww","ZH_lep_hww")] = 0.0

rho[("WH_lep_hww","ZH_lep_hww")] = -0.16

