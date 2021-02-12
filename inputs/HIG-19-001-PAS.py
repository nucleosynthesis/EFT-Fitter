from collections import OrderedDict as od

# TO CHECK:
# * qqH_REST: measurement leading to negative on diagonal of Vinv, omitting for now

name = "HIG-19-001-PAS"

# Bestfit + uncertainties
X = od()
# ggH bins
X["ggH_0J_PTH_0_10_hzz"] = {"bestfit":0.87, "Up01Sigma":0.28, "Down01Sigma":0.25, "Up01SigmaExp":0.28, "Down01SigmaExp":0.25, "merged":False}
X["ggH_0J_PTH_GT10_hzz"] = {"bestfit":1.06, "Up01Sigma":0.18, "Down01Sigma":0.17, "Up01SigmaExp":0.18, "Down01SigmaExp":0.17, "merged":False}
X["ggH_1J_PTH_0_60_hzz"] = {"bestfit":0.78, "Up01Sigma":0.48, "Down01Sigma":0.51, "Up01SigmaExp":0.48, "Down01SigmaExp":0.51, "merged":False}
X["ggH_1J_PTH_60_120_hzz"] = {"bestfit":0.82, "Up01Sigma":0.41, "Down01Sigma":0.51, "Up01SigmaExp":0.41, "Down01SigmaExp":0.51, "merged":False}
X["ggH_1J_PTH_120_200_hzz"] = {"bestfit":1.52, "Up01Sigma":1.09, "Down01Sigma":0.96, "Up01SigmaExp":1.09, "Down01SigmaExp":0.96, "merged":False}
X["ggH_GE2J_MJJ_0_350_PTH_0_60_hzz"] = {"bestfit":1.47, "Up01Sigma":1.35, "Down01Sigma":1.13, "Up01SigmaExp":1.35, "Down01SigmaExp":1.13, "merged":False}
X["ggH_GE2J_MJJ_0_350_PTH_60_120_hzz"] = {"bestfit":1.59, "Up01Sigma":0.83, "Down01Sigma":0.8, "Up01SigmaExp":0.83, "Down01SigmaExp":0.8, "merged":False}
X["ggH_GE2J_MJJ_0_350_PTH_120_200_hzz"] = {"bestfit":1.16, "Up01Sigma":0.87, "Down01Sigma":0.75, "Up01SigmaExp":0.87, "Down01SigmaExp":0.75, "merged":False}
X["ggH_PTH_GT200_hzz"] = {"bestfit":0.47, "Up01Sigma":0.51, "Down01Sigma":0.47, "Up01SigmaExp":0.51, "Down01SigmaExp":0.47, "merged":False}
X["ggH_VBFlike_hzz"] = {"bestfit":0.0, "Up01Sigma":3.28, "Down01Sigma":3.28, "Up01SigmaExp":3.28, "Down01SigmaExp":3.28, "merged":True, "STXS_fractions":{"ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25":0.315,"ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25":0.385,"ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25":0.14,"ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25":0.16}}
# qqH bins
X["qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz"] = {"bestfit":1.71, "Up01Sigma":1.91, "Down01Sigma":1.71, "Up01SigmaExp":1.91, "Down01SigmaExp":1.71, "merged":False}
X["qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz"] = {"bestfit":0.93, "Up01Sigma":1.17, "Down01Sigma":0.9, "Up01SigmaExp":1.17, "Down01SigmaExp":0.9, "merged":False}
X["qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz"] = {"bestfit":2.89, "Up01Sigma":2.88, "Down01Sigma":2.89, "Up01SigmaExp":2.88, "Down01SigmaExp":2.89, "merged":True, "STXS_fractions":{"qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25":0.49,"qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25":0.51}}
#X["qqH_REST_hzz"] = {"bestfit":0.0, "Up01Sigma":2.43, "Down01Sigma":2.43, "Up01SigmaExp":2.43, "Down01SigmaExp":2.43, "merged":False}
X["qqH_GE2J_MJJ_GT350_PTH_GT200_hzz"] = {"bestfit":0.0, "Up01Sigma":0.73, "Down01Sigma":0.73, "Up01SigmaExp":0.73, "Down01SigmaExp":0.73, "merged":False}
# VH lep bins
X["VH_lep_PTV_0_150_hzz"] = {"bestfit":3.21, "Up01Sigma":2.49, "Down01Sigma":1.85, "Up01SigmaExp":2.49, "Down01SigmaExp":1.85, "merged":True, "STXS_fractions":{"WH_lep_PTV_0_75":0.405,"WH_lep_PTV_75_150":0.255,"ZH_lep_PTV_0_75":0.203,"ZH_lep_PTV_75_150":0.137}}
X["VH_lep_PTV_GT150_hzz"] = {"bestfit":0.00, "Up01Sigma":1.57, "Down01Sigma":1.57, "Up01SigmaExp":1.57, "Down01SigmaExp":1.57, "merged":True, "STXS_fractions":{"WH_lep_PTV_150_250_0J":0.278,"WH_lep_PTV_150_250_GE1J":0.214,"WH_lep_PTV_GT250":0.160,"ZH_lep_PTV_150_250_0J":0.144,"ZH_lep_PTV_150_250_GE1J":0.118,"ZH_lep_PTV_GT250":0.086}}
# VH had
X["qqH_GE2J_MJJ_60_120_hzz"] = {"bestfit":0.57, "Up01Sigma":1.20, "Down01Sigma":1.20, "Up01SigmaExp":1.20, "Down01SigmaExp":1.20, "merged":False}
# ttH bins
X["ttH_hzz"] = {"bestfit":0.07, "Up01Sigma":0.9, "Down01Sigma":0.9, "Up01SigmaExp":0.9, "Down01SigmaExp":0.9, "merged":False}

# Correlations
rho = od()
rho[("ggH_0J_PTH_0_10_hzz","ggH_0J_PTH_GT10_hzz")] = 0.12
rho[("ggH_0J_PTH_0_10_hzz","ggH_1J_PTH_0_60_hzz")] = -0.01
rho[("ggH_0J_PTH_0_10_hzz","ggH_1J_PTH_60_120_hzz")] = 0.05
rho[("ggH_0J_PTH_0_10_hzz","ggH_1J_PTH_120_200_hzz")] = 0.04
rho[("ggH_0J_PTH_0_10_hzz","ggH_GE2J_MJJ_0_350_PTH_0_60_hzz")] = 0.02
rho[("ggH_0J_PTH_0_10_hzz","ggH_GE2J_MJJ_0_350_PTH_60_120_hzz")] = 0.04
rho[("ggH_0J_PTH_0_10_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = 0.04
rho[("ggH_0J_PTH_0_10_hzz","ggH_PTH_GT200_hzz")] = 0.02
rho[("ggH_0J_PTH_0_10_hzz","ggH_VBFlike_hzz")] = 0.00
rho[("ggH_0J_PTH_0_10_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.05
rho[("ggH_0J_PTH_0_10_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.02
rho[("ggH_0J_PTH_0_10_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.06
rho[("ggH_0J_PTH_0_10_hzz","qqH_REST_hzz")] = 0.
rho[("ggH_0J_PTH_0_10_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("ggH_0J_PTH_0_10_hzz","VH_lep_PTV_0_150_hzz")] = 0.13
rho[("ggH_0J_PTH_0_10_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("ggH_0J_PTH_0_10_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.03
rho[("ggH_0J_PTH_0_10_hzz","ttH_hzz")] = 0.01

rho[("ggH_0J_PTH_GT10_hzz","ggH_1J_PTH_0_60_hzz")] = -0.5
rho[("ggH_0J_PTH_GT10_hzz","ggH_1J_PTH_60_120_hzz")] = 0.14
rho[("ggH_0J_PTH_GT10_hzz","ggH_1J_PTH_120_200_hzz")] = 0.19
rho[("ggH_0J_PTH_GT10_hzz","ggH_GE2J_MJJ_0_350_PTH_0_60_hzz")] = 0.15
rho[("ggH_0J_PTH_GT10_hzz","ggH_GE2J_MJJ_0_350_PTH_60_120_hzz")] = 0.
rho[("ggH_0J_PTH_GT10_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = 0.04
rho[("ggH_0J_PTH_GT10_hzz","ggH_PTH_GT200_hzz")] = 0.06
rho[("ggH_0J_PTH_GT10_hzz","ggH_VBFlike_hzz")] = 0.0
rho[("ggH_0J_PTH_GT10_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.14
rho[("ggH_0J_PTH_GT10_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.04
rho[("ggH_0J_PTH_GT10_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = -0.03
rho[("ggH_0J_PTH_GT10_hzz","qqH_REST_hzz")] = -0.06
rho[("ggH_0J_PTH_GT10_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.0
rho[("ggH_0J_PTH_GT10_hzz","VH_lep_PTV_0_150_hzz")] = -0.37
rho[("ggH_0J_PTH_GT10_hzz","VH_lep_PTV_GT150_hzz")] = 0.0
rho[("ggH_0J_PTH_GT10_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.11
rho[("ggH_0J_PTH_GT10_hzz","ttH_hzz")] = 0.0

rho[("ggH_1J_PTH_0_60_hzz","ggH_1J_PTH_60_120_hzz")] = -0.02
rho[("ggH_1J_PTH_0_60_hzz","ggH_1J_PTH_120_200_hzz")] = 0.08
rho[("ggH_1J_PTH_0_60_hzz","ggH_GE2J_MJJ_0_350_PTH_0_60_hzz")] = -0.31
rho[("ggH_1J_PTH_0_60_hzz","ggH_GE2J_MJJ_0_350_PTH_60_120_hzz")] = 0.24
rho[("ggH_1J_PTH_0_60_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = 0.14
rho[("ggH_1J_PTH_0_60_hzz","ggH_PTH_GT200_hzz")] = 0.15
rho[("ggH_1J_PTH_0_60_hzz","ggH_VBFlike_hzz")] = 0.
rho[("ggH_1J_PTH_0_60_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.06
rho[("ggH_1J_PTH_0_60_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.05
rho[("ggH_1J_PTH_0_60_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.1
rho[("ggH_1J_PTH_0_60_hzz","qqH_REST_hzz")] = -0.46
rho[("ggH_1J_PTH_0_60_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("ggH_1J_PTH_0_60_hzz","VH_lep_PTV_0_150_hzz")] = -0.04
rho[("ggH_1J_PTH_0_60_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("ggH_1J_PTH_0_60_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.05 
rho[("ggH_1J_PTH_0_60_hzz","ttH_hzz")] = 0.02

rho[("ggH_1J_PTH_60_120_hzz","ggH_1J_PTH_120_200_hzz")] = 0.15
rho[("ggH_1J_PTH_60_120_hzz","ggH_GE2J_MJJ_0_350_PTH_0_60_hzz")] = 0.23
rho[("ggH_1J_PTH_60_120_hzz","ggH_GE2J_MJJ_0_350_PTH_60_120_hzz")] = -0.21
rho[("ggH_1J_PTH_60_120_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = 0.22
rho[("ggH_1J_PTH_60_120_hzz","ggH_PTH_GT200_hzz")] = 0.24
rho[("ggH_1J_PTH_60_120_hzz","ggH_VBFlike_hzz")] = 0.0
rho[("ggH_1J_PTH_60_120_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.11
rho[("ggH_1J_PTH_60_120_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.08
rho[("ggH_1J_PTH_60_120_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = -0.02
rho[("ggH_1J_PTH_60_120_hzz","qqH_REST_hzz")] = -0.8
rho[("ggH_1J_PTH_60_120_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("ggH_1J_PTH_60_120_hzz","VH_lep_PTV_0_150_hzz")] = -0.11
rho[("ggH_1J_PTH_60_120_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("ggH_1J_PTH_60_120_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.01 
rho[("ggH_1J_PTH_60_120_hzz","ttH_hzz")] = 0.04

rho[("ggH_1J_PTH_120_200_hzz","ggH_GE2J_MJJ_0_350_PTH_0_60_hzz")] = 0.07
rho[("ggH_1J_PTH_120_200_hzz","ggH_GE2J_MJJ_0_350_PTH_60_120_hzz")] = 0.13
rho[("ggH_1J_PTH_120_200_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = -0.12
rho[("ggH_1J_PTH_120_200_hzz","ggH_PTH_GT200_hzz")] = 0.26
rho[("ggH_1J_PTH_120_200_hzz","ggH_VBFlike_hzz")] = -0.01
rho[("ggH_1J_PTH_120_200_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.07
rho[("ggH_1J_PTH_120_200_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.04
rho[("ggH_1J_PTH_120_200_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = -0.01
rho[("ggH_1J_PTH_120_200_hzz","qqH_REST_hzz")] = -0.35
rho[("ggH_1J_PTH_120_200_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("ggH_1J_PTH_120_200_hzz","VH_lep_PTV_0_150_hzz")] = -0.02
rho[("ggH_1J_PTH_120_200_hzz","VH_lep_PTV_GT150_hzz")] = 0.0
rho[("ggH_1J_PTH_120_200_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.03
rho[("ggH_1J_PTH_120_200_hzz","ttH_hzz")] = 0.02

rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","ggH_GE2J_MJJ_0_350_PTH_60_120_hzz")] = 0.11
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = 0.14
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","ggH_PTH_GT200_hzz")] = 0.19
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","ggH_VBFlike_hzz")] = -0.01
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.04
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.05
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.03
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","qqH_REST_hzz")] = -0.17
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.0
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","VH_lep_PTV_0_150_hzz")] = -0.04
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","VH_lep_PTV_GT150_hzz")] = 0.0
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.08
rho[("ggH_GE2J_MJJ_0_350_PTH_0_60_hzz","ttH_hzz")] = 0.

rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","ggH_GE2J_MJJ_0_350_PTH_120_200_hzz")] = 0.18
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","ggH_PTH_GT200_hzz")] = 0.25
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","ggH_VBFlike_hzz")] = -0.01
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.07
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.07
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = -0.01
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","qqH_REST_hzz")] = -0.28
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","VH_lep_PTV_0_150_hzz")] = -0.03
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.17 
rho[("ggH_GE2J_MJJ_0_350_PTH_60_120_hzz","ttH_hzz")] = -0.03

rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","ggH_PTH_GT200_hzz")] = 0.27
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","ggH_VBFlike_hzz")] = -0.01
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.09
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.08
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.02
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","qqH_REST_hzz")] = -0.38
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.0
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","VH_lep_PTV_0_150_hzz")] = 0.01
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.21 
rho[("ggH_GE2J_MJJ_0_350_PTH_120_200_hzz","ttH_hzz")] = -0.02

rho[("ggH_PTH_GT200_hzz","ggH_VBFlike_hzz")] = 0.0
rho[("ggH_PTH_GT200_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.19
rho[("ggH_PTH_GT200_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.11
rho[("ggH_PTH_GT200_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.09
rho[("ggH_PTH_GT200_hzz","qqH_REST_hzz")] = -0.83
rho[("ggH_PTH_GT200_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = -0.06
rho[("ggH_PTH_GT200_hzz","VH_lep_PTV_0_150_hzz")] = 0.
rho[("ggH_PTH_GT200_hzz","VH_lep_PTV_GT150_hzz")] = -0.02
rho[("ggH_PTH_GT200_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.16
rho[("ggH_PTH_GT200_hzz","ttH_hzz")] = -0.04

rho[("ggH_VBFlike_hzz","qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz")] = -0.51
rho[("ggH_VBFlike_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = -0.39
rho[("ggH_VBFlike_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = -0.68
rho[("ggH_VBFlike_hzz","qqH_REST_hzz")] = 0.
rho[("ggH_VBFlike_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("ggH_VBFlike_hzz","VH_lep_PTV_0_150_hzz")] = 0.
rho[("ggH_VBFlike_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("ggH_VBFlike_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0. 
rho[("ggH_VBFlike_hzz","ttH_hzz")] = 0.03

rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz")] = 0.19
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.35
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","qqH_REST_hzz")] = -0.16
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","VH_lep_PTV_0_150_hzz")] = 0.02
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0. 
rho[("qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_hzz","ttH_hzz")] = 0.

rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz")] = 0.04
rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","qqH_REST_hzz")] = -0.19
rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","VH_lep_PTV_0_150_hzz")] = 0.01
rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","VH_lep_PTV_GT150_hzz")] = 0.00
rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.01
rho[("qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_hzz","ttH_hzz")] = 0.03

rho[("qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz","qqH_REST_hzz")] = -0.04
rho[("qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.
rho[("qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz","VH_lep_PTV_0_150_hzz")] = 0.03
rho[("qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.03 
rho[("qqH_GE2J_MJJ_GT300_PTH_0_200_PTHJJ_GT25_hzz","ttH_hzz")] = -0.2

rho[("qqH_REST_hzz","qqH_GE2J_MJJ_GT350_PTH_GT200_hzz")] = 0.00
rho[("qqH_REST_hzz","VH_lep_PTV_0_150_hzz")] = 0.02
rho[("qqH_REST_hzz","VH_lep_PTV_GT150_hzz")] = 0.00
rho[("qqH_REST_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.03
rho[("qqH_REST_hzz","ttH_hzz")] = -0.01

rho[("qqH_GE2J_MJJ_GT350_PTH_GT200_hzz","VH_lep_PTV_0_150_hzz")] = 0.01
rho[("qqH_GE2J_MJJ_GT350_PTH_GT200_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("qqH_GE2J_MJJ_GT350_PTH_GT200_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.01
rho[("qqH_GE2J_MJJ_GT350_PTH_GT200_hzz","ttH_hzz")] = 0.01

rho[("VH_lep_PTV_0_150_hzz","VH_lep_PTV_GT150_hzz")] = 0.
rho[("VH_lep_PTV_0_150_hzz","qqH_GE2J_MJJ_60_120_hzz")] = -0.02 
rho[("VH_lep_PTV_0_150_hzz","ttH_hzz")] = -0.04

rho[("VH_lep_PTV_GT150_hzz","qqH_GE2J_MJJ_60_120_hzz")] = 0.
rho[("VH_lep_PTV_GT150_hzz","ttH_hzz")] = 0.

rho[("qqH_GE2J_MJJ_60_120_hzz","ttH_hzz")] = -0.06


