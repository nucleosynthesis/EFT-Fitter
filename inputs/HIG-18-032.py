from collections import OrderedDict as od

# H->tautau (2016+2017): ggH + qqH STXS stage 1.0

name = "HIG-18-032"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["ggH_0J_htt"] = {"bestfit":-0.4, "Up01Sigma":0.76, "Down01Sigma":0.75, "Up01SigmaExp":0.76, "Down01SigmaExp":0.75, "merged":False}
X["ggH_1J_PTH_0_60_htt"] = {"bestfit":-0.34, "Up01Sigma":1.37, "Down01Sigma":1.39, "Up01SigmaExp":1.37, "Down01SigmaExp":1.39, "merged":False}
X["ggH_1J_PTH_60_120_htt"] = {"bestfit":1.26, "Up01Sigma":1.56, "Down01Sigma":1.50, "Up01SigmaExp":1.56, "Down01SigmaExp":1.50, "merged":False}
X["ggH_1J_PTH_GT120_htt"] = {"bestfit":1.80, "Up01Sigma":1.18, "Down01Sigma":1.01, "Up01SigmaExp":1.18, "Down01SigmaExp":1.01, "merged":True, "STXS_fractions":{"ggH_1J_PTH_120_200":0.8095,"ggH_1J_PTH_GT200":0.1905}}
X["ggH_GE2J_htt"] = {"bestfit":0.47, "Up01Sigma":0.91, "Down01Sigma":0.86, "Up01SigmaExp":0.91, "Down01SigmaExp":0.86, "merged":True, "STXS_fractions":{"ggH_GE2J_PTH_0_60":0.238,"ggH_GE2J_PTH_60_120":0.369,"ggH_GE2J_PTH_120_200":0.189,"ggH_GE2J_PTH_GT200":0.082,"ggH_VBFTOPO_JET3VETO":0.049,"ggH_VBFTOPO_JET3":0.074}}
# qqH Bins
X["qqH_VBFTOPO_htt"] = {"bestfit":1.0, "Up01Sigma":0.30, "Down01Sigma":0.29, "Up01SigmaExp":0.30, "Down01SigmaExp":0.29, "merged":True, "STXS_fractions":{"qqH_VBFTOPO_JET3VETO":0.728,"qqH_VBFTOPO_JET3":0.272}}
X["qqH_VH2JET_htt"] = {"bestfit":-1.17, "Up01Sigma":1.54, "Down01Sigma":1.47, "Up01SigmaExp":1.54, "Down01SigmaExp":1.47, "merged":False}
X["qqH_PTJET1_GT200_htt"] = {"bestfit":1.41, "Up01Sigma":1.03, "Down01Sigma":1.05, "Up01SigmaExp":1.03, "Down01SigmaExp":1.05, "merged":False}
X["qqH_REST_htt"] = {"bestfit":-1.06, "Up01Sigma":2.75, "Down01Sigma":2.67, "Up01SigmaExp":2.75, "Down01SigmaExp":2.67, "merged":False}

# Correlations
rho = od()
rho[("ggH_0J_htt","ggH_1J_PTH_0_60_htt")] = 0.14
rho[("ggH_0J_htt","ggH_1J_PTH_60_120_htt")] = 0.16
rho[("ggH_0J_htt","ggH_1J_PTH_GT120_htt")] = 0.06
rho[("ggH_0J_htt","ggH_GE2J_htt")] = 0.14
rho[("ggH_0J_htt","qqH_VBFTOPO_htt")] = 0.02
rho[("ggH_0J_htt","qqH_VH2JET_htt")] = -0.05
rho[("ggH_0J_htt","qqH_PTJET1_GT200_htt")] = -0.04 
rho[("ggH_0J_htt","qqH_REST_htt")] = -0.06

rho[("ggH_1J_PTH_0_60_htt","ggH_1J_PTH_60_120_htt")] = -0.28
rho[("ggH_1J_PTH_0_60_htt","ggH_1J_PTH_GT120_htt")] = 0.1
rho[("ggH_1J_PTH_0_60_htt","ggH_GE2J_htt")] = 0.08
rho[("ggH_1J_PTH_0_60_htt","qqH_VBFTOPO_htt")] = 0.02
rho[("ggH_1J_PTH_0_60_htt","qqH_VH2JET_htt")] = -0.03
rho[("ggH_1J_PTH_0_60_htt","qqH_PTJET1_GT200_htt")] = -0.06
rho[("ggH_1J_PTH_0_60_htt","qqH_REST_htt")] = -0.05

rho[("ggH_1J_PTH_60_120_htt","ggH_1J_PTH_GT120_htt")] = -0.05
rho[("ggH_1J_PTH_60_120_htt","ggH_GE2J_htt")] = 0.26
rho[("ggH_1J_PTH_60_120_htt","qqH_VBFTOPO_htt")] = 0.21
rho[("ggH_1J_PTH_60_120_htt","qqH_VH2JET_htt")] = -0.11
rho[("ggH_1J_PTH_60_120_htt","qqH_PTJET1_GT200_htt")] = -0.09
rho[("ggH_1J_PTH_60_120_htt","qqH_REST_htt")] = -0.5

rho[("ggH_1J_PTH_GT120_htt","ggH_GE2J_htt")] = -0.09
rho[("ggH_1J_PTH_GT120_htt","qqH_VBFTOPO_htt")] = 0.11
rho[("ggH_1J_PTH_GT120_htt","qqH_VH2JET_htt")] = 0.1
rho[("ggH_1J_PTH_GT120_htt","qqH_PTJET1_GT200_htt")] = -0.16
rho[("ggH_1J_PTH_GT120_htt","qqH_REST_htt")] = -0.18

rho[("ggH_GE2J_htt","qqH_VBFTOPO_htt")] = 0.3
rho[("ggH_GE2J_htt","qqH_VH2JET_htt")] = -0.61
rho[("ggH_GE2J_htt","qqH_PTJET1_GT200_htt")] = -0.48
rho[("ggH_GE2J_htt","qqH_REST_htt")] = -0.67

rho[("qqH_VBFTOPO_htt","qqH_VH2JET_htt")] = -0.16
rho[("qqH_VBFTOPO_htt","qqH_PTJET1_GT200_htt")] = -0.23 
rho[("qqH_VBFTOPO_htt","qqH_REST_htt")] = -0.45

rho[("qqH_VH2JET_htt","qqH_PTJET1_GT200_htt")] = 0.29
rho[("qqH_VH2JET_htt","qqH_REST_htt")] = 0.38

rho[("qqH_PTJET1_GT200_htt","qqH_REST_htt")] = 0.4
