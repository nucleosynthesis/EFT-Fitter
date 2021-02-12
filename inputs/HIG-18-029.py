from collections import OrderedDict as od

name = "HIG-18-029"

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
X["ggH_0J_hgg"] = {"bestfit":1.17, "Up01Sigma":0.2, "Down01Sigma":0.2, "Up01SigmaExp":0.2, "Down01SigmaExp":0.18, "merged":False}
X["ggH_1J_PTH_0_60_hgg"] = {"bestfit":1.5, "Up01Sigma":0.7, "Down01Sigma":0.6, "Up01SigmaExp":0.8, "Down01SigmaExp":0.6, "merged":False}
X["ggH_1J_PTH_60_120_hgg"] = {"bestfit":0.5, "Up01Sigma":0.5, "Down01Sigma":0.4, "Up01SigmaExp":0.6, "Down01SigmaExp":0.6, "merged":False}
X["ggH_1J_PTH_120_200_hgg"] = {"bestfit":2.0, "Up01Sigma":1.0, "Down01Sigma":0.7, "Up01SigmaExp":0.9, "Down01SigmaExp":0.8, "merged":False}
X["ggH_1J_PTH_GT200_hgg"] = {"bestfit":1.8, "Up01Sigma":1.7, "Down01Sigma":1.5, "Up01SigmaExp":1.6, "Down01SigmaExp":1.5, "merged":False}
X["ggH_GE2J_PTH_0_60_hgg"] = {"bestfit":0.3, "Up01Sigma":1.5, "Down01Sigma":1.3, "Up01SigmaExp":1.4, "Down01SigmaExp":1.3, "merged":False}
X["ggH_GE2J_PTH_60_120_hgg"] = {"bestfit":2.6, "Up01Sigma":1.1, "Down01Sigma":1.1, "Up01SigmaExp":1.0, "Down01SigmaExp":0.9, "merged":False}
X["ggH_GE2J_PTH_120_200_hgg"] = {"bestfit":0.6, "Up01Sigma":0.8, "Down01Sigma":0.8, "Up01SigmaExp":0.9, "Down01SigmaExp":0.8, "merged":False}
X["ggH_GE2J_PTH_GT200_hgg"] = {"bestfit":2.8, "Up01Sigma":1.1, "Down01Sigma":1.2, "Up01SigmaExp":1.3, "Down01SigmaExp":1.2, "merged":False}
X["ggH_VBFlike_hgg"] = {"bestfit":0.0, "Up01Sigma":2.8, "Down01Sigma":2.5, "Up01SigmaExp":2.8, "Down01SigmaExp":2.5, "merged":True, "STXS_fractions":{"ggH_VBFTOPO_JET3VETO":0.455,"ggH_VBFTOPO_JET3":0.545}}
X["qqH_VBFTOPO_JET3VETO_hgg"] = {"bestfit":1.3, "Up01Sigma":0.6, "Down01Sigma":0.5, "Up01SigmaExp":0.9, "Down01SigmaExp":0.8, "merged":False}
X["qqH_VBFTOPO_JET3_hgg"] = {"bestfit":0.0, "Up01Sigma":2.2, "Down01Sigma":2.2, "Up01SigmaExp":2.2, "Down01SigmaExp":2.2, "merged":False}
X["qqH_other_hgg"] = {"bestfit":0.0, "Up01Sigma":1.8, "Down01Sigma":2.3, "Up01SigmaExp":2.4, "Down01SigmaExp":2.3, "merged":True, "STXS_fractions":{"qqH_PTJET1_GT200":0.04,"qqH_VH2JET":0.14,"qqH_REST":0.82}}

# Correlations
rho = od()
rho[("ggH_0J_hgg","ggH_1J_PTH_0_60_hgg")] = -0.24 # Correlation coefficients
rho[("ggH_0J_hgg","ggH_1J_PTH_60_120_hgg")] = 0.03
rho[("ggH_0J_hgg","ggH_1J_PTH_120_200_hgg")] = 0.16
rho[("ggH_0J_hgg","ggH_1J_PTH_GT200_hgg")] = 0.11
rho[("ggH_0J_hgg","ggH_GE2J_PTH_0_60_hgg")] = -0.07
rho[("ggH_0J_hgg","ggH_GE2J_PTH_60_120_hgg")] = 0.0
rho[("ggH_0J_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.0
rho[("ggH_0J_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.13
rho[("ggH_0J_hgg","ggH_VBFlike_hgg")] = 0.29
rho[("ggH_0J_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = -0.22
rho[("ggH_0J_hgg","qqH_VBFTOPO_JET3_hgg")] = -0.26
rho[("ggH_0J_hgg","qqH_other_hgg")] = -0.11

rho[("ggH_1J_PTH_0_60_hgg","ggH_1J_PTH_60_120_hgg")] = 0.29
rho[("ggH_1J_PTH_0_60_hgg","ggH_1J_PTH_120_200_hgg")] = 0.26
rho[("ggH_1J_PTH_0_60_hgg","ggH_1J_PTH_GT200_hgg")] = 0.21
rho[("ggH_1J_PTH_0_60_hgg","ggH_GE2J_PTH_0_60_hgg")] = 0.03
rho[("ggH_1J_PTH_0_60_hgg","ggH_GE2J_PTH_60_120_hgg")] = 0.33
rho[("ggH_1J_PTH_0_60_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.36
rho[("ggH_1J_PTH_0_60_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.34
rho[("ggH_1J_PTH_0_60_hgg","ggH_VBFlike_hgg")] = -0.10
rho[("ggH_1J_PTH_0_60_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.12
rho[("ggH_1J_PTH_0_60_hgg","qqH_VBFTOPO_JET3_hgg")] = 0.12
rho[("ggH_1J_PTH_0_60_hgg","qqH_other_hgg")] = -0.41

rho[("ggH_1J_PTH_60_120_hgg","ggH_1J_PTH_120_200_hgg")] = 0.43
rho[("ggH_1J_PTH_60_120_hgg","ggH_1J_PTH_GT200_hgg")] = 0.37
rho[("ggH_1J_PTH_60_120_hgg","ggH_GE2J_PTH_0_60_hgg")] = 0.39
rho[("ggH_1J_PTH_60_120_hgg","ggH_GE2J_PTH_60_120_hgg")] = 0.3
rho[("ggH_1J_PTH_60_120_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.56
rho[("ggH_1J_PTH_60_120_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.56
rho[("ggH_1J_PTH_60_120_hgg","ggH_VBFlike_hgg")] = -0.1
rho[("ggH_1J_PTH_60_120_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.16
rho[("ggH_1J_PTH_60_120_hgg","qqH_VBFTOPO_JET3_hgg")] = 0.12
rho[("ggH_1J_PTH_60_120_hgg","qqH_other_hgg")] = -0.68

rho[("ggH_1J_PTH_120_200_hgg","ggH_1J_PTH_GT200_hgg")] = 0.4
rho[("ggH_1J_PTH_120_200_hgg","ggH_GE2J_PTH_0_60_hgg")] = 0.21
rho[("ggH_1J_PTH_120_200_hgg","ggH_GE2J_PTH_60_120_hgg")] = 0.33
rho[("ggH_1J_PTH_120_200_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.33
rho[("ggH_1J_PTH_120_200_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.57
rho[("ggH_1J_PTH_120_200_hgg","ggH_VBFlike_hgg")] = 0.21
rho[("ggH_1J_PTH_120_200_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = -0.1
rho[("ggH_1J_PTH_120_200_hgg","qqH_VBFTOPO_JET3_hgg")] = -0.18
rho[("ggH_1J_PTH_120_200_hgg","qqH_other_hgg")] = -0.64

rho[("ggH_1J_PTH_GT200_hgg","ggH_GE2J_PTH_0_60_hgg")] = 0.22
rho[("ggH_1J_PTH_GT200_hgg","ggH_GE2J_PTH_60_120_hgg")] = 0.31
rho[("ggH_1J_PTH_GT200_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.4
rho[("ggH_1J_PTH_GT200_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.36
rho[("ggH_1J_PTH_GT200_hgg","ggH_VBFlike_hgg")] = 0.08
rho[("ggH_1J_PTH_GT200_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.0
rho[("ggH_1J_PTH_GT200_hgg","qqH_VBFTOPO_JET3_hgg")] = -0.06
rho[("ggH_1J_PTH_GT200_hgg","qqH_other_hgg")] = -0.54

rho[("ggH_GE2J_PTH_0_60_hgg","ggH_GE2J_PTH_60_120_hgg")] = 0.47
rho[("ggH_GE2J_PTH_0_60_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.53
rho[("ggH_GE2J_PTH_0_60_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.37
rho[("ggH_GE2J_PTH_0_60_hgg","ggH_VBFlike_hgg")] = -0.52
rho[("ggH_GE2J_PTH_0_60_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.48
rho[("ggH_GE2J_PTH_0_60_hgg","qqH_VBFTOPO_JET3_hgg")] = 0.51
rho[("ggH_GE2J_PTH_0_60_hgg","qqH_other_hgg")] = -0.48

rho[("ggH_GE2J_PTH_60_120_hgg","ggH_GE2J_PTH_120_200_hgg")] = 0.64
rho[("ggH_GE2J_PTH_60_120_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.58
rho[("ggH_GE2J_PTH_60_120_hgg","ggH_VBFlike_hgg")] = -0.34
rho[("ggH_GE2J_PTH_60_120_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.35
rho[("ggH_GE2J_PTH_60_120_hgg","qqH_VBFTOPO_JET3_hgg")] = 0.35
rho[("ggH_GE2J_PTH_60_120_hgg","qqH_other_hgg")] = -0.64

rho[("ggH_GE2J_PTH_120_200_hgg","ggH_GE2J_PTH_GT200_hgg")] = 0.66
rho[("ggH_GE2J_PTH_120_200_hgg","ggH_VBFlike_hgg")] = -0.32
rho[("ggH_GE2J_PTH_120_200_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.36
rho[("ggH_GE2J_PTH_120_200_hgg","qqH_VBFTOPO_JET3_hgg")] = 0.33
rho[("ggH_GE2J_PTH_120_200_hgg","qqH_other_hgg")] = -0.77

rho[("ggH_GE2J_PTH_GT200_hgg","ggH_VBFlike_hgg")] = 0.03
rho[("ggH_GE2J_PTH_GT200_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = 0.07
rho[("ggH_GE2J_PTH_GT200_hgg","qqH_VBFTOPO_JET3_hgg")] = -0.01
rho[("ggH_GE2J_PTH_GT200_hgg","qqH_other_hgg")] = -0.82

rho[("ggH_VBFlike_hgg","qqH_VBFTOPO_JET3VETO_hgg")] = -0.83
rho[("ggH_VBFlike_hgg","qqH_VBFTOPO_JET3_hgg")] = -0.96
rho[("ggH_VBFlike_hgg","qqH_other_hgg")] = 0.00

rho[("qqH_VBFTOPO_JET3VETO_hgg","qqH_VBFTOPO_JET3_hgg")] = 0.74
rho[("qqH_VBFTOPO_JET3VETO_hgg","qqH_other_hgg")] = -0.14

rho[("qqH_VBFTOPO_JET3_hgg","qqH_other_hgg")] = -0.03
