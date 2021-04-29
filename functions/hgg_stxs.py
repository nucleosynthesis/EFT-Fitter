from collections import OrderedDict as od
import sys
import json

# Functions inferred from the purity matrix of the HIG-19-015 anaylsis: saved in json file
#with open("functions/hgg_stxs/stage1p2_extended_pruned.json","r") as jf: fracs = json.load(jf)
with open("functions/hgg_stxs/stage1p2_extended.json","r") as jf: fracs = json.load(jf)

# Define function to calc
def STXS_func(pois,cat="0J_PTH_0_10_Tag0"):
  eq = 0
  for p,f in iter(fracs[cat].items()):
    eq += pois[p]*f
  return eq

def RECO_0J_PTH_0_10_Tag0(pois):
  return STXS_func(pois,cat="0J_PTH_0_10_Tag0")

def RECO_0J_PTH_0_10_Tag0(pois):
  return STXS_func(pois,cat="0J_PTH_0_10_Tag0")

def RECO_0J_PTH_0_10_Tag1(pois):
  return STXS_func(pois,cat="0J_PTH_0_10_Tag1")

def RECO_0J_PTH_0_10_Tag2(pois):
  return STXS_func(pois,cat="0J_PTH_0_10_Tag2")

def RECO_0J_PTH_GT10_Tag0(pois):
  return STXS_func(pois,cat="0J_PTH_GT10_Tag0")

def RECO_0J_PTH_GT10_Tag1(pois):
  return STXS_func(pois,cat="0J_PTH_GT10_Tag1")

def RECO_0J_PTH_GT10_Tag2(pois):
  return STXS_func(pois,cat="0J_PTH_GT10_Tag2")

def RECO_1J_PTH_0_60_Tag0(pois):
  return STXS_func(pois,cat="1J_PTH_0_60_Tag0")

def RECO_1J_PTH_0_60_Tag1(pois):
  return STXS_func(pois,cat="1J_PTH_0_60_Tag1")

def RECO_1J_PTH_0_60_Tag2(pois):
  return STXS_func(pois,cat="1J_PTH_0_60_Tag2")

def RECO_1J_PTH_120_200_Tag0(pois):
  return STXS_func(pois,cat="1J_PTH_120_200_Tag0")

def RECO_1J_PTH_120_200_Tag1(pois):
  return STXS_func(pois,cat="1J_PTH_120_200_Tag1")

def RECO_1J_PTH_120_200_Tag2(pois):
  return STXS_func(pois,cat="1J_PTH_120_200_Tag2")

def RECO_1J_PTH_60_120_Tag0(pois):
  return STXS_func(pois,cat="1J_PTH_60_120_Tag0")

def RECO_1J_PTH_60_120_Tag1(pois):
  return STXS_func(pois,cat="1J_PTH_60_120_Tag1")

def RECO_1J_PTH_60_120_Tag2(pois):
  return STXS_func(pois,cat="1J_PTH_60_120_Tag2")

def RECO_GE2J_PTH_0_60_Tag0(pois):
  return STXS_func(pois,cat="GE2J_PTH_0_60_Tag0")

def RECO_GE2J_PTH_0_60_Tag1(pois):
  return STXS_func(pois,cat="GE2J_PTH_0_60_Tag1")

def RECO_GE2J_PTH_0_60_Tag2(pois):
  return STXS_func(pois,cat="GE2J_PTH_0_60_Tag2")

def RECO_GE2J_PTH_120_200_Tag0(pois):
  return STXS_func(pois,cat="GE2J_PTH_120_200_Tag0")

def RECO_GE2J_PTH_120_200_Tag1(pois):
  return STXS_func(pois,cat="GE2J_PTH_120_200_Tag1")

def RECO_GE2J_PTH_120_200_Tag2(pois):
  return STXS_func(pois,cat="GE2J_PTH_120_200_Tag2")

def RECO_GE2J_PTH_60_120_Tag0(pois):
  return STXS_func(pois,cat="GE2J_PTH_60_120_Tag0")

def RECO_GE2J_PTH_60_120_Tag1(pois):
  return STXS_func(pois,cat="GE2J_PTH_60_120_Tag1")

def RECO_GE2J_PTH_60_120_Tag2(pois):
  return STXS_func(pois,cat="GE2J_PTH_60_120_Tag2")

def RECO_PTH_200_300_Tag0(pois):
  return STXS_func(pois,cat="PTH_200_300_Tag0")

def RECO_PTH_200_300_Tag1(pois):
  return STXS_func(pois,cat="PTH_200_300_Tag1")

def RECO_PTH_300_450_Tag0(pois):
  return STXS_func(pois,cat="PTH_300_450_Tag0")

def RECO_PTH_300_450_Tag1(pois):
  return STXS_func(pois,cat="PTH_300_450_Tag1")

def RECO_PTH_450_650_Tag0(pois):
  return STXS_func(pois,cat="PTH_450_650_Tag0")

def RECO_PTH_GT650_Tag0(pois):
  return STXS_func(pois,cat="PTH_GT650_Tag0")

def RECO_THQ_LEP(pois):
  return STXS_func(pois,cat="THQ_LEP")

def RECO_TTH_HAD_PTH_0_60_Tag0(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_0_60_Tag0")

def RECO_TTH_HAD_PTH_0_60_Tag1(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_0_60_Tag1")

def RECO_TTH_HAD_PTH_0_60_Tag2(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_0_60_Tag2")

def RECO_TTH_HAD_PTH_120_200_Tag0(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_120_200_Tag0")

def RECO_TTH_HAD_PTH_120_200_Tag1(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_120_200_Tag1")

def RECO_TTH_HAD_PTH_120_200_Tag2(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_120_200_Tag2")

def RECO_TTH_HAD_PTH_120_200_Tag3(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_120_200_Tag3")

def RECO_TTH_HAD_PTH_200_300_Tag0(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_200_300_Tag0")

def RECO_TTH_HAD_PTH_200_300_Tag1(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_200_300_Tag1")

def RECO_TTH_HAD_PTH_200_300_Tag2(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_200_300_Tag2")

def RECO_TTH_HAD_PTH_60_120_Tag0(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_60_120_Tag0")

def RECO_TTH_HAD_PTH_60_120_Tag1(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_60_120_Tag1")

def RECO_TTH_HAD_PTH_60_120_Tag2(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_60_120_Tag2")

def RECO_TTH_HAD_PTH_GT300_Tag0(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_GT300_Tag0")

def RECO_TTH_HAD_PTH_GT300_Tag1(pois):
  return STXS_func(pois,cat="TTH_HAD_PTH_GT300_Tag1")

def RECO_TTH_LEP_PTH_0_60_Tag0(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_0_60_Tag0")

def RECO_TTH_LEP_PTH_0_60_Tag1(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_0_60_Tag1")

def RECO_TTH_LEP_PTH_0_60_Tag2(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_0_60_Tag2")

def RECO_TTH_LEP_PTH_120_200_Tag0(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_120_200_Tag0")

def RECO_TTH_LEP_PTH_120_200_Tag1(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_120_200_Tag1")

def RECO_TTH_LEP_PTH_200_300_Tag0(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_200_300_Tag0")

def RECO_TTH_LEP_PTH_60_120_Tag0(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_60_120_Tag0")

def RECO_TTH_LEP_PTH_60_120_Tag1(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_60_120_Tag1")

def RECO_TTH_LEP_PTH_60_120_Tag2(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_60_120_Tag2")

def RECO_TTH_LEP_PTH_GT300_Tag0(pois):
  return STXS_func(pois,cat="TTH_LEP_PTH_GT300_Tag0")

def RECO_VBFLIKEGGH_Tag0(pois):
  return STXS_func(pois,cat="VBFLIKEGGH_Tag0")

def RECO_VBFLIKEGGH_Tag1(pois):
  return STXS_func(pois,cat="VBFLIKEGGH_Tag1")

def RECO_VBFTOPO_BSM_Tag0(pois):
  return STXS_func(pois,cat="VBFTOPO_BSM_Tag0")

def RECO_VBFTOPO_BSM_Tag1(pois):
  return STXS_func(pois,cat="VBFTOPO_BSM_Tag1")

def RECO_VBFTOPO_JET3VETO_HIGHMJJ_Tag0(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3VETO_HIGHMJJ_Tag0")

def RECO_VBFTOPO_JET3VETO_HIGHMJJ_Tag1(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3VETO_HIGHMJJ_Tag1")

def RECO_VBFTOPO_JET3VETO_LOWMJJ_Tag0(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3VETO_LOWMJJ_Tag0")

def RECO_VBFTOPO_JET3VETO_LOWMJJ_Tag1(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3VETO_LOWMJJ_Tag1")

def RECO_VBFTOPO_JET3_HIGHMJJ_Tag0(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3_HIGHMJJ_Tag0")

def RECO_VBFTOPO_JET3_HIGHMJJ_Tag1(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3_HIGHMJJ_Tag1")

def RECO_VBFTOPO_JET3_LOWMJJ_Tag0(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3_LOWMJJ_Tag0")

def RECO_VBFTOPO_JET3_LOWMJJ_Tag1(pois):
  return STXS_func(pois,cat="VBFTOPO_JET3_LOWMJJ_Tag1")

def RECO_VBFTOPO_VHHAD_Tag0(pois):
  return STXS_func(pois,cat="VBFTOPO_VHHAD_Tag0")

def RECO_VBFTOPO_VHHAD_Tag1(pois):
  return STXS_func(pois,cat="VBFTOPO_VHHAD_Tag1")

def RECO_VH_MET_Tag0(pois):
  return STXS_func(pois,cat="VH_MET_Tag0")

def RECO_VH_MET_Tag1(pois):
  return STXS_func(pois,cat="VH_MET_Tag1")

def RECO_VH_MET_Tag2(pois):
  return STXS_func(pois,cat="VH_MET_Tag2")

def RECO_WH_LEP_PTV_0_75_Tag0(pois):
  return STXS_func(pois,cat="WH_LEP_PTV_0_75_Tag0")

def RECO_WH_LEP_PTV_0_75_Tag1(pois):
  return STXS_func(pois,cat="WH_LEP_PTV_0_75_Tag1")

def RECO_WH_LEP_PTV_75_150_Tag0(pois):
  return STXS_func(pois,cat="WH_LEP_PTV_75_150_Tag0")

def RECO_WH_LEP_PTV_75_150_Tag1(pois):
  return STXS_func(pois,cat="WH_LEP_PTV_75_150_Tag1")

def RECO_WH_LEP_PTV_GT150_Tag0(pois):
  return STXS_func(pois,cat="WH_LEP_PTV_GT150_Tag0")

def RECO_ZH_LEP_Tag0(pois):
  return STXS_func(pois,cat="ZH_LEP_Tag0")

def RECO_ZH_LEP_Tag1(pois):
  return STXS_func(pois,cat="ZH_LEP_Tag1")

functions = od()
functions["0J_PTH_0_10_Tag0"] = RECO_0J_PTH_0_10_Tag0
functions["0J_PTH_0_10_Tag1"] = RECO_0J_PTH_0_10_Tag1
functions["0J_PTH_0_10_Tag2"] = RECO_0J_PTH_0_10_Tag2
functions["0J_PTH_GT10_Tag0"] = RECO_0J_PTH_GT10_Tag0
functions["0J_PTH_GT10_Tag1"] = RECO_0J_PTH_GT10_Tag1
functions["0J_PTH_GT10_Tag2"] = RECO_0J_PTH_GT10_Tag2
functions["1J_PTH_0_60_Tag0"] = RECO_1J_PTH_0_60_Tag0
functions["1J_PTH_0_60_Tag1"] = RECO_1J_PTH_0_60_Tag1
functions["1J_PTH_0_60_Tag2"] = RECO_1J_PTH_0_60_Tag2
functions["1J_PTH_120_200_Tag0"] = RECO_1J_PTH_120_200_Tag0
functions["1J_PTH_120_200_Tag1"] = RECO_1J_PTH_120_200_Tag1
functions["1J_PTH_120_200_Tag2"] = RECO_1J_PTH_120_200_Tag2
functions["1J_PTH_60_120_Tag0"] = RECO_1J_PTH_60_120_Tag0
functions["1J_PTH_60_120_Tag1"] = RECO_1J_PTH_60_120_Tag1
functions["1J_PTH_60_120_Tag2"] = RECO_1J_PTH_60_120_Tag2
functions["GE2J_PTH_0_60_Tag0"] = RECO_GE2J_PTH_0_60_Tag0
functions["GE2J_PTH_0_60_Tag1"] = RECO_GE2J_PTH_0_60_Tag1
functions["GE2J_PTH_0_60_Tag2"] = RECO_GE2J_PTH_0_60_Tag2
functions["GE2J_PTH_120_200_Tag0"] = RECO_GE2J_PTH_120_200_Tag0
functions["GE2J_PTH_120_200_Tag1"] = RECO_GE2J_PTH_120_200_Tag1
functions["GE2J_PTH_120_200_Tag2"] = RECO_GE2J_PTH_120_200_Tag2
functions["GE2J_PTH_60_120_Tag0"] = RECO_GE2J_PTH_60_120_Tag0
functions["GE2J_PTH_60_120_Tag1"] = RECO_GE2J_PTH_60_120_Tag1
functions["GE2J_PTH_60_120_Tag2"] = RECO_GE2J_PTH_60_120_Tag2
functions["PTH_200_300_Tag0"] = RECO_PTH_200_300_Tag0
functions["PTH_200_300_Tag1"] = RECO_PTH_200_300_Tag1
functions["PTH_300_450_Tag0"] = RECO_PTH_300_450_Tag0
functions["PTH_300_450_Tag1"] = RECO_PTH_300_450_Tag1
functions["PTH_450_650_Tag0"] = RECO_PTH_450_650_Tag0
functions["PTH_GT650_Tag0"] = RECO_PTH_GT650_Tag0
functions["THQ_LEP"] = RECO_THQ_LEP
functions["TTH_HAD_PTH_0_60_Tag0"] = RECO_TTH_HAD_PTH_0_60_Tag0
functions["TTH_HAD_PTH_0_60_Tag1"] = RECO_TTH_HAD_PTH_0_60_Tag1
functions["TTH_HAD_PTH_0_60_Tag2"] = RECO_TTH_HAD_PTH_0_60_Tag2
functions["TTH_HAD_PTH_120_200_Tag0"] = RECO_TTH_HAD_PTH_120_200_Tag0
functions["TTH_HAD_PTH_120_200_Tag1"] = RECO_TTH_HAD_PTH_120_200_Tag1
functions["TTH_HAD_PTH_120_200_Tag2"] = RECO_TTH_HAD_PTH_120_200_Tag2
functions["TTH_HAD_PTH_120_200_Tag3"] = RECO_TTH_HAD_PTH_120_200_Tag3
functions["TTH_HAD_PTH_200_300_Tag0"] = RECO_TTH_HAD_PTH_200_300_Tag0
functions["TTH_HAD_PTH_200_300_Tag1"] = RECO_TTH_HAD_PTH_200_300_Tag1
functions["TTH_HAD_PTH_200_300_Tag2"] = RECO_TTH_HAD_PTH_200_300_Tag2
functions["TTH_HAD_PTH_60_120_Tag0"] = RECO_TTH_HAD_PTH_60_120_Tag0
functions["TTH_HAD_PTH_60_120_Tag1"] = RECO_TTH_HAD_PTH_60_120_Tag1
functions["TTH_HAD_PTH_60_120_Tag2"] = RECO_TTH_HAD_PTH_60_120_Tag2
functions["TTH_HAD_PTH_GT300_Tag0"] = RECO_TTH_HAD_PTH_GT300_Tag0
functions["TTH_HAD_PTH_GT300_Tag1"] = RECO_TTH_HAD_PTH_GT300_Tag1
functions["TTH_LEP_PTH_0_60_Tag0"] = RECO_TTH_LEP_PTH_0_60_Tag0
functions["TTH_LEP_PTH_0_60_Tag1"] = RECO_TTH_LEP_PTH_0_60_Tag1
functions["TTH_LEP_PTH_0_60_Tag2"] = RECO_TTH_LEP_PTH_0_60_Tag2
functions["TTH_LEP_PTH_120_200_Tag0"] = RECO_TTH_LEP_PTH_120_200_Tag0
functions["TTH_LEP_PTH_120_200_Tag1"] = RECO_TTH_LEP_PTH_120_200_Tag1
functions["TTH_LEP_PTH_200_300_Tag0"] = RECO_TTH_LEP_PTH_200_300_Tag0
functions["TTH_LEP_PTH_60_120_Tag0"] = RECO_TTH_LEP_PTH_60_120_Tag0
functions["TTH_LEP_PTH_60_120_Tag1"] = RECO_TTH_LEP_PTH_60_120_Tag1
functions["TTH_LEP_PTH_60_120_Tag2"] = RECO_TTH_LEP_PTH_60_120_Tag2
functions["TTH_LEP_PTH_GT300_Tag0"] = RECO_TTH_LEP_PTH_GT300_Tag0
functions["VBFLIKEGGH_Tag0"] = RECO_VBFLIKEGGH_Tag0
functions["VBFLIKEGGH_Tag1"] = RECO_VBFLIKEGGH_Tag1
functions["VBFTOPO_BSM_Tag0"] = RECO_VBFTOPO_BSM_Tag0
functions["VBFTOPO_BSM_Tag1"] = RECO_VBFTOPO_BSM_Tag1
functions["VBFTOPO_JET3VETO_HIGHMJJ_Tag0"] = RECO_VBFTOPO_JET3VETO_HIGHMJJ_Tag0
functions["VBFTOPO_JET3VETO_HIGHMJJ_Tag1"] = RECO_VBFTOPO_JET3VETO_HIGHMJJ_Tag1
functions["VBFTOPO_JET3VETO_LOWMJJ_Tag0"] = RECO_VBFTOPO_JET3VETO_LOWMJJ_Tag0
functions["VBFTOPO_JET3VETO_LOWMJJ_Tag1"] = RECO_VBFTOPO_JET3VETO_LOWMJJ_Tag1
functions["VBFTOPO_JET3_HIGHMJJ_Tag0"] = RECO_VBFTOPO_JET3_HIGHMJJ_Tag0
functions["VBFTOPO_JET3_HIGHMJJ_Tag1"] = RECO_VBFTOPO_JET3_HIGHMJJ_Tag1
functions["VBFTOPO_JET3_LOWMJJ_Tag0"] = RECO_VBFTOPO_JET3_LOWMJJ_Tag0
functions["VBFTOPO_JET3_LOWMJJ_Tag1"] = RECO_VBFTOPO_JET3_LOWMJJ_Tag1
functions["VBFTOPO_VHHAD_Tag0"] = RECO_VBFTOPO_VHHAD_Tag0
functions["VBFTOPO_VHHAD_Tag1"] = RECO_VBFTOPO_VHHAD_Tag1
functions["VH_MET_Tag0"] = RECO_VH_MET_Tag0
functions["VH_MET_Tag1"] = RECO_VH_MET_Tag1
functions["VH_MET_Tag2"] = RECO_VH_MET_Tag2
functions["WH_LEP_PTV_0_75_Tag0"] = RECO_WH_LEP_PTV_0_75_Tag0
functions["WH_LEP_PTV_0_75_Tag1"] = RECO_WH_LEP_PTV_0_75_Tag1
functions["WH_LEP_PTV_75_150_Tag0"] = RECO_WH_LEP_PTV_75_150_Tag0
functions["WH_LEP_PTV_75_150_Tag1"] = RECO_WH_LEP_PTV_75_150_Tag1
functions["WH_LEP_PTV_GT150_Tag0"] = RECO_WH_LEP_PTV_GT150_Tag0
functions["ZH_LEP_Tag0"] = RECO_ZH_LEP_Tag0
functions["ZH_LEP_Tag1"] = RECO_ZH_LEP_Tag1
