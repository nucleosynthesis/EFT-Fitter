from collections import OrderedDict as od

name = "HIG19015"

from tools.likelihoods import rbf_spline,splinesum

params = [
'0J_PTH_0_10_Tag0', '0J_PTH_0_10_Tag1', '0J_PTH_0_10_Tag2', '0J_PTH_GT10_Tag0', '0J_PTH_GT10_Tag1', '0J_PTH_GT10_Tag2', '1J_PTH_0_60_Tag0', '1J_PTH_0_60_Tag1', '1J_PTH_0_60_Tag2', '1J_PTH_120_200_Tag0', '1J_PTH_120_200_Tag1', '1J_PTH_120_200_Tag2', '1J_PTH_60_120_Tag0', '1J_PTH_60_120_Tag1', '1J_PTH_60_120_Tag2', 'GE2J_PTH_0_60_Tag0', 'GE2J_PTH_0_60_Tag1', 'GE2J_PTH_0_60_Tag2', 'GE2J_PTH_120_200_Tag0', 'GE2J_PTH_120_200_Tag1', 'GE2J_PTH_120_200_Tag2', 'GE2J_PTH_60_120_Tag0', 'GE2J_PTH_60_120_Tag1', 'GE2J_PTH_60_120_Tag2', 'PTH_200_300_Tag0', 'PTH_200_300_Tag1', 'PTH_300_450_Tag0', 'PTH_300_450_Tag1', 'PTH_450_650_Tag0', 'PTH_GT650_Tag0', 'THQ_LEP', 'TTH_HAD_PTH_0_60_Tag0', 'TTH_HAD_PTH_0_60_Tag1', 'TTH_HAD_PTH_0_60_Tag2', 'TTH_HAD_PTH_120_200_Tag0', 'TTH_HAD_PTH_120_200_Tag1', 'TTH_HAD_PTH_120_200_Tag2', 'TTH_HAD_PTH_120_200_Tag3', 'TTH_HAD_PTH_200_300_Tag0', 'TTH_HAD_PTH_200_300_Tag1', 'TTH_HAD_PTH_200_300_Tag2', 'TTH_HAD_PTH_60_120_Tag0', 'TTH_HAD_PTH_60_120_Tag1', 'TTH_HAD_PTH_60_120_Tag2', 'TTH_HAD_PTH_GT300_Tag0', 'TTH_HAD_PTH_GT300_Tag1', 'TTH_LEP_PTH_0_60_Tag0', 'TTH_LEP_PTH_0_60_Tag1', 'TTH_LEP_PTH_0_60_Tag2', 'TTH_LEP_PTH_120_200_Tag0', 'TTH_LEP_PTH_120_200_Tag1', 'TTH_LEP_PTH_200_300_Tag0', 'TTH_LEP_PTH_60_120_Tag0', 'TTH_LEP_PTH_60_120_Tag1', 'TTH_LEP_PTH_60_120_Tag2', 'TTH_LEP_PTH_GT300_Tag0', 'VBFLIKEGGH_Tag0', 'VBFLIKEGGH_Tag1', 'VBFTOPO_BSM_Tag0', 'VBFTOPO_BSM_Tag1', 'VBFTOPO_JET3VETO_HIGHMJJ_Tag0', 'VBFTOPO_JET3VETO_HIGHMJJ_Tag1', 'VBFTOPO_JET3VETO_LOWMJJ_Tag0', 'VBFTOPO_JET3VETO_LOWMJJ_Tag1', 'VBFTOPO_JET3_HIGHMJJ_Tag0', 'VBFTOPO_JET3_HIGHMJJ_Tag1', 'VBFTOPO_JET3_LOWMJJ_Tag0', 'VBFTOPO_JET3_LOWMJJ_Tag1', 'VBFTOPO_VHHAD_Tag0', 'VBFTOPO_VHHAD_Tag1', 'VH_MET_Tag0', 'VH_MET_Tag1', 'VH_MET_Tag2', 'WH_LEP_PTV_0_75_Tag0', 'WH_LEP_PTV_0_75_Tag1', 'WH_LEP_PTV_75_150_Tag0', 'WH_LEP_PTV_75_150_Tag1', 'WH_LEP_PTV_GT150_Tag0', 'ZH_LEP_Tag0', 'ZH_LEP_Tag1'
]

splines = []
for p in params : 
 spline = rbf_spline(1,use_scipy_interp=False)
 # we're missing quite a few points so better to use scipy's interpolate
 spline.initialise("inputs/HIG_19_015/root/output_ic_crab_full_wotheory_RECO_"+p+"_combine.txt","chi2",2,parameter_rename=[p],rescaleAxis=True)
 splines.append(spline)
ch2 = splinesum(splines)
X = {}
X['likelihood'] = {"likelihood":ch2}


