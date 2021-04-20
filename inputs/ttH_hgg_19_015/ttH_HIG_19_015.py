from collections import OrderedDict as od

name = "ttH_hgg_theory"

from tools.likelihoods import rbf_spline,splinesum

params = [
"TTH_HAD_PTH_60_120_Tag2"  
, "THQ_LEP"                  
, "TTH_HAD_PTH_GT300_Tag0"  
, "TTH_HAD_PTH_0_60_Tag0"    
, "TTH_HAD_PTH_GT300_Tag1" 
, "TTH_HAD_PTH_0_60_Tag1"    
, "TTH_LEP_PTH_0_60_Tag0" 
, "TTH_HAD_PTH_0_60_Tag2"    
, "TTH_LEP_PTH_0_60_Tag1" 
, "TTH_HAD_PTH_120_200_Tag0" 
, "TTH_LEP_PTH_0_60_Tag2" 
, "TTH_HAD_PTH_120_200_Tag1" 
, "TTH_LEP_PTH_120_200_Tag0" 
, "TTH_HAD_PTH_120_200_Tag2" 
, "TTH_LEP_PTH_120_200_Tag1" 
, "TTH_HAD_PTH_120_200_Tag3" 
, "TTH_LEP_PTH_200_300_Tag0" 
, "TTH_HAD_PTH_200_300_Tag0" 
, "TTH_LEP_PTH_60_120_Tag0" 
, "TTH_HAD_PTH_200_300_Tag1" 
, "TTH_LEP_PTH_60_120_Tag1" 
, "TTH_HAD_PTH_200_300_Tag2" 
, "TTH_LEP_PTH_60_120_Tag2" 
, "TTH_HAD_PTH_60_120_Tag0"  
, "TTH_LEP_PTH_GT300_Tag0" 
, "TTH_HAD_PTH_60_120_Tag1"
]

splines = []
for p in params : 
 spline = rbf_spline(1,use_scipy_interp=False)
 spline.initialise("inputs/ttH_hgg_19_015/"+p+"_combine.txt","chi2",3)
 splines.append(spline)
ch2 = splinesum(splines)
X = {}
X['likelihood'] = {"likelihood":ch2}


