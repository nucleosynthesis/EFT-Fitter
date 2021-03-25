from collections import OrderedDict as od

name = "test_lh"

from tools.likelihoods import rbf_spline,splinesum

spline = rbf_spline(2)
#spline = rbf_spline(1,use_scipy_interp=True)

spline.initialise("out-2D.txt","chi2",600)
# can also put  numbers as columns in text file and read from there 
# param f 
# 1.0     0.2
# 3.4     0.5
# ...
# spline.initialise("output.txt","chi2")
X = {}
chi2 = splinesum([spline])
X['likelihood'] = {"likelihood":chi2}

