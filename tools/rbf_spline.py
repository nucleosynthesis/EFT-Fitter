import numpy as np
import sys

import pandas as pd
# object that returns a radial basis spline 

# eps is the scale of the basis function for the distance metric. Axes should be of a similar scale. reccomend to use 
# rescaleAxis=True, in this case eps will be in terms of number of nearest points, so this should be an integer > 1 and 
# < the total number of basis points

class rbf_spline:
    def __init__(self,ndim=1):
        self._ndim = ndim
        self._initialised = False
        self.radialFuncs = dict([("gaussian", self.radialGauss),
                                 ("multiquadric", self.radialMultiQuad),
                                 ("inversemultiquadric", self.radialInverseMultiQuad),
                                 ("cubic", self.radialCubic)])
    
    def _initialise(self,input_data,target_col,eps,rescaleAxis):
        self._input_data = input_data
        self._target_col = target_col
        self._input_points = input_data.drop(target_col, 
                                             axis="columns").to_numpy()
        
        self._eps = eps  
        self._rescaleAxis = rescaleAxis

        self._M = len(input_data) # Num points
        if self._M < 1 : sys.exit("Error - rbf_spline must be initialised with at least one basis point")
        
        self._parameter_keys = list(input_data.columns)
        self._parameter_keys.remove(target_col)

        if self._ndim!=len(self._parameter_keys): 
            sys.exit("Error - initialise given points with more dimensions (%g) than ndim (%g)"%(len(self._parameter_keys),self._ndim))

        self._axis_pts = self._M**(1./self._ndim)
       
        self.calculateWeights()

    def initialise(self,input_points,target_col,radial_func="gaussian",eps=10,rescaleAxis=True):
        try:
            self.radialFunc = self.radialFuncs[radial_func]
        except KeyError:
            sys.exit("Error - function '%s' not in '%s'"%(radial_func, list(self.radialFuncs.keys())))
        self._initialise(input_points,target_col,eps,rescaleAxis)
    
    def initialise_text(self,input_file,target_col,radial_func="gaussian",eps=10,rescaleAxis=True):
        df = pd.read_csv(input_file,index_col=False,delimiter=' ')
        self.initialise(df,target_col,radial_func,eps,rescaleAxis)
    
    def diff2(self, points1, points2):
        # The interpolator must have been initialised on points2
        v = points1[:, np.newaxis, :] - points2[np.newaxis, :, :]
        if self._rescaleAxis: v=self._axis_pts*v/(np.max(points2, axis=0) - np.min(points2, axis=0))
        return np.power(v, 2)

    def getDistSquare(self, col):
        return self.diff2(col, col)
        
    def getDistFromSquare(self, point, inp):
        dk2 = np.sum(self.diff2(point, inp), axis=-1).flatten()
        return dk2

    def getRadialArg(self, d2):
        return (d2/(self._eps*self._eps))

    def radialGauss(self,d2):
        return np.e**(-1*self.getRadialArg(d2))
    
    def radialMultiQuad(self,d2):
        return np.sqrt(1+self.getRadialArg(d2)) 
        
    def radialInverseMultiQuad(self, d2):
        return 1/self.radialMultiQuad(self.getRadialArg(d2))

    def radialCubic(self, d2):
        return np.power(self.getRadialArg(d2), 3/2)

    def evaluate(self,point):
        if not self._initialised:
            print("Error - must first initialise spline with set of points before calling evaluate()") 
            return np.nan
        if not set(point.keys())==set(self._parameter_keys): 
            print ("Error - must have same variable labels, you provided - ",point.keys(),", I only know about - ",self._parameter_keys)
            return np.nan
        vals = self.radialFunc(self.getDistFromSquare(point.to_numpy(), self._input_points))
        weighted_vals = self._weights * vals
        return sum(weighted_vals)

<<<<<<< HEAD
    def calculateWeights(self) : 
        inp = self._input_points
        B = self._input_data[self._target_col].to_numpy()
        d2 = np.sum(self.diff2(inp, inp), axis=2)
        A = self.radialFunc(d2) 
        np.fill_diagonal(A, 1)
    
=======
        p   = np.array([np.array(point[k]) for k in self._parameter_keys])
        dx  = np.array((np.array(self._v)-p)/self._k)

        dx2 = np.array([dx[i].dot(dx[i]) for i in range(self._M)])
        phi = np.array(self.vectorized_radialFunc(dx2)) 
        vals_sum = self._weights.dot(phi)
        return self._max_f_vec*vals_sum

    def evaluate_grad(self,point,param):
        if not self._initialised:
            print("Error - must first initialise spline with set of points before calling evaluate_grad()") 
            return NaN
        if param not in point.keys(): return 0 # this is so I can be lazy later

        if self._use_scipy_interp: 
            sys.exit("no gradient for scipy interpolate (yet?)")

        p   = np.array([np.array(point[k]) for k in self._parameter_keys])
        dx  = np.array((np.array(self._v)-p)/self._k)
        dx2 = np.array([dx[i].dot(dx[i]) for i in range(self._M)]) #(self.vectorized_squarepoint(dx))
        parameter_index = self._parameter_keys.index(param)

        vpar = np.array([self._v[i][parameter_index] for i in range(self._M)])
        ddx  = 2*(p[parameter_index]-vpar)/(self._sk[parameter_index])
        # for generic radial function, change below to use "self.vectorised_radialFunc_grad(dx2)" once implemented for all radial functions,  
        # and remove -1 from the return at the end! 
        dphi   = ddx*self.vectorized_radialFunc(dx2)  # this is true ONLY for the case of the gaussian radial function! 
        vals_sum = self._weights.dot(dphi)
        return -1./(self._eps*self._eps)*self._max_f_vec*vals_sum

    def calculateWeights(self,f) : 
        A = np.array([np.zeros(self._M) for i in range(self._M)])

        for i in range(self._M):
            A[i][i]=1.
            for j in range(i+1,self._M):
                d2  = self.getDistSquare(i,j)
                rad = self.vectorized_radialFunc(d2)
                A[i][j] = rad
                A[j][i] = rad
        
        B = np.array(f)
>>>>>>> 91ac1237e47d374bb15182aa7e0d36e09b8f6943
        self._weights = np.linalg.solve(A,B)
        self._initialised=True


"""
# very simple example of it 
spline = rbf_spline(1)
import matplotlib.pyplot as plt
# from a dataframe

# Two types on input, dataframe or text file. The latter converts to a dataframe first. 

data = {'chi2':[0,0.500583,0.864236,1.38561,2.12717,0.0942076,0.195518,0.325861,0.481119,0.657985,0.853566,1.06564,1.29224,1.53174,1.78274,2.04419,2.31472,2.59386,2.8806,3.17434,3.47449,3.78047,4.0918,4.40799,4.72866,5.05341,5.38184,5.7136,6.04834,6.38599,6.72614,7.06876,7.41364,7.7604,8.10892,0.255522,0.102403,0.0218385,0.000177219,0.0269786]
        ,'r':[1.6 , 0.9 , 0.7, 0.5 ,0.3 ,2.1 ,2.3 ,2.5 ,2.7 ,2.9 ,3.1 ,3.3 ,3.5 ,3.7 ,3.9 ,4.1 ,4.3 ,4.5 ,4.7 ,4.9 ,5.1 ,5.3 ,5.5 ,5.7 ,5.9 ,6.1 ,6.3 ,6.5 ,6.7 ,6.9 ,7.1 ,7.3 ,7.5 ,7.7 ,7.9 ,1.1 ,1.3 ,1.5 ,1.7 ,1.9]}


df = pd.DataFrame(data=data)

#spline.initialise(df,'chi2', radial_func="gaussian", eps=2)
spline.initialise_text('test.txt','chi2', radial_func="gaussian", eps=2)

x = data["r"] 
yfix = data["chi2"] 
xx = np.linspace(min(x)+0.01,max(x)-0.01,num=100)
yint = [spline.evaluate(pd.DataFrame({"r":[xi]})) for xi in xx]
plt.plot(x,yfix,marker=".",linestyle="None")
plt.plot(xx,yint,color="red")
plt.xlabel("r")
plt.ylabel("chi2")
plt.show()
"""
