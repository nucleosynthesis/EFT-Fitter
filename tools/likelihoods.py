#from functools import lru_cache
import numpy as np
import sys
from numpy.core.function_base import linspace

from scipy import interpolate
from numpy.core.numeric import NaN 
# object that returns a radial basis spline 

class rbf_spline:
    def __init__(self,ndim=1,use_scipy_interp=False):
        self._ndim = ndim
        self._initialised = False 
        self._use_scipy_interp = use_scipy_interp
    
    def _initialise(self,input_points,target_col,eps,rescaleAxis):
        # This is the basic function 
        # Expects a list of dictionaries (inputs_points) 
        # eg [{"x":1,"y":2,"f":4},{"x":3,"y":1,"f":6} ... ]

        # Each dictionary should have ndim+1 keys
        # (probably a better way to structure as a dataframe)

        self._eps = eps  
        self._rescaleAxis = rescaleAxis

        self._M = len(input_points)
        if self._M < 1 : sys.exit("Error - rbf_spline must be initialised with at least one basis point")
        
        self._parameter_keys = list(filter(lambda k: k!=target_col, input_points[0].keys()))
        if self._ndim=="auto" : self._ndim = len(self._parameter_keys)
        if self._ndim!=len(self._parameter_keys): 
            sys.exit("Error - initialise given points with more dimensions (%g) than ndim (%g)"%(len(self._parameter_keys),self._ndim))

        self._axis_pts = self._M**(1./self._ndim)

        self._v_map = [ {k:v for k,v in a.items() if k!=target_col} for a in input_points ]
        self._r_map =  {k: [min([input_points[i][k] for i in range(self._M)]), \
                            max([input_points[i][k] for i in range(self._M)])] \
                            for k in self._parameter_keys}
        f_vec = [input_points[i][target_col] for i in range(self._M)]
        self._f_vec = f_vec
        
        max_f_vec = max([abs(f) for f in f_vec])
        f_vec = np.array(f_vec)/max_f_vec
        self.calculateWeights(f_vec)
        self._max_f_vec = max_f_vec
    
    def _initialise_scipy(self,input_points,target_col): 
        self._M = len(input_points)
        self._v_map = [ {k:v for k,v in a.items() if k!=target_col} for a in input_points ]
        f_vec = [input_points[i][target_col] for i in range(self._M)]
        self._f_vec = f_vec
        self._parameter_keys = list(filter(lambda k: k!=target_col, input_points[0].keys()))
        if self._ndim=="auto" : self._ndim = len(self._parameter_keys)
        if self._ndim!=len(self._parameter_keys): 
            sys.exit("Error - initialise given points with more dimensions (%g) than ndim (%g)"%(len(self._parameter_keys),self._ndim))
        if self._ndim!=1 : sys.exit("Error - can only use scipy interpolate for 1D spline")
        points = [ [input_points[i][self._parameter_keys[0]],input_points[i][target_col]] for i in range(self._M) ]
        points.sort()
        f_vec = [p[1] for p in points]
        p_vec = [p[0] for p in points] 
        self._f = interpolate.interp1d(p_vec,f_vec,"cubic")
        #self._f = interpolate.Rbf(p_vec,f_vec,function='gaussian',epsilon = eps) <- seems to work less well than the implementation above
        self._r_map =  {k: [min([input_points[i][k] for i in range(self._M)]), \
                            max([input_points[i][k] for i in range(self._M)])] \
                            for k in self._parameter_keys}
        self._initialised = True

    def initialise(self,input_points,target_col,eps=10,rescaleAxis=True,parameter_rename=[]):
        if type(input_points)==str: 
            fi = open(input_points,"r")
            keys = []
            input_points = []
            for i,line in enumerate(fi.readlines()): 
                vals = line.split()
                if not len(vals): continue
                if i==0 : 
                  keys = vals
                  if len(parameter_rename): 
                   for i in range(len(parameter_rename)): keys[i]=parameter_rename[i] 
                else: 
                  pt = {keys[j]:float(vals[j])  for j in range(len(keys))}
                  input_points.append(pt)
        if self._use_scipy_interp : self._initialise_scipy(input_points,target_col)
        else: self._initialise(input_points,target_col,eps,rescaleAxis)
    
    def diff(self,a,b,k): # note this is not quite sqrt(diff) but really d(diff2)/dX
        v=a-b
        c = 1. 
        if self._rescaleAxis: c = (self._axis_pts)*(self._axis_pts)/((self._r_map[k][1]-self._r_map[k][0])*(self._r_map[k][1]-self._r_map[k][0]))
        return v*c

    def diff2(self,a,b,k):
        if self._rescaleAxis: v=(self._axis_pts)*(a-b)/(self._r_map[k][1]-self._r_map[k][0])
        else: v=a-b
        return v*v  

    def getDistSquare(self,i, j):
        dk2 = np.array([ self.diff2(self._v_map[i][k],self._v_map[j][k],k) for k in self._parameter_keys ])
        return sum(dk2)

    def getDistFromSquare(self,point, i):
        dk2 = np.array([ self.diff2(self._v_map[i][k],point[k],k) for k in self._parameter_keys ])
        return sum(dk2)

    def getGradDistFrom(self,point,i,param): 
        return 2*self.diff(point[param],self._v_map[i][param],param)

    """
    def radialFunc(self,d2): 
        return np.sqrt((d2/self._eps)**2 + 1)
    """
    def radialFunc(self,d2):
        expo = (d2/(self._eps*self._eps))
        ret_val = np.exp(-1.*expo)  
        #if ret_val < 1e-12 : return 0
        return ret_val
    

    def evaluate(self,point):
        if not self._initialised:
            print("Error - must first initialise spline with set of points before calling evaluate()") 
            return NaN
        #if not set(point.keys())==set(self._parameter_keys): 
        #    print ("Error - must have same variable labels, you provided - ",point.keys(),", I only know about - ",self._parameter_keys)
        #    return NaN
        # check bounds of points and set to edge if its there 
        """
        for p in point.keys():
            if   point[p] < self._r_map[p][0]: 
              #print("ERROR - out of range (<) for ",p,"=", point[p], "bounds=",self._r_map[p])
              #return 1e3 
              point[p] = self._r_map[p][0]+1e-3
            elif point[p] > self._r_map[p][1]: 
              #print("ERROR - out of range (>) for ",p,"=", point[p], "bounds=",self._r_map[p])
              #return 1e3 
              point[p] = self._r_map[p][1]-1e-3
        """
        if self._use_scipy_interp: 
            return self._f(point[self._parameter_keys[0]])
        vals_sum = self._weights.dot(np.array([self.radialFunc(self.getDistFromSquare(point,i)) for i in range(self._M)]))
        #vals*=self._max_f_vec
        return self._max_f_vec*vals_sum

    def evaluate_grad(self,point,param):
        if not self._initialised:
            print("Error - must first initialise spline with set of points before calling evaluate_grad()") 
            return NaN
        if param not in point.keys(): return 0 # this is so I can be lazy later
        """
        for p in point.keys():
            if   point[p] < self._r_map[p][0]: 
              #print("ERROR - out of range (<) for ",p,"=", point[p], "bounds=",self._r_map[p])
              #return 1e3 
              point[p] = self._r_map[p][0]+1e-3
            elif point[p] > self._r_map[p][1]: 
              #print("ERROR - out of range (>) for ",p,"=", point[p], "bounds=",self._r_map[p])
              #return 1e3 
              point[p] = self._r_map[p][1]-1e-3
        """
        if self._use_scipy_interp: 
            sys.exit("no gradient for scipy interpolate (yet?)")

        #vals = self._weights * np.array([self.radialFunc(self.getDistFromSquare(point,i)) for i in range(self._M)])
        #vals*=self._max_f_vec
        vals_sum = self._weights.dot(np.array([ self.radialFunc(self.getDistFromSquare(point,i))*self.getGradDistFrom(point,i,param) for i in range(self._M)] ) )
        return -1./(self._eps*self._eps)*self._max_f_vec*vals_sum

    def calculateWeights(self,f) : 
        A = np.array([np.zeros(self._M) for i in range(self._M)])

        for i in range(self._M):
            A[i][i]=1.
            for j in range(i+1,self._M):
                d2  = self.getDistSquare(i,j)
                rad = self.radialFunc(d2)
                A[i][j] = rad
                A[j][i] = rad
        
        B = np.array(f)
        self._weights = np.linalg.solve(A,B)
        self._initialised=True
    
    def getParameters(self):
        return self._parameter_keys[:]

    def getMinMax(self,param):
        return self._r_map[param][0],self._r_map[param][1]

    def getF(self):
        return self._f_vec[:]

    def getPoints(self,axis):
        # return a grid of points up to 2D
        if self._ndim > 2 : 
          print("Error, can only return grids for ndim<=2 splines")
          return []
        # return a set of points given axis ["x","y"]
        if self._ndim==1:
          return [self._v_map[i][axis[0]] for i in range(self._M)]
        else:  
          return [[self._v_map[i][axis[0]],self._v_map[i][axis[1]]] for i in range(self._M)]

class splinesum:
    def __init__(self,splines): 
        self._parameter_keys = []
        [self._parameter_keys.extend(sp.getParameters()) for sp in splines]
        if len(set(self._parameter_keys)) != len(self._parameter_keys): 
            print("Error - splinesum cannot have same parameter in multiple splines! ...")
            for i,sp in enumerate(splines): print("spline %d - "%i,sp.getParameters())
            sys.exit()

        self._splines = splines 
        self._N = len(self._splines)

        self._pmap = {i: self._splines[i].getParameters() for i in range(self._N)}

    def getParameters(self):
        return self._parameter_keys[:]
    
    def getNSplines(self):
        return self._N
    
    def getSpline(self,i):
        return self._splines[i]
    #@lru_cache(maxsize=20)
    def evaluate(self,point):
        return sum(np.asarray([ sp.evaluate({k:v for k,v in point.items() if k in sp.getParameters()}) for sp in self._splines ]))
    
    def evaluate_grad(self,point):
        ret_map={p:sum([sp.evaluate_grad({k:v for k,v in point.items() if k in sp.getParameters()},p) \
            for sp in filter(lambda x : p in x.getParameters(),self._splines)]) for p in self._parameter_keys}
        return ret_map
        #return {p: sum([sp.evaluate_grad({k:v for k,v in point.items() if k in sp.getParameters()},p) for sp in self._splines ])

def quad(r):
        return 3*((r-2)**2)
def quad_g(r):
        return (6*r)-12
if __name__ == "__main__":


    # very simple example of it 
    spline = rbf_spline(1,use_scipy_interp=False)
    # do a very basic quadratic - 3*(r-2)^2
    ip  = [{'chi2': quad(r), 'r': r} for r in np.linspace(-3,5,30)  ]
    #spline.initialise(ip,"chi2")
    spline.initialise(ip,"chi2",5)
    ch2 = splinesum([spline])
    import matplotlib.pyplot as plt
    x = [ip[i]["r"] for i in range(len(ip))]
    yfix = [ip[i]["chi2"] for i in range(len(ip))]
    xx = np.linspace(min(x)+0.01,max(x)-0.01,num=100)
    yint = [spline.evaluate({"r":xi}) for xi in xx]
    yint_g = [ch2.evaluate_grad({"r":xi})["r"] for xi in xx]
    
    #[ spline.evaluate_grad({"r":xi},"r") for xi in xx ]
    yint_ga = [ quad_g(r) for r in xx ]
    
    plt.plot(x,yfix,marker=".",linestyle="None")
    plt.plot(xx,yint,color="red")
    plt.plot(xx,yint_g,color="blue")
    plt.plot(xx,yint_ga,color="green")

    plt.grid()
    plt.xlabel("r")
    plt.ylabel("chi2")
    plt.show()


