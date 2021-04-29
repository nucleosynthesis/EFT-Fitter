import numpy as np
import sys

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
    
    
    def diff2(self,a,b,k):
        v=0
        if self._rescaleAxis: v=(self._axis_pts)*(a-b)/(self._r_map[k][1]-self._r_map[k][0])
        #if self._rescaleAxis: v=(a-b)/abs((self._r_map[k][1]-self._r_map[k][0]))
        else: v=a-b
        return v*v  

    def getDistSquare(self,i, j):
        dk2 = [ self.diff2(self._v_map[i][k],self._v_map[j][k],k) for k in self._parameter_keys ]
        return sum(dk2)

    def getDistFromSquare(self,point, i):
        dk2 = [ self.diff2(self._v_map[i][k],point[k],k) for k in self._parameter_keys ]
        return sum(dk2)

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
        if not set(point.keys())==set(self._parameter_keys): 
            print ("Error - must have same variable labels, you provided - ",point.keys(),", I only know about - ",self._parameter_keys)
            return NaN
        # check bounds of points and set to edge if its there 
        for p in point.keys():
            if   point[p] < self._r_map[p][0]: 
              #print("ERROR - out of range (<) for ",p,"=", point[p])
              #return 1e3 
              point[p] = self._r_map[p][0]+1e-3
            elif point[p] > self._r_map[p][1]: 
              #print("ERROR - out of range (>) for ",p,"=", point[p])
              #return 1e3 
              point[p] = self._r_map[p][1]-1e-3
        if self._use_scipy_interp: 
            return self._f(point[self._parameter_keys[0]])
        vals = self._weights * np.array([self.radialFunc(self.getDistFromSquare(point,i)) for i in range(self._M)])
        vals*=self._max_f_vec
        return sum(vals)

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

    def evaluate(self,point):
        #sm = 0
        #for i,sp in enumerate(splines): 
        #    spoint = {k:v for k,v in point.items() if k in sp.getParameters()} 
        #    sm += sp.evaluate(spoint)
        #print("splinesum evaluate point",point)
        #print([ sp.evaluate({k:v for k,v in point.items() if k in sp.getParameters()}) for sp in self._splines ])
        #print("evaluating at point",point)
        #print("chi2vals in sum  - ", [ sp.evaluate({k:v for k,v in point.items() if k in sp.getParameters()}) for sp in self._splines ])
        return sum([ sp.evaluate({k:v for k,v in point.items() if k in sp.getParameters()}) for sp in self._splines ])

"""
# very simple example of it 
spline = rbf_spline(1,use_scipy_interp=True)
ip  = [{'chi2': 0.0, 'r': 0.15625189244747162}, {'chi2': 22.59496307373047, 'r': -0.0949999988079071}, {'chi2': 21.038896560668945, 'r': -0.08500000089406967}, {'chi2': 19.512514114379883, 'r': -0.07500000298023224}, {'chi2': 18.018226623535156, 'r': -0.06499999761581421}, {'chi2': 16.55982780456543, 'r': -0.054999999701976776}, {'chi2': 15.140303611755371, 'r': -0.044999998062849045}, {'chi2': 13.763101577758789, 'r': -0.03500000014901161}, {'chi2': 12.431910514831543, 'r': -0.02499999850988388}, {'chi2': 11.150689125061035, 'r': -0.014999998733401299}, {'chi2': 9.92322063446045, 'r': -0.004999998491257429}, {'chi2': 8.753515243530273, 'r': 0.005000002216547728}, {'chi2': 7.645490646362305, 'r': 0.015000002458691597}, {'chi2': 6.602964878082275, 'r': 0.02500000223517418}, {'chi2': 5.62955379486084, 'r': 0.03500000387430191}, {'chi2': 4.728558540344238, 'r': 0.04500000178813934}, {'chi2': 3.9028549194335938, 'r': 0.055000003427267075}, {'chi2': 3.154820442199707, 'r': 0.0650000050663948}, {'chi2': 2.4862303733825684, 'r': 0.07500000298023224}, {'chi2': 1.8981651067733765, 'r': 0.08500000834465027}, {'chi2': 1.3911701440811157, 'r': 0.0950000062584877}, {'chi2': 0.965129554271698, 'r': 0.10500000417232513}, {'chi2': 0.6190124750137329, 'r': 0.11500000208616257}, {'chi2': 0.3516332507133484, 'r': 0.125}, {'chi2': 0.16055263578891754, 'r': 0.13500000536441803}, {'chi2': 0.04442401975393295, 'r': 0.14500001072883606}, {'chi2': 0.000573080382309854, 'r': 0.1550000011920929}, {'chi2': 0.026336580514907837, 'r': 0.16500000655651093}, {'chi2': 0.11885059624910355, 'r': 0.17500001192092896}, {'chi2': 0.27496013045310974, 'r': 0.1850000023841858}, {'chi2': 0.4920821785926819, 'r': 0.19500000774860382}, {'chi2': 0.7674822211265564, 'r': 0.20500001311302185}, {'chi2': 1.0982145071029663, 'r': 0.2150000035762787}, {'chi2': 1.4816529750823975, 'r': 0.22500000894069672}, {'chi2': 1.9153952598571777, 'r': 0.23500001430511475}, {'chi2': 2.397061824798584, 'r': 0.24500000476837158}, {'chi2': 2.9244658946990967, 'r': 0.2550000250339508}, {'chi2': 3.4955008029937744, 'r': 0.26500001549720764}, {'chi2': 4.108160018920898, 'r': 0.2750000059604645}, {'chi2': 4.760610103607178, 'r': 0.2850000262260437}, {'chi2': 5.4509968757629395, 'r': 0.29500001668930054}]
#spline.initialise(ip,"chi2")
spline.initialise("output.txt","chi2",10)
import matplotlib.pyplot as plt
x = [ip[i]["r"] for i in range(len(ip))]
yfix = [ip[i]["chi2"] for i in range(len(ip))]
xx = np.linspace(min(x)+0.01,max(x)-0.01,num=100)
yint = [spline.evaluate({"r":xi}) for xi in xx]
plt.plot(x,yfix,marker=".",linestyle="None")
plt.plot(xx,yint,color="red")
plt.xlabel("r")
plt.ylabel("chi2")
plt.show()
"""

