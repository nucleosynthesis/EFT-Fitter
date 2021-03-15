import numpy 
import sys 
# object that returns a radial basis spline 

class rbf_spline(ndim=1):
    self._ndim = ndim
    self._initialised = False 
    
    
    def initialise(self,input_points,target_col,eps=10):
        # This is the basic function 
        # Expects a list of dictionaries (inputs_points) 
        # eg [{"x":1,"y":2,"f":4},{"x":3,"y":1,"f":6} ... ]

        # Each dictionary should have ndim+1 keys
        # (probably a better way to structure as a dataframe)

        self._eps = eps;  

        self._M = len(input_points)
        if self._M < 1 : sys.exit("Error - rbf_spline must be initialised with at least one basis point")
        
        axis_pts_ = self.M_**1./self._ndim;

        self._f_vec = [input_points[i][target_col] for i in range(self._M)]
        self._parameter_keys = list(filter(lambda k: k!=target_col, input_points[0].keys()))
        self.r_map =  [{k: [min([input_points[i][k] for i in range(self._M)]), max([input_points[i][k] for i in range(self._M)])] for k in self._parameter_keys} ]
    
    def initialise(self,input_file,target_col):

