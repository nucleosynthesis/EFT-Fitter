import sys
import numpy as np
import pandas as pd
import numpy.typing as npt

# -----------------
# Basis functions
# -----------------
class radialGauss():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.exp(-input)
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return -self.evaluate(input)

class radialMultiQuad():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.sqrt(1+input)
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return 1/(2*self.evaluate(input))
    
class radialInverseMultiQuad():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.divide(1, np.sqrt(1+input))
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return -1/(2*np.power(1+input, 3/2))

class radialLinear():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.sqrt(input)
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return 1/(2*self.evaluate(input))

class radialCubic():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.power(input, 3/2)
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return 3*np.sqrt(input)/2
    
class radialQuintic():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.power(input, 5/2)
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return 5*np.power(input, 3/2)/2

class radialThinPlate():
    def __init__(self) -> None:
        return
    def evaluate(self, 
                 input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return np.multiply(input, np.log(np.sqrt(input)))
    def getDeltaPhi(self, 
                    input: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        return (np.log(input)+1)/2
# -----------------

class rbf_spline:
    def __init__(self, ndim=1) -> None:
        self._ndim = ndim
        self._initialised = False
        self._radialFuncs = dict(
            [("gaussian", radialGauss),
             ("multiquadric", radialMultiQuad),
             ("inversemultiquadric", radialInverseMultiQuad),
             ("linear", radialLinear),
             ("cubic", radialCubic),
             ("quintic", radialQuintic),
             ("thinplate", radialThinPlate)
            ])

    def _initialise(self, input_data: pd.DataFrame, target: str, 
                    eps: float, rescaleAxis: bool) -> None:
        # Parse args
        self._input_data = input_data
        self._target_col = target
        self._input_pts = input_data.drop(target, axis="columns").to_numpy()
        self._eps = eps  
        self._rescaleAxis = rescaleAxis
        self._parameter_keys = list(input_data.columns)
        self._parameter_keys.remove(target)

        # Check number of basis points
        self._M = len(input_data)
        if self._M < 1 : 
            sys.exit("Error - At least one basis point is required")
        
        # Check dimensions
        if self._ndim!=len(self._parameter_keys): 
            sys.exit(f"Error - initialise given points with more dimensions " +
                     f"({len(self._parameter_keys)}) than ndim ({self._ndim})")

        # Get scalings by axis (column)
        self._axis_pts = np.power(self._M, 1./self._ndim)
        if self._rescaleAxis:
            self._scale = np.divide(self._axis_pts, 
                                    (np.max(self._input_pts, axis=0) -
                                     np.min(self._input_pts, axis=0)))
        else:
            self._scale = 1

        self.calculateWeights()

    def initialise(self, input_data: pd.DataFrame, target_col: str, 
                   radial_func: str="gaussian", eps: float=10.,
                   rescaleAxis: bool=True) -> None:
        # Get basis function and initialise
        try:
            self.radialFunc = self._radialFuncs[radial_func]()
        except KeyError:
            sys.exit(f"Error - function '{radial_func}' not in " +
                     f"'{list(self._radialFuncs.keys())}'")
        self._initialise(input_data, target_col, eps, rescaleAxis)
        
    def initialise_text(self, input_file: str, target_col, 
                        radial_func: str="gaussian", eps: float=10.,
                        rescaleAxis: bool=True) -> None:
        df = pd.read_csv(input_file, index_col=False, delimiter=' ')
        self.initialise(df,target_col,radial_func,eps,rescaleAxis)
        
    def diff(self, points1: npt.NDArray[np.float32],
             points2: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        # Get diff between two sets of points, pairwise
        v = np.multiply(self._scale, (points1[:, np.newaxis, :] - 
                                      points2[np.newaxis, :, :]))
        return v    
    
    def diff2(self, points1: npt.NDArray[np.float32], 
              points2: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        # Get squared diff between two sets of points, pairwise
        return np.power(self.diff(points1, points2), 2)
    
    def getDistFromSquare(self, point: npt.NDArray[np.float32]):
        # Get distance between a point and the basis points, per axis
        return self.diff2(point, self._input_pts)
    
    def getRadialArg(self, 
                     d2: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        # Get arg to pass to basis functions
        return np.divide(d2, self._eps*self._eps)

    def grad_r2(self, point) -> npt.NDArray[np.float32]:
        # Calculates grad(|r|^2)
        return (2*self.diff(point, self._input_pts)*self._scale/(self._eps*self._eps))

    def evaluate(self, point: pd.DataFrame) -> npt.NDArray[np.float32]:
        # Check input is okay (can be turned off for perfomance)
        if not self._initialised:
            print("Error - must first initialise spline with set of points " + 
                  "before calling evaluate()") 
            return np.array(np.nan)
        if not set(point.keys()) == set(self._parameter_keys): 
            print(f"Error - {point.keys()} must match {self._parameter_keys}")
            return np.array(np.nan)
        
        # Evaluate spline at point
        point_arr = point.to_numpy()
        radial_arg = self.getRadialArg(np.sum(self.getDistFromSquare(point_arr), axis=-1))
        vals = self.radialFunc.evaluate(radial_arg).flatten()
        
        # Get val and grads
        weighted_vals = self._weights * vals
        ret_val = np.sum(weighted_vals)
        
        return ret_val.astype(float)
    
    def evaluate_grad(self, point: pd.DataFrame) -> npt.NDArray[np.float32]:
        # Check input is okay (can be turned off for perfomance)
        if not self._initialised:
            print("Error - must first initialise spline with set of points " + 
                  "before calling evaluate()") 
            return np.array(np.nan)
        if not set(point.keys()) == set(self._parameter_keys): 
            print(f"Error - {point.keys()} must match {self._parameter_keys}")
            return np.array(np.nan)

        # Evaluate spline at point
        point_arr = point.to_numpy()
        radial_arg = self.getRadialArg(self.getDistFromSquare(point_arr))
        delta_phi = self.radialFunc.getDeltaPhi(radial_arg)
        grad_phi = np.linalg.norm(delta_phi, axis=-1)
        grads = self.grad_r2(point_arr) * grad_phi.reshape(1, self._M, 1) * np.sign(delta_phi)
        
        # Get val and grads
        weighted_grads = np.multiply(self._weights.reshape(1, self._M, 1), grads)
        ret_grad = np.sum(weighted_grads, axis=1)
        
        return ret_grad.astype(float)
        
    def calculateWeights(self) -> None: 
        # Solve interpolation matrix equation for weights
        inp = self._input_pts
        B = self._input_data[self._target_col].to_numpy()
        d2 = np.sum(self.diff2(inp, inp), axis=2)
        A = self.radialFunc.evaluate(self.getRadialArg(d2)) 
        np.fill_diagonal(A, 1)
    
        self._interp_mat = A
        self._inv_interp_mat = np.linalg.inv(A)
        self._weights = np.dot(self._inv_interp_mat, B)
        self._initialised = True
        
    def calculateLOOCV(self) -> float:
        # Get leave-one-out cross-validation error, implementing
        # https://doi.org/10.1023/A:1018975909870]
        if not self._initialised:
            print("Error - must first initialise spline with set of points " + 
                  "before calling evaluate()") 
            return np.nan
        
        cost_vec = self._weights / self._inv_interp_mat.diagonal()       
        return np.linalg.norm(cost_vec)
     
"""
# Simple example
import matplotlib.pyplot as plt
import numpy as np
data = {'chi2':[0,0.500583,0.864236,1.38561,2.12717,0.0942076,0.195518,0.325861,0.481119,0.657985,0.853566,1.06564,1.29224,1.53174,1.78274,2.04419,2.31472,2.59386,2.8806,3.17434,3.47449,3.78047,4.0918,4.40799,4.72866,5.05341,5.38184,5.7136,6.04834,6.38599,6.72614,7.06876,7.41364,7.7604,8.10892,0.255522,0.102403,0.0218385,0.000177219,0.0269786]
        ,'r':[1.6 , 0.9 , 0.7, 0.5 ,0.3 ,2.1 ,2.3 ,2.5 ,2.7 ,2.9 ,3.1 ,3.3 ,3.5 ,3.7 ,3.9 ,4.1 ,4.3 ,4.5 ,4.7 ,4.9 ,5.1 ,5.3 ,5.5 ,5.7 ,5.9 ,6.1 ,6.3 ,6.5 ,6.7 ,6.9 ,7.1 ,7.3 ,7.5 ,7.7 ,7.9 ,1.1 ,1.3 ,1.5 ,1.7 ,1.9]}
df = pd.DataFrame(data=data)
spline = rbf_spline(1)
# spline.initialise_text('test.txt','chi2', radial_func="gaussian", eps=2)

# Get the optimal epsilon by minimising the LOOCV error
eps_range = np.linspace(0.01, 5, num=25)
best_eps = None
best_loocv = np.inf
for eps in eps_range:
    spline.initialise(df,'chi2', radial_func="gaussian", eps=eps, rescaleAxis=True)
    loocv = spline.calculateLOOCV()
    if loocv < best_loocv:
        best_loocv = loocv
        best_eps = eps
    print(f"LOOCV: {loocv} for epsilon: {eps}")
print(f"Best epsilon: {best_eps}")
spline.initialise(df,'chi2', radial_func="gaussian", eps=best_eps, rescaleAxis=True)

x = data["r"] 
yfix = data["chi2"] 
xx = np.linspace(min(x)+0.01,max(x)-0.01,num=100)

yint = np.empty(xx.shape)
ygrad_int = np.empty(xx.shape)
for i, xi in enumerate(xx):
    val = spline.evaluate(pd.DataFrame({"r":[xi]}))
    grad = spline.evaluate_grad(pd.DataFrame({"r":[xi]}))
    yint[i] = val
    ygrad_int[i] = grad.flatten()[0]
ygrad_fdiff = np.gradient(yint, xx)

fig, ax = plt.subplots(2, 1, sharex=True)
fig.subplots_adjust(hspace=0)
ax[0].plot(x,yfix,marker=".",linestyle="None")
ax[0].plot(xx,yint,color="red")
ax[0].set_ylabel("chi2")

ax[1].plot(xx,ygrad_fdiff)
ax[1].plot(xx,ygrad_int,color="red")
ax[1].set_xlabel("r")
ax[1].set_ylabel("grad chi2")
plt.show()
"""
