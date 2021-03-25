from inputs.test_lhcscan_2D import *

import numpy
import matplotlib.pyplot as plt
from pylab import meshgrid, contour, show, imshow, clabel

x = numpy.arange(0,4,0.1)
y = numpy.arange(0,10,0.1)
X,Y = meshgrid(x,y)

@numpy.vectorize
def simple_eval(xx,yy):
  return chi2.evaluate({"r_ggH":xx,"r_qqH":yy})


Z = simple_eval(X,Y)
#im   = imshow(Z,aspect = 1)
cset = contour(X,Y,Z,levels=100)
clabel(cset,inline=True,fmt='%1.1f',fontsize=10)
plt.xlabel("$\mu_{F}$")
plt.ylabel("$\mu_{V}$")
#colorbar(im)


from scipy.optimize import minimize

cGvals = numpy.arange(0,4,0.1)
def getchi2(ca,args):
  cg  = args[0]
  return simple_eval(cg,ca)

init = [5.]
cAvals = []
cAvals_n = []
cAvals_2 = []
for cg in cGvals:
 res =  minimize(getchi2,init,args=[cg])
 mincax2 = 9999.
 camin = 5.
 for ca in numpy.arange(0,10,0.1):
   ch2 = getchi2(ca,[cg])
   if ch2<mincax2: 
    camin = ca
    mincax2 = ch2
 cAvals.append(res.x[0])
 cAvals_n.append(camin)
 res2 =  minimize(getchi2,init,args=[cg],method='Nelder-mead')
 cAvals_2.append(res2.x[0])

plt.plot(cGvals,cAvals)
plt.plot(cGvals,cAvals_n)
plt.plot(cGvals,cAvals_2)

show()


