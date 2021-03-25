# Simple python fitting object
from scipy.optimize import minimize
from scipy import linalg

import array
import numpy as np
import sys
import re

from collections import OrderedDict as od
verbose=False

class fitter:

  def __init__(self,parametersOfInterest,functions,inputs,doAsimov=False):
    self.POIS = parametersOfInterest
    self.FUNCTIONS = functions
    self.INPUTS = []
    self.doAsimov = doAsimov
    self.linearOnly = False

    self.prepareInputs(inputs)
    self.preparePOIS()
    self.preparePTerms()

  def prepareInputs(self,inputMeasurements):
    for im in inputMeasurements:
      self.INPUTS.append( INPUT(im,self.FUNCTIONS,self.doAsimov) )

  def preparePOIS(self):
    self.PList = []
    self.P0 = []
    self.PBounds = []
    for p, vals in self.POIS.items():
      self.PList.append(p)
      self.P0.append(vals['nominal'])
      self.PBounds.append( (vals['range'][0],vals['range'][1]) )
    self.P0 = np.array( self.P0 )
    # Initially freeze all POIS: change state in minimizer function
    self.PToFitList = []

  def preparePTerms(self):
    # Extract all pois which enter scaling functions: add to dict
    self.PTerms = od()
    for _input in self.INPUTS:
      for i in range(_input.nbins):
        sfs = [_input.ProdScaling[i],_input.DecayScaling[i][0],_input.DecayScaling[i][1]]
        for sf in sfs:
          for term in sf:
            if term == "const": continue
            # Ignore prefactor (A/B)
            for poi in term.split("_")[1:]:
              # If poi not already in self.PTerms then add
              if poi not in self.PTerms.keys(): self.PTerms[poi] = 0.

  def setPOIS(self,poiDict):
    P = []
    for ip,ipoi in enumerate(self.PList):
      if ipoi in poiDict: P.append( poiDict[ipoi] )
      else: P.append( self.P0[ip] )
    self.P0 = np.array(P)
    # Re-evaluate the PTerms
    self.evaluatePTerms()

  def resetPOIS(self):
    P = []
    for p, vals in self.POIS.items(): P.append(vals['nominal'])
    self.P0 = np.array(P)
    # Re-evaluate the PTerms
    self.evaluatePTerms()

  def getPOIStr(self):
    pstr = "("
    for ip,p in enumerate(self.PList): pstr += " %s = %.4f,"%(p,self.P0[ip])
    pstr = "%s )"%pstr[:-1]
    return pstr 

  def setLinearOnly(self,_linearOnly=True): self.linearOnly = _linearOnly
    
  # Function evaluate all pois in scaling functions
  def evaluatePTerms(self):
    pterms = self.PTerms.keys()
    for pterm in pterms:
      x = 0
      if pterm in self.PList:
        x += self.P0[self.PList.index(pterm)]*float(self.POIS[pterm]['factor'])*float(self.POIS[pterm]['multiplier']
)
      # For rotated bases: e.g cWW-cB
      else:
        if pterm == "cWW":
          if "cWWPluscB" in self.PList:
            x += 0.5*self.P0[self.PList.index("cWWPluscB")]*float(self.POIS['cWWPluscB']['factor'])*float(self.POIS['cWWPluscB']['multiplier'])
          if "cWWMinuscB" in self.PList:
            x += 0.5*self.P0[self.PList.index("cWWMinuscB")]*float(self.POIS['cWWMinuscB']['factor'])*float(self.POIS['cWWMinuscB']['multiplier'])

        elif pterm == "cB":
          if "cWWPluscB" in self.PList:
            x += 0.5*self.P0[self.PList.index("cWWPluscB")]*float(self.POIS['cWWPluscB']['factor'])*float(self.POIS['cWWPluscB']['multiplier'])
          if "cWWMinuscB" in self.PList:
            x -= 0.5*self.P0[self.PList.index("cWWMinuscB")]*float(self.POIS['cWWMinuscB']['factor'])*float(self.POIS['cWWMinuscB']['multiplier'])

      # Update term value in pterm list
      self.PTerms[pterm] = x
  
  # Evaluate scaling functions for current set of POIS
  # this should realy be moved to the function definition themselves
  # to allow more flexibility in defining the map of params -> predictions 
  def evaluateScalingFunctions(self,terms):   
    # Calculate scaling function
    mu = terms['const']
    for term, coeff in terms.items():
      if term == "const": continue
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Linear terms
      if term[0] == "A":
        ipoi = term.split("_")[-1]
        # Add term to scaling function
        mu += coeff*self.PTerms[ipoi]

      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # Quadratic terms
      if not self.linearOnly:
        if term[0] == "B":
          # Square terms
          if len(term.split("_")) == 2:
            ipoi = term.split("_")[-1]
            # Add term to scaling function
            mu += coeff*self.PTerms[ipoi]*self.PTerms[ipoi]
             
          # Cross terms
          elif len(term.split("_")) == 3:
            ipoi = term.split("_")[1]
            jpoi = term.split("_")[2]
            # Add term to scaling function
            mu += coeff*self.PTerms[ipoi]*self.PTerms[jpoi]

    # Confine mu to be positive: rounding issue with prefactors
    #return max(0.,mu)

    return mu


  # Function to calculate chi2 for current set of POIS
  def getChi2(self,verbose=True):
    return GetChi2([],self)


  # Minimizer function
  def minimize(self,freezePOIS=[],verbose=True):

    # Define pois to profile
    self.PToFitList, PToFit, PToFitBounds = [], [], []
    for ip, ipoi in enumerate(self.PList):
      if ipoi in freezePOIS: continue
      else: 
        self.PToFitList.append(ipoi)
        PToFit.append(self.P0[ip])
        PToFitBounds.append(self.PBounds[ip])

    # Do minimisation
    #print ("need to fit",PToFit)
    self.FitResult = minimize(GetChi2,PToFit,args=self,bounds=PToFitBounds,method='Nelder-mead')
    #print(self.FitResult)
    # Set POI values for those profiled
    for ip, ipoi in enumerate(self.PToFitList): self.setPOIS({ipoi:self.FitResult.x[ip]})

    # Run getChi2 with verbose messages
    if verbose: self.getChi2(verbose=True) 

    # Reset all POIs to fixed state
    self.PToFitList = []


  # Function to perform chi2 scan for single param
  def scan_fixed(self,poi,npoints=1000,reset=True):
    # Reset all pois to nominal values
    if reset: self.resetPOIS()
    # Loop over range of pois and calc chi2
    pvals = np.linspace( self.POIS[poi]['range'][0], self.POIS[poi]['range'][1], npoints )
    chi2 = []
    for p in pvals:
      self.setPOIS({poi:p})
      chi2.append(self.getChi2(verbose=False))
    return pvals, np.array(chi2)

  # Function to perform chi2 scan when profiling other parameters
  def scan_profiled(self,poi,npoints=100,freezeOtherPOIS=[],reset=True,resetEachStep=False,reverseScan=False,verbose=False):
    # Reset all pois to nominal values
    if reset: self.resetPOIS()

    # Add scanned parameter into list of params to freeze
    if poi not in freezeOtherPOIS: freezeOtherPOIS.append(poi)    

    # Loop over range of pois and minimize, keeping track of other parameter values
    if reverseScan: pvals = np.linspace( self.POIS[poi]['range'][1], self.POIS[poi]['range'][0], npoints )
    else: pvals = np.linspace( self.POIS[poi]['range'][0], self.POIS[poi]['range'][1], npoints )
    chi2 = []
    allpvals = []
    for ip,p in enumerate(pvals):
      if resetEachStep: self.resetPOIS()
      self.setPOIS({poi:p})
      self.minimize(freezePOIS=freezeOtherPOIS,verbose=False)
      if verbose: print(" --> [VERBOSE] Finished point (%g/%g): %s = %.3f | %s | chi2 = %.4f"%(ip,npoints,poi,p,self.getPOIStr(),self.FitResult.fun))
      chi2.append(self.FitResult.fun)
      allpvals.append(self.P0)

    # Drop first element: minimizer not yet steady
    pvals = pvals[1:]
    chi2 = np.array(chi2)[1:]
    allpvals = np.array(allpvals)[1:]

    if reverseScan: return np.flip(pvals,0), np.flip(chi2,0), np.flip(allpvals,0)
    else: return pvals, np.array(chi2), np.array(allpvals)

  # Function to extract value of scling fnct for poi
  def scaling1D(self,poi,func,npoints=1000,reset=True):
    # Reset all pois to nominal values
    if reset: self.resetPOIS()
    pvals = np.linspace( self.POIS[poi]['range'][0], self.POIS[poi]['range'][1], npoints )
    mu = []
    for p in pvals:
      self.setPOIS({poi:p})
      mu.append(self.evaluateScalingFunctions(extractTerms(self.FUNCTIONS[func])))
    mu = np.array(mu)
    return pvals, mu

  # Function to extract value of scling fnct for poi
  def scaling2D(self,xpoi,ypoi,func,npoints=[100,100],reset=True):
    # Reset all pois to nominal values
    if reset: self.resetPOIS()
    xpvals = np.linspace( self.POIS[xpoi]['range'][0], self.POIS[xpoi]['range'][1], npoints[0] )
    ypvals = np.linspace( self.POIS[ypoi]['range'][0], self.POIS[ypoi]['range'][1], npoints[1] )
    x = []
    y = []
    mu = []
    for xp in xpvals:
      for yp in ypvals:
        self.setPOIS({xpoi:xp,ypoi:yp})
        x.append(xp)
        y.append(yp)
        mu.append(self.evaluateScalingFunctions(extractTerms(self.FUNCTIONS[func])))
    return np.array([x,y]).transpose(), np.array(mu)


  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chi2: function to minimize
#  * Takes as input array of POIS to fit, if list is empty then will not reset
def GetChi2(P,FIT):
  chi2 = 0

  # Set POIS for minimizer
  if len(P) != 0:
    for ip,ipoi in enumerate(FIT.PToFitList): FIT.setPOIS({ipoi:P[ip]})

  # Calculate chi2 terms
  for _input in FIT.INPUTS:
    X0 = _input.X0
    X = np.asarray( [ FIT.evaluateScalingFunctions(_input.ProdScaling[i])*(FIT.evaluateScalingFunctions(_input.DecayScaling[i][0])/FIT.evaluateScalingFunctions(_input.DecayScaling[i][1])) for i in range(_input.nbins)] )
    
    if _input.type == "spline":
          #print("to fit -> ",FIT.PToFitList)
          #for ip,ipoi in enumerate(FIT.PToFitList): 
            #print({ipoi:P[ip]})
          #print("evaluated -> ",{"%s"%_input.XList[j]:X[j] for j in range(len(X))})
          chi2 = X0.evaluate({"%s"%_input.XList[j]:X[j] for j in range(len(X))})
          #print("...ch2  = ",chi2)
          return chi2
    else: 
          xarr = X-X0
          xarrT = xarr.T
          chi2 += xarrT.dot( _input.Vinv.dot( xarr ) )

  if verbose:
    print(" ----------------------------------------------------------------------------")
    print(" --> Scaling functions: Linear" if FIT.linearOnly else " --> Scaling functions: Quadratic")
    print("")
    print(" --> Parameters of interest:")
    for ip,ipoi in enumerate(FIT.PList): 
     if ipoi in FIT.PToFitList: print("   * %s = %.3f"%(ipoi,FIT.P0[ip]))
     else: print("   * %s = %.3f (fixed)"%(ipoi,FIT.P0[ip]))

    for _input in FIT.INPUTS:
      print("\n --> Inputs: %s"%_input.name)
      if _input.type == "spline" : continue
      X0 = _input.X0
      xvals_list = [ FIT.evaluateScalingFunctions(_input.ProdScaling[i])*(FIT.evaluateScalingFunctions(_input.DecayScaling[i][0])/FIT.evaluateScalingFunctions(_input.DecayScaling[i][1])) for i in range(_input.nbins)]
      X = np.asarray(xvals_list)
      for ix in range(_input.nbins):
        print("   * %-50s : X0 = %.6f, X(p) = %.6f"%(_input.XList[ix],X0[ix],X[ix]))
      print("   * Vij = ")
      for i in range(len(_input.Vinv[0])):
        vstr =  "           ("
        for j in range(len(_input.Vinv[0])): vstr += " %-5.2f "%_input.Vinv[i][j]
        vstr += ")"
        print(vstr)

    print("\n --> Result: chi2 = %.8f"%chi2)
      
  return chi2
 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract all required info for input measurement
class INPUT:

  def __init__(self,inputMeasurement,functions,doAsimov):

    # Name
    self.name = inputMeasurement['name']

    # X0
    self.X0 = []
    self.XList = []
    
    if (len(inputMeasurement['X'].keys()))>1:
      self.type = "measurement"
      for x, vals in inputMeasurement['X'].items(): 
        if doAsimov: self.X0.append(1.)
        else: self.X0.append(float(vals['bestfit']))
        self.XList.append(x)
      self.X0 = np.asarray(self.X0)
    else:
      self.type = "spline"
      self.X0 = inputMeasurement['X']["likelihood"]['likelihood']
      self.XList = self.X0.getParameters()
      inputMeasurement['X']['likelihood']["merged"] = False
      for k in self.XList  : 
        inputMeasurement['X'][k] = {"likelihood":1} # could use this better
        inputMeasurement['X'][k]["merged"] = False

    # Scaling functions
    self.ProdScaling = []
    self.DecayScaling = []
    for x, vals in inputMeasurement["X"].items(): 
      if x == "likelihood": continue
      
      # first check if the function is already provided 
      if x in functions.keys(): 
       self.ProdScaling.append( extractTerms(functions[x]) )
       self.DecayScaling.append( [extractTerms("1."),extractTerms(functions['tot'])] )
       continue
      # Then appeal to prod*decay

      # Extract production and decay strings
      prod = "_".join(x.split("_")[:-1])
      dec = x.split("_")[-1]

      # For merged STXS bins
      if vals['merged']:
        terms = od()
        for stxs_bin,frac in vals['STXS_fractions'].items():
          _terms = extractTerms(functions[stxs_bin],multiplier=frac)
          for k,v in _terms.items():
            if k in terms: terms[k] = terms[k]+v
            else: terms[k] = v
        self.ProdScaling.append(terms)
        # Add function to list
        functions[prod] = termsToFunction(terms)

      # For unmerged STXS bin/process
      else: 
        self.ProdScaling.append( extractTerms(functions[prod]) )

      # Add list of partial width and total width scaling
      self.DecayScaling.append( [extractTerms(functions[dec]),extractTerms(functions['tot'])] )

    if self.type=="spline": 
      nbins = len( self.XList )
      self.nbins = nbins
      return

    # Error matrix
    corr = []
    for ix in inputMeasurement['X']:
      for jx in inputMeasurement['X']:
        if (ix,jx) in inputMeasurement['rho'].keys(): p = inputMeasurement['rho'][(ix,jx)]
        elif (jx,ix) in inputMeasurement['rho'].keys(): p = inputMeasurement['rho'][(jx,ix)]
        else:
          if ix == jx: p = 1.
          else: 
            p = 0.
            print(" --> [WARNING] No correlation info given for (%s,%s). Assuming = 0"%(ix,jx))
        corr.append(p)
    corr = array.array('d',corr)

    # Extract symmetrized errors
    # FIXME: add something more appropriate for observed!
    err = []
    for x,vals in inputMeasurement['X'].items():
      if doAsimov: err.append( 0.5*(vals['Up01SigmaExp']+vals['Down01SigmaExp']) )
      else: err.append( 0.5*(vals['Up01Sigma']+vals['Down01Sigma']) )
    
    # Do some squarification and inverting
    nbins = len( inputMeasurement['X'].keys() )
    corr_sq = [ corr[i:i+nbins] for i in range(0,len(corr),nbins)]
    cov_sq = [ [ corr_sq[i][j]*err[i]*err[j] for i in range(nbins)] for j in range(nbins) ]
    self.err_mat = np.array(cov_sq)
    self.Vinv = linalg.inv(self.err_mat)
    self.nbins = nbins
    
# Define functions
def extractTerms(_func,multiplier=1):
  terms = od()
  terms['const'] = 0
  f = re.sub(" ","",_func)
  f = re.sub("-","+-",f)
  for t in f.split("+"):
    c = t.split("*")
    if len(c) == 1: terms['const'] = float(c[0])*multiplier
    elif len(c) == 2: terms['A_%s'%c[1]] = float(c[0])*multiplier
    else:
      if c[1] == c[2]: terms['B_%s'%c[1]] = float(c[0])*multiplier
      else: terms['B_%s_%s'%(c[1],c[2])] = float(c[0])*multiplier
  return terms

# Terms to function
def termsToFunction(_terms):
  f = ""
  for k,v in _terms.items():
    if k == "const": f += "%.4f"%v
    elif k[0] == "A": f += "+%.4f*%s"%(v,k.split("_")[-1])
    elif k[0] == "B": 
      if len(k.split("_"))==2: f += "+%.4f*%s*%s"%(v,k.split("_")[-1],k.split("_")[-1])
      else: f += "+%.4f*%s*%s"%(v,k.split("_")[1],k.split("_")[2])
  f = re.sub("\+-","-",f)
  return f

# Function for printing matrix
def printMatrix(_mat):
  for i in range(len(_mat[0])):
    vstr =  "("
    for j in range(len(_mat[0])): vstr += " %-5.2f "%_mat[i][j]
    vstr += ")"
    print(vstr)
