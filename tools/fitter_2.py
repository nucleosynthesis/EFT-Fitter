# Simple python fitting object
SCIPY_MINIMIZE=False
if SCIPY_MINIMIZE : from scipy.optimize import minimize
else: import iminuit.minimize as minimize

from scipy import linalg

import array
import numpy as np
import sys
import re

from collections import OrderedDict as od

from tools.input import INPUT
from tools.helpers import *
from tools.Chi2 import GetChi2

class fitter:

  def __init__(self,parametersOfInterest,functions,inputs,doAsimov=False,theory_uncerts=None):
    self.POIS = parametersOfInterest
    self.FUNCTIONS = functions
    self.INPUTS = []
    self.doAsimov = doAsimov
    self.linearOnly = False

    self.functions = {}
    self.thuncerts = {}
    self.nps = []
    self.npindex = {}
    self.has_uncerts = False
    self.iX = 0

    self.prepareInputs(inputs)
    if theory_uncerts: self.prepareTHU(theory_uncerts)
    self.preparePOIS()
    self.preparePTerms(self.FUNCTIONS)

    self.global_min_chi2 = 0 

  def getFreePOIs(self):
    # All POIs can be included, but we can choose only some of them to float ever, important so we don't have to re-write functions 
    return list(filter(lambda x :  not self.POIS[x]["freeze"], self.POIS.keys()))
  
  def getFrozenPOIs(self):
    # All POIs can be included, but we can choose only some of them to float ever, important so we don't have to re-write functions 
    return list(filter(lambda x : self.POIS[x]["freeze"], self.POIS.keys()))

  def prepareInputs(self,inputMeasurements):
    for im in inputMeasurements:
      self.INPUTS.append( INPUT(im,self.FUNCTIONS,self.doAsimov) )

  def prepareTHU(self,thinput):
    # add the uncertainties
    print("Adding theory uncertainties")
    self.has_uncerts = True
    for x,v in thinput['TH'].items():
      self.thuncerts[x] = v
      self.nps.append(0)
      self.npindex[x]  = self.iX
      self.iX+=1
    # correlations 
    corr = []
    nbins = self.iX
    for ix in thinput['TH']:
      for jx in thinput['TH']:
        if (ix,jx) in thinput['rhoTH'].keys(): p = thinput['rhoTH'][(ix,jx)]
        elif (jx,ix) in thinput['rhoTH'].keys(): p = thinput['rhoTH'][(jx,ix)]
        else:
          if ix == jx: p = 1.
          else: 
            p = 0.
            print(" --> [WARNING] No correlation info given for (%s,%s). Assuming = 0"%(ix,jx))
        corr.append(p)
    corr = array.array('d',corr)
    corr_sq = [ corr[i:i+nbins] for i in range(0,len(corr),nbins)]
    self.tHVinv = linalg.inv(corr_sq)
    
  def preparePOIS(self):
    self.PList = []
    self.P0 = []
    self.PBounds = []
    for p, vals in self.POIS.items():
      self.PList.append(p)
      self.P0.append(vals['nominal'])
      self.PBounds.append( (vals['range'][0],vals['range'][1]) )
      if "freeze" not in vals: vals["freeze"]=0
    self.P0 = np.array( self.P0 )
    # Initially freeze all POIS: change state in minimizer function
    self.PToFitList = []

  def preparePTerms(self,functions):
    self.PTerms = od()
    for p,func in functions.items():
      self.functions[p]=func
      self.PTerms[p]=1.

  def getPOIS(self):
   pois = {}
   for i,poi in enumerate(self.PList):
      pois[poi]=self.P0[i]
   return pois

  def setPOIS(self,poiDict):
    P = od({self.PList[ip]:self.P0[ip] for ip in range(len(self.PList)) })
    P.update(poiDict)
    self.P0 = np.array([v for k,v in P.items()])
    self.evaluatePTerms()
  

  def setNuisances(self,key_vals):
    #print(key_vals)
    #print(self.npindex)
    for k,v in key_vals.items(): self.nps[k]=v # can we not always assume this runs in the right order?
    #self.nps = [v for k,v in key_vals.items()] #<- assume this is always in the right order?
    #self.nps.update(key_vals) <-  might be faster if it was in a dict?

  def evaluateTHUncertainty(self,expr):
    if not self.has_uncerts: return 0
    return self.nps[self.npindex[expr]]*self.thuncerts[expr]

  def resetNuisances(self):
    if self.has_uncerts:
     self.setNuisances({iK:0 for iK in range(len(self.nps))})
     #for iK in range(len(self.nps)): self.nps[iK] = 0
 
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
  
  def evaluatePTerms(self):
    poi_map = { p:self.P0[self.PList.index(p)] for p in self.PList }
    for p in self.PTerms.keys():
      x = self.functions[p](poi_map)
      self.PTerms[p] = x
  
  # Evaluate scaling functions for current set of POIS
  def evaluateScalingFunctions(self,term):
  # All we need to do is to look for the functio name in functions
    if term in self.PTerms.keys(): return self.PTerms[term]
    else: sys.exit("ERROR - call to evaluateScalingFunctions(%s), no known function %s "%(term,term))

  # Function to calculate chi2 for current set of POIS
  def getChi2(self,verbose=True):
    return GetChi2([],self)

  # Function to set the parameters to the global minimum
  def setGlobalMinimum(self,setParamsToNominal=False): 

    #for p in self.POIS.keys(): print(p,self.POIS[p]['nominal'])
    #print(self.P0)
    freezePOIS = self.getFrozenPOIs()
    self.minimize(freezePOIS=freezePOIS,verbose=False)
    self.global_min_chi2=self.FitResult.fun
    self.evaluatePTerms()
    print("Best-fit at global min ...")
    for i,p in enumerate(self.POIS.keys()): 
     if p in self.getFreePOIs(): print(" ",p,"%.3f"%self.P0[i])
    if setParamsToNominal :
      i=0 
      for  p,vals in self.POIS.items(): 
        vals['nominal'] = self.P0[i] 
        print("Setting as nominal , ",p,"=",vals['nominal'])
        i+=1
      print("Set nominal values of POIs to global best-fit point")
    #return np.array([self.global_min_chi2]), np.array([self.P0])

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

    if self.has_uncerts: 
      PToFit.extend([np for np in self.nps])
      PToFitBounds.extend([[-4,4] for i in self.nps])

    if SCIPY_MINIMIZE: self.FitResult = minimize(GetChi2,PToFit,args=self,bounds=PToFitBounds)
    else:  self.FitResult = minimize(GetChi2,PToFit,args=[self],bounds=PToFitBounds,options={"stra":0})
    self.setPOIS({ipoi:self.FitResult.x[ip] for ip, ipoi in enumerate(self.PToFitList)})

    # Run getChi2 with verbose messages
    if verbose: self.getChi2(verbose=True) 

    # Reset all POIs to fixed state
    self.PToFitList = []

  # Function to perform chi2 scan for single param
  def scan_fixed(self,poi,npoints=1000,reset=True):
    # Reset all pois to nominal values
    if reset: self.resetPOIS()
    print("Re-Setting to nominal values for fixed scan, ",self.P0)
    self.resetNuisances()
    
    # step one is to do a global fit to find the minimum 
    toFreezePOIs = list(self.getFreePOIs())
    toFreezePOIs.remove(poi)
    self.minimize(freezePOIS=toFreezePOIs,verbose=False) 
    nll_global = self.FitResult.fun
    
    # Loop over range of pois and calc chi2
    pvals = np.linspace( self.POIS[poi]['range'][0], self.POIS[poi]['range'][1], npoints )
    pvals
    chi2 = []
    allparams = []
    allpredictions = []

    # How can this be vectorized (issue since pois stored in class :/)
    for p in pvals:
      self.setPOIS({poi:p})
      if self.has_uncerts: 
        self.minimize(freezePOIS=self.getFreePOIs(),verbose=False) 
        chi2.append(self.FitResult.fun-nll_global) 
      else: chi2.append(self.getChi2(verbose=False)-nll_global)
      allparams.append(self.getPOIS())
      allpredictions.append(self.getPredictions())

    return pvals, np.array(chi2), allparams, allpredictions

  # Function to perform chi2 scan when profiling other parameters
  def scan_profiled(self,poi,npoints=100,freezeOtherPOIS=[],reset=True,resetEachStep=False,reverseScan=False,verbose=False):
    # Reset all pois to nominal values
    if reset: self.resetPOIS()
    self.resetNuisances()
    # Add scanned parameter into list of params to freeze
    if poi not in freezeOtherPOIS: freezeOtherPOIS.append(poi)    
    for p in self.getFrozenPOIs(): 
      if p not in freezeOtherPOIS: freezeOtherPOIS.append(p)
    # Loop over range of pois and minimize, keeping track of other parameter values
    if reverseScan: pvals = np.linspace( self.POIS[poi]['range'][1], self.POIS[poi]['range'][0], npoints )
    else: pvals = np.linspace( self.POIS[poi]['range'][0], self.POIS[poi]['range'][1], npoints )
    chi2 = []
    allpvals = []
    allparams = []
    allpredictions = []
    #print("Scanning ",poi," freezing ", freezeOtherPOIS)
    for ip,p in enumerate(pvals):
      if resetEachStep: self.resetPOIS()
      self.setPOIS({poi:p})
      #print("nuisances before minimization",self.nps)
      self.minimize(freezePOIS=freezeOtherPOIS,verbose=False)
      #print("nuisances after minimization",self.nps)
      if verbose: print(" --> [VERBOSE] Finished point (%g/%g): %s = %.3f | %s | chi2 = %.4f"%(ip,npoints,poi,p,self.getPOIStr(),self.FitResult.fun-self.global_min_chi2))
      chi2.append(self.FitResult.fun-self.global_min_chi2)
      allpvals.append(self.P0)
      allparams.append(self.getPOIS())
      allpredictions.append(self.getPredictions())

    if reverseScan: return np.flip(pvals,0), np.flip(chi2,0), np.flip(allpvals,0), allparams.reverse(), allpredictions.reverse()
    else: return pvals, np.array(chi2), np.array(allpvals), allparams, allpredictions

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

  # Function to return a map of the predictions for the current set of parameters 
  def getPredictions(self): 
    allX={}
    for _input in self.INPUTS:
      X0 = _input.X0
      X = { _input.name+_input.XList[i]: (1+self.evaluateTHUncertainty(_input.XList[i]))*self.evaluateScalingFunctions(_input.XList[i]) for i in range(_input.nbins) } 
      allX.update(X)
    return allX
