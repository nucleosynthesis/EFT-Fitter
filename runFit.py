import os, sys
import json
import re
from optparse import OptionParser
from collections import OrderedDict as od

from importlib import import_module
import pickle

def get_options():
  parser = OptionParser()
  parser.add_option('--pois', dest='pois', default='params.HEL', help="Name of json file storing pois")
  parser.add_option('--scanpois', dest='scanpois', default='', help="comma separated list of pois to scan (if empty, default is all pois)")
  parser.add_option('--output', dest='outputstr', default='', help="Identifier string for output results")
  parser.add_option('--functions', dest='functions', default='functions.HEL_STXS', help="Name of json file storing functions")
  parser.add_option('--inputs', dest='inputs', default='', help="Comma separated list of input files")
  parser.add_option('--npoints', dest='npoints', default=20,type=int, help="number of points in the scan")
  parser.add_option('--theory_uncert', dest='theory_uncerts', default='', help="config for theory uncertainties")
  parser.add_option('--doAsimov', dest='doAsimov', default=False, action="store_true", help="Do asimov fit (i.e. set all best-fit to nominal)")
  parser.add_option('--doReset', dest='doReset', default=False, action="store_true", help="Reset poi values each step in profiled scan")
  parser.add_option('--doFlip', dest='doFlip', default=False, action="store_true", help="Start scan from max val of poi")
  parser.add_option('--doLinear', dest='doLinear', default=False, action="store_true", help="Also run the scan using linear terms of functions (defined in --functions) -- only appropriate for EFT models with SM+linear+BSM terms)")
  parser.add_option('--setParamsToNominal', dest='setParamsToNominal', default=False, action="store_true", help="Set nominal values of the POIs to those at the global minimum")
  return parser.parse_args()
(opt,args) = get_options()

# Load parameters of interest
pois = import_module(opt.pois).pois
if len(opt.scanpois) : opt.scanpois = opt.scanpois.split(",")
else: opt.scanpois = []
# Load functions
functions = import_module(opt.functions).functions
try:  
  grad_functions = import_module(opt.functions).grad_functions
except:
  print ("No gradient functions found in ",opt.functions)
  grad_functions=functions

# Load input measurements
inputs = []
for i in opt.inputs.split(","):
  _cfg = import_module(i)
  _input = od()
  _input['name'] = _cfg.name
  _input['X'] = _cfg.X
  if hasattr(_cfg,"rho"): _input['rho'] = _cfg.rho
  inputs.append(_input)

if len(opt.theory_uncerts):
  _cfg = import_module(opt.theory_uncerts)
  _th_input  = od()
  _th_input['TH'] = _cfg.TH
  if hasattr(_cfg,"rhoTH"):  _th_input['rhoTH'] = _cfg.rhoTH
  opt.theory_uncerts = _th_input
from tools.fitter_2 import *

fit = fitter(pois,functions,grad_functions,inputs,opt.doAsimov,opt.theory_uncerts)
# Perform scans
results = od()

fit.setGlobalMinimum(opt.setParamsToNominal)

for poi in fit.getFreePOIs():

  if (len(opt.scanpois)>0) and (poi not in opt.scanpois): continue 
  print(" --> Running fits for: %s"%poi)
  results[poi] = od()

  # Quadratic
  fit.setLinearOnly(False)

  # Profiled scan (full)
  print("    * Quadratic: profiled")
  result = fit.scan_profiled(poi,npoints=opt.npoints,freezeOtherPOIS=[],resetEachStep=opt.doReset,reverseScan=opt.doFlip,verbose=True)
  
  results[poi]["profiled"] = od()
  results[poi]["profiled"]['pvals'] = result.pvals
  results[poi]["profiled"]['chi2'] = result.chi2
  results[poi]["profiled"]['allpvals'] = result.allpvals
  results[poi]["profiled"]['dchi2'] = result.chi2-result.chi2.min()
  
  results[poi]["profiled"]['predictions']   = {}
  for pred_l in result.allpredictions[0].keys() : 
    results[poi]["profiled"]['predictions'][pred_l] = [ result.allpredictions[i][pred_l] for i in range(len(result.allpredictions)) ]
  
  results[poi]["profiled"]['otherpoi']   = {}
  for other_poi_l in result.allparams[0].keys() : 
    results[poi]["profiled"]['otherpoi'][other_poi_l] = [ result.allparams[i][other_poi_l] for i in range(len(result.allparams)) ]

extStr = opt.outputstr
if opt.doAsimov: extStr += "_asimov"
else: extStr += "_observed"
if opt.doReset: extStr += "_reset"
if opt.doFlip: extStr += "_flip"

with open("results%s.pkl"%extStr,"wb") as fpkl: pickle.dump(results,fpkl)
#from scipy import interpolate
#import matplotlib.pyplot as plt
#f = interpolate.interp1d(cg,dchi2)
#cg_new = np.linspace(-8,10,1000)
#dchi2_new = f(cg_new)
#plt.plot(cg,dchi2,'o',cg_new,dchi2_new,'-')

