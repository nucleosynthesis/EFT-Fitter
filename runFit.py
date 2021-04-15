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
  parser.add_option('--output', dest='outputstr', default='', help="Identifier string for output results")
  parser.add_option('--functions', dest='functions', default='functions.HEL_STXS', help="Name of json file storing functions")
  parser.add_option('--inputs', dest='inputs', default='', help="Comma separated list of input files")
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

# Load functions
functions = import_module(opt.functions).functions

# Load input measurements
inputs = []
for i in opt.inputs.split(","):
  _cfg = import_module(i)
  _input = od()
  _input['name'] = _cfg.name
  _input['X'] = _cfg.X
  if "rho" in _input.keys(): _input['rho'] = _cfg.rho
  inputs.append(_input)

if len(opt.theory_uncerts):
  _cfg = import_module(opt.theory_uncerts)
  _th_input  = od()
  _th_input['TH'] = _cfg.TH
  if hasattr(_cfg,"rhoTH"):  _th_input['rhoTH'] = _cfg.rhoTH
  opt.theory_uncerts = _th_input
from tools.fitter import *

fit = fitter(pois,functions,inputs,opt.doAsimov,opt.theory_uncerts)

# Perform scans
results = od()

fit.setGlobalMinimum(opt.setParamsToNominal)
for poi in pois.keys():

  print(" --> Running fits for: %s"%poi)
  results[poi] = od()

  # Quadratic
  fit.setLinearOnly(False)

  # Fixed scan
  print("    * Quadratic: fixed")
  pvals_f, chi2_f = fit.scan_fixed(poi,npoints=100)
  results[poi]["fixed"] = od()
  results[poi]["fixed"]['pvals'] = pvals_f
  results[poi]["fixed"]['chi2'] = chi2_f
  results[poi]["fixed"]['dchi2'] = chi2_f-chi2_f.min()

  # Profiled scan (full)
  print("    * Quadratic: profiled")
  pvals_p, chi2_p, all_pvals_p = fit.scan_profiled(poi,npoints=50,freezeOtherPOIS=[],resetEachStep=opt.doReset,reverseScan=opt.doFlip,verbose=True)
  results[poi]["profiled"] = od()
  results[poi]["profiled"]['pvals'] = pvals_p
  results[poi]["profiled"]['chi2'] = chi2_p
  results[poi]["profiled"]['allpvals'] = all_pvals_p
  results[poi]["profiled"]['dchi2'] = chi2_p-chi2_p.min()

  # Linear
  if not opt.doLinear: continue
  fit.setLinearOnly(True)

  # Fixed scan
  print("    * Linear: fixed")
  pvals_f_lin, chi2_f_lin = fit.scan_fixed(poi,npoints=100)
  results[poi]["fixed_linear"] = od()
  results[poi]["fixed_linear"]['pvals'] = pvals_f_lin
  results[poi]["fixed_linear"]['chi2'] = chi2_f_lin
  results[poi]["fixed_linear"]['dchi2'] = chi2_f_lin-chi2_f_lin.min()

  # Profiled scan (full)
  print("    * Linear: profiled")
  pvals_p_lin, chi2_p_lin, all_pvals_p_lin = fit.scan_profiled(poi,npoints=50,freezeOtherPOIS=[],resetEachStep=opt.doReset,reverseScan=opt.doFlip,verbose=True)
  results[poi]["profiled_linear"] = od()
  results[poi]["profiled_linear"]['pvals'] = pvals_p_lin
  results[poi]["profiled_linear"]['chi2'] = chi2_p_lin
  results[poi]["profiled_linear"]['allpvals'] = all_pvals_p_lin
  results[poi]["profiled_linear"]['dchi2'] = chi2_p_lin-chi2_p_lin.min() 
  

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

