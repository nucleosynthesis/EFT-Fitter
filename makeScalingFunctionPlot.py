import os, sys
import json
import re
from optparse import OptionParser
from collections import OrderedDict as od
from importlib import import_module
import pickle
import ROOT
import numpy as np

def get_options():
  parser = OptionParser()
  parser.add_option('--pois', dest='pois', default='params.HEL', help="Name of json file storing pois")
  parser.add_option('--functions', dest='functions', default='functions.HEL_STXS', help="Name of json file storing functions")
  parser.add_option('--inputs', dest='inputs', default='', help="Comma separated list of input files")
  parser.add_option("--translateBins", dest="translateBins", default=None, help="Translate STXS bins")
  return parser.parse_args()
(opt,args) = get_options()

# Functions for translations
def Translate(name, ndict):
    return ndict[name] if name in ndict else name
def LoadTranslations(jsonfilename):
    with open(jsonfilename) as jsonfile:
        return json.load(jsonfile)
translateBins = {} if opt.translateBins is None else LoadTranslations(opt.translateBins)

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
  _input['rho'] = _cfg.rho
  inputs.append(_input)

from tools.fitter import *

fit = fitter(pois,functions,inputs,False)

#stxs_bins = ['ttH']
stxs_bins = ['ZH_lep_PTV_0_75','ZH_lep_PTV_75_150','ZH_lep_PTV_150_250_0J','ZH_lep_PTV_150_250_GE1J','ZH_lep_PTV_GT250','ZH_lep']

scaling = od()
for stxs_bin in stxs_bins:

  scaling[stxs_bin] = od()
  for poi in pois.keys(): scaling[stxs_bin][poi] = od()

  # Quadratic
  fit.setLinearOnly(False)  
  for poi in pois.keys(): 
    scaling[stxs_bin][poi]['quad'] = od()
    c, mu = fit.scaling1D(poi,stxs_bin,npoints=1000)
    scaling[stxs_bin][poi]['quad']['c'] = c
    scaling[stxs_bin][poi]['quad']['mu'] = mu
 
  # Linear
  fit.setLinearOnly()
  for poi in pois.keys():
    scaling[stxs_bin][poi]['lin'] = od()
    c,mu = fit.scaling1D(poi,stxs_bin,npoints=1000)
    scaling[stxs_bin][poi]['lin']['c'] = c
    scaling[stxs_bin][poi]['lin']['mu'] = mu

# Mage graphs
grs = od()
for stxs_bin in stxs_bins:
  for poi in pois.keys():
    grs['%s_vs_%s_quad'%(stxs_bin,poi)] = ROOT.TGraph()
    grs['%s_vs_%s_lin'%(stxs_bin,poi)] = ROOT.TGraph()
    for i in range(len(scaling[stxs_bin][poi]['quad']['c'])): grs['%s_vs_%s_quad'%(stxs_bin,poi)].SetPoint( grs['%s_vs_%s_quad'%(stxs_bin,poi)].GetN(),scaling[stxs_bin][poi]['quad']['c'][i], scaling[stxs_bin][poi]['quad']['mu'][i] )
    for i in range(len(scaling[stxs_bin][poi]['lin']['c'])): grs['%s_vs_%s_lin'%(stxs_bin,poi)].SetPoint( grs['%s_vs_%s_lin'%(stxs_bin,poi)].GetN(),scaling[stxs_bin][poi]['lin']['c'][i], scaling[stxs_bin][poi]['lin']['mu'][i] )

# Make plot
styleMap = od()
styleMap['quad'] = {'LineWidth':3,'LineStyle':1,'MarkerSize':0}
styleMap['quad_dummy'] = {'LineWidth':3,'LineStyle':1,'MarkerSize':0}
styleMap['lin'] = {'LineWidth':2, 'LineStyle':2,'MarkerSize':0}
styleMap['lin_dummy'] = {'LineColor':12, 'LineWidth':2, 'LineStyle':2,'MarkerSize':0}
#styleMap['lin_dummy'] = {'LineColor':ROOT.kMagenta-7, 'LineWidth':2, 'LineStyle':2,'MarkerSize':0}

colorMap = od()
colorMap['ZH_lep'] = {'LineColor':ROOT.kRed-4,'MarkerColor':ROOT.kRed-4}
colorMap['ZH_lep_PTV_0_75'] = {'LineColor':ROOT.kGreen-8,'MarkerColor':ROOT.kGreen-8}
colorMap['ZH_lep_PTV_75_150'] = {'LineColor':ROOT.kGreen-7,'MarkerColor':ROOT.kGreen-7}
colorMap['ZH_lep_PTV_150_250_0J'] = {'LineColor':ROOT.kGreen+1,'MarkerColor':ROOT.kGreen+1}
colorMap['ZH_lep_PTV_150_250_GE1J'] = {'LineColor':ROOT.kGreen+3,'MarkerColor':ROOT.kGreen+3}
colorMap['ZH_lep_PTV_GT250'] = {'LineColor':ROOT.kBlack,'MarkerColor':ROOT.kBlack}
colorMap['ttH'] = {'LineColor':ROOT.kMagenta-7,'MarkerColor':ROOT.kMagenta-7}

# POI str
poi = "cWWMinuscB"
hmax = 2.5

import math
m = "%g"%math.log(1/pois[poi]['multiplier'],10)
if m == '1': m = ''
if poi == "cWWMinuscB":
  pstr_stripped = "c_{WW} #minus c_{B}"
  pstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
else:
  pstr_stripped = "c_{%s}"%poi.split("c")[-1]
  pstr = "c_{%s} x 10^{%s}"%(poi.split("c")[-1],m)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

canv = ROOT.TCanvas("canv_%s"%poi,"canv_%s"%poi,700,500)
#canv = ROOT.TCanvas("canv_%s"%poi,"canv_%s"%poi,900,500)
canv.SetBottomMargin(0.15)
canv.SetTickx()
canv.SetTicky()

prange = pois[poi]['range'][1]-pois[poi]['range'][0]

h_axes = ROOT.TH1F("haxes","",100, pois[poi]['range'][0]-0.1*prange, pois[poi]['range'][1]+0.1*prange )
h_axes.SetMaximum(hmax)
h_axes.SetMinimum(-0.2)
h_axes.SetTitle("")
h_axes.GetXaxis().SetTitle(pstr)
h_axes.GetXaxis().SetTitleSize(0.05)
h_axes.GetXaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetTitle("#mu^{i}_{prod}(%s)"%pstr_stripped)
h_axes.GetYaxis().SetTitleSize(0.05)
h_axes.GetYaxis().SetTitleOffset(0.8)
h_axes.GetYaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetLabelOffset(0.007)
h_axes.GetYaxis().CenterTitle()
h_axes.SetLineWidth(0)
h_axes.Draw()

for stxs_bin in stxs_bins:
  for k, v in colorMap[stxs_bin].iteritems():
    getattr(grs["%s_vs_%s_quad"%(stxs_bin,poi)],"Set%s"%k)(v)
    getattr(grs["%s_vs_%s_lin"%(stxs_bin,poi)],"Set%s"%k)(v)
  for k, v in styleMap['quad'].iteritems(): getattr(grs["%s_vs_%s_quad"%(stxs_bin,poi)],"Set%s"%k)(v)
  for k, v in styleMap['lin'].iteritems(): getattr(grs["%s_vs_%s_lin"%(stxs_bin,poi)],"Set%s"%k)(v)
  grs["%s_vs_%s_quad"%(stxs_bin,poi)].Draw("Same C")
  grs["%s_vs_%s_lin"%(stxs_bin,poi)].Draw("Same C")
  
# Lines
hlines = {}
yvals = [0,1]
for i in range(len(yvals)):
  yval = yvals[i]
  hlines['hline_%g'%i] = ROOT.TLine(pois[poi]['range'][0]-0.1*prange,yval,pois[poi]['range'][1]+0.1*prange,yval)
  hlines['hline_%g'%i].SetLineColorAlpha(15,0.5)
  hlines['hline_%g'%i].SetLineStyle(2)
  hlines['hline_%g'%i].SetLineWidth(1)
  hlines['hline_%g'%i].Draw("SAME")

vlines = {}
xvals = [pois[poi]['range'][0],0,pois[poi]['range'][1]]
for i in range(len(xvals)):
  xval = xvals[i]
  vlines['vline_%g'%i] = ROOT.TLine(xval,-0.2,xval,hmax)
  vlines['vline_%g'%i].SetLineColorAlpha(15,0.5)
  vlines['vline_%g'%i].SetLineStyle(2)
  vlines['vline_%g'%i].SetLineWidth(1)
  vlines['vline_%g'%i].Draw("SAME")


# Text
lat0 = ROOT.TLatex()
lat0.SetTextFont(42)
lat0.SetTextAlign(11)
lat0.SetNDC()
lat0.SetTextSize(0.045)
lat0.DrawLatex(0.1,0.92,"HEL UFO")

lat1 = ROOT.TLatex()
lat1.SetTextFont(42)
lat1.SetTextAlign(23)
lat1.SetTextSize(0.03)
xpos = pois[poi]['range'][0]-0.05*prange
lat1.DrawLatex(xpos,1.,"'#color[15]{#sigma = #sigma_{SM}}")
lat1.DrawLatex(xpos,0.,"#color[15]{#sigma = 0}")

lat2 = ROOT.TLatex()
lat2.SetTextFont(42)
lat2.SetTextAlign(23)
lat2.SetTextAngle(90)
lat2.SetTextSize(0.045)
lat2.SetTextAlign(21)
lat2.DrawLatex(pois[poi]['range'][0]-0.02*prange,0.9*hmax,"#color[15]{c_{min}}")
lat2.SetTextAlign(23)
lat2.DrawLatex(pois[poi]['range'][1]+0.01*prange,0.9*hmax,"#color[15]{c_{max}}")


# Legend

# Create dummy graph for linear
gr_lin_dummy = ROOT.TGraph()
for k,v in styleMap['lin_dummy'].iteritems(): getattr(gr_lin_dummy,"Set%s"%k)(v)

leg = ROOT.TLegend(0.55,0.22,0.8,0.48)
#leg = ROOT.TLegend(0.63,0.28,0.8,0.38)
leg.SetFillStyle(0)
leg.SetLineColor(0)
leg.SetTextSize(0.0275)
#leg.SetTextSize(0.035)
for stxs_bin in stxs_bins:  leg.AddEntry( grs["%s_vs_%s_quad"%(stxs_bin,poi)], Translate(stxs_bin,translateBins), "L")
leg.AddEntry(gr_lin_dummy,"(Lin. terms only)","L")
leg.Draw("Same")

canv.Update()
canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/ZH_lep_vs_%s.png"%poi)
canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/ZH_lep_vs_%s.pdf"%poi)
#canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/ttH_vs_%s.png"%poi)
#canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/ttH_vs_%s.pdf"%poi)
