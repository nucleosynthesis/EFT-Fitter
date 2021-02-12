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
  parser.add_option('--poi', dest='poi', default='cG', help="POI to plot")
  parser.add_option('--npoints', dest='npoints', default=1000, type='int', help="Number of points")
  parser.add_option('--hmax', dest='hmax', default=2.5, type='float', help="Histogram max")
  parser.add_option('--hmin', dest='hmin', default=-0.2, type='float', help="Histogram min")
  parser.add_option("--translateChannels", dest="translateChannels", default=None, help="Translate STXS bins")
  parser.add_option("--leg-pos", dest="leg_pos", default="bottom_right", help="Legend position")
  return parser.parse_args()
(opt,args) = get_options()

# Functions for translations
def Translate(name, ndict):
    return ndict[name] if name in ndict else name
def LoadTranslations(jsonfilename):
    with open(jsonfilename) as jsonfile:
        return json.load(jsonfile)
translateChannels = {} if opt.translateChannels is None else LoadTranslations(opt.translateChannels)

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

channels = ['hgg','hzz','hww','htt','hbb']

scaling = od()
# Extract total width scaling
scaling['tot'] = od()

# Quadratic
fit.setLinearOnly(False)
scaling['tot']['quad'] = od()
c,mu = fit.scaling1D(opt.poi,'tot',npoints=opt.npoints)
scaling['tot']['quad']['c'] = c
scaling['tot']['quad']['mu'] = mu

# Linear
fit.setLinearOnly()
scaling['tot']['lin'] = od()
c,mu = fit.scaling1D(opt.poi,'tot',npoints=opt.npoints)
scaling['tot']['lin']['c'] = c
scaling['tot']['lin']['mu'] = mu

# Extract BR and partial width scaling for channels
ch_unaffected = []
for ch in channels:

  scaling[ch] = od()
  scaling[ch]['partial'] = od()
  scaling[ch]['br'] = od()

  # Quadratic
  fit.setLinearOnly(False)  
  scaling[ch]['partial']['quad'] = od()
  scaling[ch]['br']['quad'] = od()
  c, mu = fit.scaling1D(opt.poi,ch,npoints=opt.npoints)
  # If partial width always 1 then add to list
  if (mu == np.ones(len(c))).sum() == opt.npoints: ch_unaffected.append(ch)
  scaling[ch]['partial']['quad']['c'] = c
  scaling[ch]['partial']['quad']['mu'] = mu
  scaling[ch]['br']['quad']['c'] = c
  scaling[ch]['br']['quad']['mu'] = mu/scaling['tot']['quad']['mu']

   # Linear
  fit.setLinearOnly()  
  scaling[ch]['partial']['lin'] = od()
  scaling[ch]['br']['lin'] = od()
  c, mu = fit.scaling1D(opt.poi,ch,npoints=opt.npoints)
  scaling[ch]['partial']['lin']['c'] = c
  scaling[ch]['partial']['lin']['mu'] = mu
  scaling[ch]['br']['lin']['c'] = c
  scaling[ch]['br']['lin']['mu'] = mu/scaling['tot']['lin']['mu']
 
# Make graphs
shifter =0.005
grs = od()
grs["tot_vs_%s_quad"%opt.poi] = ROOT.TGraph()
grs["tot_vs_%s_lin"%opt.poi] = ROOT.TGraph()
for i in range(len(scaling['tot']['quad']['c'])): grs['tot_vs_%s_quad'%opt.poi].SetPoint( grs['tot_vs_%s_quad'%opt.poi].GetN(),scaling['tot']['quad']['c'][i], scaling['tot']['quad']['mu'][i] )
for i in range(len(scaling['tot']['lin']['c'])): grs['tot_vs_%s_lin'%opt.poi].SetPoint( grs['tot_vs_%s_lin'%opt.poi].GetN(),scaling['tot']['lin']['c'][i], scaling['tot']['lin']['mu'][i] )
for ch in channels:
  for mode in ['partial','br']:
    grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)] = ROOT.TGraph()
    for i in range(len(scaling[ch][mode]['quad']['c'])): 
      mu = scaling[ch][mode]['quad']['mu'][i]
      if ch in ch_unaffected: 
        # Find index
        chidx = ch_unaffected.index(ch)
        if chidx % 2 == 0: mu = mu+shifter*(int(0.5*chidx)+1)
        else: mu = mu-shifter*(int(0.5*chidx)+1)
      grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)].SetPoint( grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)].GetN(),scaling[ch][mode]['quad']['c'][i], mu )
    grs["%s_vs_%s_%s_lin"%(ch,opt.poi,mode)] = ROOT.TGraph()
    for i in range(len(scaling[ch][mode]['lin']['c'])): 
      mu = scaling[ch][mode]['lin']['mu'][i]
      if ch in ch_unaffected: 
        # Find index
        chidx = ch_unaffected.index(ch)
        if chidx % 2 == 0: mu = mu+shifter*(int(0.5*chidx)+1)
        else: mu = mu-shifter*(int(0.5*chidx)+1)
      grs["%s_vs_%s_%s_lin"%(ch,opt.poi,mode)].SetPoint( grs["%s_vs_%s_%s_lin"%(ch,opt.poi,mode)].GetN(),scaling[ch][mode]['lin']['c'][i], mu )

# Make plot
styleMap = od()
styleMap['br_quad'] = {'LineWidth':2,'LineStyle':1,'MarkerSize':0}
styleMap['partial_quad'] = {'LineWidth':2,'LineStyle':2,'MarkerSize':0}
styleMap['partial_dummy'] = {'LineColor':12,'LineWidth':2,'LineStyle':2,'MarkerSize':0}
styleMap['tot_quad'] = {'LineColor':1,'LineWidth':2,'LineStyle':1,'MarkerSize':0}

colorMap = od()
colorMap['hgg'] = {'LineColor':ROOT.kRed-4,'MarkerColor':ROOT.kRed-4}
colorMap['hzz'] = {'LineColor':ROOT.kGreen-3,'MarkerColor':ROOT.kGreen-3}
colorMap['hww'] = {'LineColor':ROOT.kOrange-3,'MarkerColor':ROOT.kOrange-3}
colorMap['htt'] = {'LineColor':ROOT.kCyan-7,'MarkerColor':ROOT.kCyan-7}
colorMap['hbb'] = {'LineColor':ROOT.kViolet+6,'MarkerColor':ROOT.kViolet+6}

# POI str
hmax = opt.hmax
hmin = opt.hmin

import math
m = "%g"%math.log(1/pois[opt.poi]['multiplier'],10)
if m == '1': m = ''
if opt.poi == "cWWMinuscB":
  pstr_stripped = "c_{WW} #minus c_{B}"
  pstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
else:
  pstr_stripped = "c_{%s}"%opt.poi.split("c")[-1]
  pstr = "c_{%s} x 10^{%s}"%(opt.poi.split("c")[-1],m)

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

canv = ROOT.TCanvas("canv_%s"%opt.poi,"canv_%s"%opt.poi,900,500)
canv.SetBottomMargin(0.15)
canv.SetTickx()
canv.SetTicky()

prange = pois[opt.poi]['range'][1]-pois[opt.poi]['range'][0]

h_axes = ROOT.TH1F("haxes","",100, pois[opt.poi]['range'][0]-0.1*prange, pois[opt.poi]['range'][1]+0.1*prange )
h_axes.SetMaximum(hmax)
h_axes.SetMinimum(hmin)
h_axes.SetTitle("")
h_axes.GetXaxis().SetTitle(pstr)
h_axes.GetXaxis().SetTitleSize(0.05)
h_axes.GetXaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetTitle("#mu^{f}_{decay}(%s)"%pstr_stripped)
h_axes.GetYaxis().SetTitleSize(0.05)
h_axes.GetYaxis().SetTitleOffset(0.8)
h_axes.GetYaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetLabelOffset(0.007)
h_axes.GetYaxis().CenterTitle()
h_axes.SetLineWidth(0)
h_axes.Draw()


# Loop over channels and plot partial width and br
for ch in channels:
  for mode in ['partial','br']:
    for k,v in styleMap['%s_quad'%mode].iteritems(): getattr(grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)],"Set%s"%k)(v)
    if ch in ch_unaffected: 
      getattr(grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)],"SetLineWidth")(1)
    for k,v in colorMap[ch].iteritems(): getattr(grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)],"Set%s"%k)(v)
    grs["%s_vs_%s_%s_quad"%(ch,opt.poi,mode)].Draw("Same C")

# Plot total width 
for k,v in styleMap['tot_quad'].iteritems(): getattr(grs["tot_vs_%s_quad"%opt.poi],"Set%s"%k)(v)
grs["tot_vs_%s_quad"%opt.poi].Draw("Same C")


# Lines
hlines = {}
yvals = [0,1]
for i in range(len(yvals)):
  yval = yvals[i]
  if( yval > hmax )|( yval < hmin ): continue
  hlines['hline_%g'%i] = ROOT.TLine(pois[opt.poi]['range'][0]-0.1*prange,yval,pois[opt.poi]['range'][1]+0.1*prange,yval)
  hlines['hline_%g'%i].SetLineColorAlpha(15,0.5)
  hlines['hline_%g'%i].SetLineStyle(2)
  hlines['hline_%g'%i].SetLineWidth(1)
  hlines['hline_%g'%i].Draw("SAME")

vlines = {}
xvals = [pois[opt.poi]['range'][0],0,pois[opt.poi]['range'][1]]
for i in range(len(xvals)):
  xval = xvals[i]
  vlines['vline_%g'%i] = ROOT.TLine(xval,hmin,xval,hmax)
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
xpos = pois[opt.poi]['range'][0]-0.04*prange
if hmin < -0.01: lat1.DrawLatex(xpos,-0.01,"#color[15]{#bf{#it{#Beta}} = 0}")
if hmax > 1.: lat1.DrawLatex(xpos,0.99,"#color[15]{#bf{#it{#Beta}} = #bf{#it{#Beta}}_{SM}}")

lat2 = ROOT.TLatex()
lat2.SetTextFont(42)
lat2.SetTextAlign(23)
lat2.SetTextAngle(90)
lat2.SetTextSize(0.045)
lat2.SetTextAlign(21)
lat2.DrawLatex(pois[opt.poi]['range'][0]-0.02*prange,0.9*(hmax-hmin)+hmin,"#color[15]{c_{min}}")
lat2.SetTextAlign(23)
lat2.DrawLatex(pois[opt.poi]['range'][1]+0.01*prange,0.9*(hmax-hmin)+hmin,"#color[15]{c_{max}}")


# Legend

# Create dummy graph for linear
gr_dummy = ROOT.TGraph()
for k,v in styleMap['partial_dummy'].iteritems(): getattr(gr_dummy,"Set%s"%k)(v)

if opt.leg_pos == "bottom_right": leg = ROOT.TLegend(0.7,0.21,0.8,0.47)
elif opt.leg_pos == "top_right": leg = ROOT.TLegend(0.7,0.59,0.8,0.85)
elif opt.leg_pos == "top_rightv2": leg = ROOT.TLegend(0.73,0.59,0.83,0.85)
elif opt.leg_pos == "top_left": leg = ROOT.TLegend(0.2,0.59,0.3,0.85)
else: leg = ROOT.TLegend(0.7,0.21,0.8,0.47)
leg.SetFillStyle(0)
leg.SetLineColor(0)
leg.SetTextSize(0.03)
leg.AddEntry(grs["tot_vs_%s_quad"%opt.poi],"#Gamma^{H}","L")
for ch in channels:  leg.AddEntry( grs["%s_vs_%s_br_quad"%(ch,opt.poi)], Translate(ch,translateChannels), "L")
leg.AddEntry(gr_dummy,"#Gamma^{f}","L")
leg.Draw("Same")

canv.Update()
canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/decay_vs_%s.png"%opt.poi)
canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/decay_vs_%s.pdf"%opt.poi)
