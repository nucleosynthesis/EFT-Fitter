import os, sys
import json
import re
from optparse import OptionParser
from collections import OrderedDict as od
from scipy.interpolate import griddata
from tools.shanePalette import set_color_palette

from importlib import import_module
import pickle
import ROOT
import numpy as np

def get_options():
  parser = OptionParser()
  parser.add_option('--pois', dest='pois', default='params.HEL', help="Name of json file storing pois")
  parser.add_option('--functions', dest='functions', default='functions.HEL_STXS', help="Name of json file storing functions")
  parser.add_option('--inputs', dest='inputs', default='', help="Comma separated list of input files")
  parser.add_option('--doLinear', dest='doLinear', default=False, action="store_true", help="Linear only")
  parser.add_option('--doStripCross', dest='doStripCross', default=False, action="store_true", help="Add contours for dropping Bjk terms (j!=k)")
  parser.add_option('--doProd', dest='doProd', default=False, action="store_true", help="Add contours for production scaling")
  parser.add_option('--doDec', dest='doDec', default=False, action="store_true", help="Add contours for decay scaling")
  parser.add_option('--xpoi', dest='xpoi', default='cG', help="x POI to plot")
  parser.add_option('--ypoi', dest='ypoi', default='cA', help="y POI to plot")
  parser.add_option('--xnpoints', dest='xnpoints', default=100, type='int', help="Number of points")
  parser.add_option('--ynpoints', dest='ynpoints', default=100, type='int', help="Number of points")
  parser.add_option('--nInterpolatePoints', dest='nInterpolatePoints', default=1000, type='int', help="Number of points for interpolation")
  parser.add_option('--nBins', dest='nBins', default=200, type='int', help="Number of bins in plot")
  parser.add_option('--proc', dest='proc', default="", help="STXS bin")
  parser.add_option('--dec', dest='dec', default="", help="Decay channel")
  parser.add_option("--translateBins", dest="translateBins", default=None, help="Translate STXS bins")
  parser.add_option("--translateChannels", dest="translateChannels", default=None, help="Translate decay channels")
  return parser.parse_args()
(opt,args) = get_options()

# Set color palette
#set_color_palette('jonno_flip')
set_color_palette('jonno_flip_qqh')
ROOT.gStyle.SetLineStyleString(2,"5 5")


# Functions for translations
def Translate(name, ndict):
    return ndict[name] if name in ndict else name
def LoadTranslations(jsonfilename):
    with open(jsonfilename) as jsonfile:
        return json.load(jsonfile)
translateBins = {} if opt.translateBins is None else LoadTranslations(opt.translateBins)
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

# Load fit
from tools.fitter import *
fit = fitter(pois,functions,inputs,False)


mu = np.ones( opt.xnpoints*opt.ynpoints )
# Production scaling
if opt.proc != '':
  xypoi, mu_prod = fit.scaling2D( opt.xpoi, opt.ypoi, opt.proc, npoints=[opt.xnpoints,opt.ynpoints])
  mu = mu*mu_prod
# Dec scaling
if opt.dec != '':
  xypoi, partial = fit.scaling2D( opt.xpoi, opt.ypoi, opt.dec, npoints=[opt.xnpoints,opt.ynpoints])
  xypoi, total = fit.scaling2D( opt.xpoi, opt.ypoi, 'tot', npoints=[opt.xnpoints,opt.ynpoints])
  mu_dec = partial/total
  mu = mu*mu_dec

if opt.doStripCross:
  mu_stripcross = np.ones( opt.xnpoints*opt.ynpoints )
  if opt.proc != '':
    xypoi, mu_prod_stripcross = fit.scaling2D( opt.xpoi, opt.ypoi, "%s_dropcross"%opt.proc, npoints=[opt.xnpoints,opt.ynpoints])
    mu_stripcross = mu_stripcross*mu_prod_stripcross
  if opt.dec != '':
    xypoi, partial_stripcross = fit.scaling2D( opt.xpoi, opt.ypoi, "%s_dropcross"%opt.dec, npoints=[opt.xnpoints,opt.ynpoints])
    xypoi, total_stripcross = fit.scaling2D( opt.xpoi, opt.ypoi, 'tot_dropcross', npoints=[opt.xnpoints,opt.ynpoints])
    mu_dec_stripcross = partial_stripcross/total_stripcross
    mu_stripcross = mu_stripcross*mu_dec_stripcross

if opt.doLinear: 
  fit.setLinearOnly()
  mu_lin = np.ones( opt.xnpoints*opt.ynpoints )
  if opt.proc != '':
    xypoi, mu_prod_lin = fit.scaling2D( opt.xpoi, opt.ypoi, opt.proc, npoints=[opt.xnpoints,opt.ynpoints])
    mu_lin = mu_lin*mu_prod_lin
  if opt.dec != '':
    xypoi, partial_lin = fit.scaling2D( opt.xpoi, opt.ypoi, opt.dec, npoints=[opt.xnpoints,opt.ynpoints])
    xypoi, total_lin = fit.scaling2D( opt.xpoi, opt.ypoi, 'tot', npoints=[opt.xnpoints,opt.ynpoints])
    mu_dec_lin = partial_lin/total_lin
    mu_lin = mu_lin*mu_dec_lin

# Create grid and interpolate
grid_x, grid_y = np.mgrid[ pois[opt.xpoi]['range'][0]:pois[opt.xpoi]['range'][1]:opt.nInterpolatePoints*1j, pois[opt.ypoi]['range'][0]:pois[opt.ypoi]['range'][1]:opt.nInterpolatePoints*1j ]
grid_vals = griddata( xypoi, mu, (grid_x,grid_y), method="cubic")

# Remove NANS
if opt.doProd: 
  grid_vals_prod = griddata( xypoi, mu_prod, (grid_x,grid_y), method="cubic")
  grid_vals_prod = grid_vals_prod[grid_vals==grid_vals]
if opt.doDec: 
  grid_vals_dec = griddata( xypoi, mu_dec, (grid_x,grid_y), method="cubic")
  grid_vals_dec = grid_vals_dec[grid_vals==grid_vals]

if opt.doStripCross: 
  grid_vals_stripcross = griddata( xypoi, mu_stripcross, (grid_x,grid_y), method="cubic")
  grid_vals_stripcross = grid_vals_stripcross[grid_vals==grid_vals]
if opt.doLinear: 
  grid_vals_linear = griddata( xypoi, mu_lin, (grid_x,grid_y), method="cubic")
  grid_vals_linear = grid_vals_linear[grid_vals==grid_vals]

grid_x = grid_x[grid_vals==grid_vals]
grid_y = grid_y[grid_vals==grid_vals]
grid_vals = grid_vals[grid_vals==grid_vals]

# Define Profile2D histogram
h2D = ROOT.TProfile2D("h","h",opt.nBins,pois[opt.xpoi]['range'][0],pois[opt.xpoi]['range'][1],opt.nBins,pois[opt.ypoi]['range'][0],pois[opt.ypoi]['range'][1])
for i in range(len(grid_vals)): h2D.Fill( grid_x[i], grid_y[i], grid_vals[i] )

if opt.doProd: 
  h2D_prod = ROOT.TProfile2D("h_prod","h_prod",opt.nBins,pois[opt.xpoi]['range'][0],pois[opt.xpoi]['range'][1],opt.nBins,pois[opt.ypoi]['range'][0],pois[opt.ypoi]['range'][1])
  for i in range(len(grid_vals_prod)): h2D_prod.Fill( grid_x[i], grid_y[i], grid_vals_prod[i] )

if opt.doDec: 
  h2D_dec = ROOT.TProfile2D("h_dec","h_dec",opt.nBins,pois[opt.xpoi]['range'][0],pois[opt.xpoi]['range'][1],opt.nBins,pois[opt.ypoi]['range'][0],pois[opt.ypoi]['range'][1])
  for i in range(len(grid_vals_dec)): h2D_dec.Fill( grid_x[i], grid_y[i], grid_vals_dec[i] )

if opt.doStripCross: 
  h2D_stripcross = ROOT.TProfile2D("h_stripcross","h_stripcross",opt.nBins,pois[opt.xpoi]['range'][0],pois[opt.xpoi]['range'][1],opt.nBins,pois[opt.ypoi]['range'][0],pois[opt.ypoi]['range'][1])
  for i in range(len(grid_vals_stripcross)): h2D_stripcross.Fill( grid_x[i], grid_y[i], grid_vals_stripcross[i] )
if opt.doLinear: 
  h2D_linear = ROOT.TProfile2D("h_linear","h_linear",opt.nBins,pois[opt.xpoi]['range'][0],pois[opt.xpoi]['range'][1],opt.nBins,pois[opt.ypoi]['range'][0],pois[opt.ypoi]['range'][1])
  for i in range(len(grid_vals_linear)): h2D_linear.Fill( grid_x[i], grid_y[i], grid_vals_linear[i] )

# Loop over bins: if content = 0 then set 999
for ibin in range(1,h2D.GetNbinsX()+1):
  for jbin in range(1,h2D.GetNbinsY()+1):
    if h2D.GetBinContent(ibin,jbin)==0:
      xc, yc = h2D.GetXaxis().GetBinCenter(ibin), h2D.GetYaxis().GetBinCenter(jbin)
      h2D.Fill(xc,yc,999)

# Set up canvas
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

canv = ROOT.TCanvas("canv","canv",600,600)
canv.SetTickx()
canv.SetTicky()
canv.SetLeftMargin(0.115)
canv.SetBottomMargin(0.115)
canv.SetRightMargin(0.18)
# Extract binwidth
xw = abs(pois[opt.xpoi]['range'][0]-pois[opt.xpoi]['range'][1])/opt.nBins
yw = abs(pois[opt.ypoi]['range'][0]-pois[opt.ypoi]['range'][1])/opt.nBins

# POI str
import math
m = "%g"%math.log(1/pois[opt.xpoi]['multiplier'],10)
if m == '1': m = ''
if opt.xpoi == "cWWMinuscB":
  xpstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
  xpstr_stripped = "c_{WW} #minus c_{B}"
else:
  xpstr = "c_{%s} x 10^{%s}"%(opt.xpoi.split("c")[-1],m)
  xpstr_stripped = "c_{%s}"%opt.xpoi.split("c")[-1]

m = "%g"%math.log(1/pois[opt.ypoi]['multiplier'],10)
if m == '1': m = ''
if opt.ypoi == "cWWMinuscB":
  ypstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
  ypstr_stripped = "c_{WW} #minus c_{B}"
else:
  ypstr = "c_{%s} x 10^{%s}"%(opt.ypoi.split("c")[-1],m)
  ypstr_stripped = "c_{%s}"%opt.ypoi.split("c")[-1]

h2D.SetContour(999)
h2D.SetTitle("")
h2D.GetXaxis().SetTitle(xpstr)
h2D.GetXaxis().CenterTitle()
h2D.GetXaxis().SetTitleSize(0.035)
h2D.GetXaxis().SetTitleOffset(1.2)
h2D.GetXaxis().SetRangeUser(pois[opt.xpoi]['range'][0]+xw,pois[opt.xpoi]['range'][1]-xw)
#h2D.GetXaxis().SetRangeUser(-5,5)
h2D.GetYaxis().SetTitle(ypstr)
h2D.GetYaxis().CenterTitle()
h2D.GetYaxis().SetTitleSize(0.035)
h2D.GetYaxis().SetTitleOffset(1.2)
h2D.GetYaxis().SetRangeUser(pois[opt.ypoi]['range'][0]+yw,pois[opt.ypoi]['range'][1]-yw)
#h2D.GetYaxis().SetRangeUser(-10,7.5)

if( opt.proc != '' )&( opt.dec != '' ): h2D.GetZaxis().SetTitle("#mu^{i,f}(%s,%s)"%(xpstr_stripped,ypstr_stripped))
elif( opt.proc != '' ): h2D.GetZaxis().SetTitle("#mu^{i}_{prod}(%s,%s)"%(xpstr_stripped,ypstr_stripped))
elif( opt.dec != '' ): h2D.GetZaxis().SetTitle("#mu^{f}_{dec}(%s,%s)"%(xpstr_stripped,ypstr_stripped))
else: h2D.GetZaxis().SetTitle("NULL")
h2D.GetZaxis().SetTitleSize(0.035)
h2D.GetZaxis().SetTitleOffset(1.2)

# Set maximum
h2D.SetMaximum(10.)
h2D.Draw("COLZ")

# Add lines
lines = {}
lines['v'] = ROOT.TLine(0,pois[opt.ypoi]['range'][0]+yw,0,pois[opt.ypoi]['range'][1]-yw)
lines['v'].SetLineColorAlpha(12,0.5)
lines['v'].SetLineStyle(2)
lines['v'].SetLineWidth(1)
lines['v'].Draw("SAME")

lines['h'] = ROOT.TLine(pois[opt.xpoi]['range'][0]+xw,0,pois[opt.xpoi]['range'][1]-xw,0)
lines['h'].SetLineColorAlpha(12,0.5)
lines['h'].SetLineStyle(2)
lines['h'].SetLineWidth(1)
lines['h'].Draw("SAME")


# Make CI contours
#contours = [0.25,0.5,0.6667,1.0,1.5,2.0,4.0]
#contours = [0.25,0.5,1.0,2.0,4.0]
contours = [0.667,1.0,2.0,4.0]
contours_dec = [0.667,1.0,1.5]
hcontours = od()
if opt.doProd:
  for cidx, c in enumerate(contours):
    hcontours["cprod_%g"%cidx] = h2D_prod.Clone()
    hcontours["cprod_%g"%cidx].SetContour(2)
    hcontours["cprod_%g"%cidx].SetContourLevel(1,c)
    hcontours["cprod_%g"%cidx].SetLineWidth(2)
    if c != 1.: hcontours["cprod_%g"%cidx].SetLineStyle(2)
    hcontours["cprod_%g"%cidx].SetLineColor(ROOT.kRed)
    hcontours["cprod_%g"%cidx].Draw("cont3same")

if opt.doDec:
  for cidx, c in enumerate(contours_dec):
    hcontours["cdec_%g"%cidx] = h2D_dec.Clone()
    hcontours["cdec_%g"%cidx].SetContour(2)
    hcontours["cdec_%g"%cidx].SetContourLevel(1,c)
    hcontours["cdec_%g"%cidx].SetLineWidth(2)
    if c != 1.: hcontours["cdec_%g"%cidx].SetLineStyle(2)
    hcontours["cdec_%g"%cidx].SetLineColor(ROOT.kAzure+1)
    hcontours["cdec_%g"%cidx].Draw("cont3same")

if opt.doLinear:
  for cidx, c in enumerate(contours):
    hcontours["clin_%g"%cidx] = h2D_linear.Clone()
    hcontours["clin_%g"%cidx].SetContour(2)
    hcontours["clin_%g"%cidx].SetContourLevel(1,c)
    hcontours["clin_%g"%cidx].SetLineWidth(1)
    if c != 1.: hcontours["clin_%g"%cidx].SetLineStyle(2)
    hcontours["clin_%g"%cidx].SetLineColor(ROOT.kAzure+10)
    hcontours["clin_%g"%cidx].Draw("cont3same")

if opt.doStripCross:
  for cidx, c in enumerate(contours):
    hcontours["cstr_%g"%cidx] = h2D_stripcross.Clone()
    hcontours["cstr_%g"%cidx].SetContour(2)
    hcontours["cstr_%g"%cidx].SetContourLevel(1,c)
    hcontours["cstr_%g"%cidx].SetLineWidth(1)
    if c != 1.: hcontours["cstr_%g"%cidx].SetLineStyle(2)
    hcontours["cstr_%g"%cidx].SetLineColor(ROOT.kBlack)
    hcontours["cstr_%g"%cidx].Draw("cont3same")

for cidx, c in enumerate(contours):
  hcontours["c_%g"%cidx] = h2D.Clone()
  hcontours["c_%g"%cidx].SetContour(2)
  hcontours["c_%g"%cidx].SetContourLevel(1,c)
  hcontours["c_%g"%cidx].SetLineWidth(3)
  if c != 1.: hcontours["c_%g"%cidx].SetLineStyle(2)
  #hcontours["c_%g"%cidx].SetLineColor(ROOT.kRed)
  hcontours["c_%g"%cidx].SetLineColor(ROOT.kBlack)
  hcontours["c_%g"%cidx].Draw("cont3same")

leg = ROOT.TLegend(0.54,0.7,0.79,0.85)
leg.SetFillStyle(1001)
leg.SetLineColor(1)
leg.SetLineWidth(1)
leg.SetTextSize(0.0275)
leg.AddEntry(hcontours["c_%g"%contours.index(1.0)],"Total scaling","L")
#if opt.doStripCross: leg.AddEntry(hcontours["cstr_%g"%contours.index(1.0)],"Drop B_{pr} (p #neq r)","L")
if opt.doProd: leg.AddEntry(hcontours["cprod_%g"%contours.index(1.0)],"Production only","L")
if opt.doDec: leg.AddEntry(hcontours["cdec_%g"%contours_dec.index(1.0)],"Decay only","L")
leg.Draw("Same")
    
lat0 = ROOT.TLatex()
lat0.SetTextFont(42)
lat0.SetTextAlign(11)
lat0.SetNDC()
lat0.SetTextSize(0.035)
lat0.DrawLatex(0.115,0.92,"HEL UFO")

lat1 = ROOT.TLatex()
lat1.SetTextFont(42)
lat1.SetTextAlign(31)
lat1.SetNDC()
lat1.SetTextSize(0.03)
if( opt.proc != '' )&( opt.dec != '' ): lat1.DrawLatex(0.82,0.92,"( %s , %s )"%(Translate(opt.proc,translateBins),Translate(opt.dec,translateChannels)))
elif( opt.proc != '' ): lat1.DrawLatex(0.82,0.92,"%s"%Translate(opt.proc,translateBins))
elif( opt.dec != '' ): lat1.DrawLatex(0.82,0.92,"%s"%Translate(opt.dec,translateChannels))

lat2 = ROOT.TLatex()
lat2.SetTextFont(42)
lat2.SetTextAlign(23)
lat2.SetTextSize(0.025)
lat2.SetTextAlign(22)
#lat2.SetTextAngle(-30)
#lat2.DrawLatex(-13,4.4,"#color[2]{0.25}")
#lat2.DrawLatex(-12,6.2,"#color[2]{#mu = 0.5}")
#lat2.DrawLatex(-11,7.7,"#color[2]{#mu = 1.0}")
#lat2.DrawLatex(-10,10.,"#color[2]{#mu = 2.0}")
#lat2.DrawLatex(-8,12.4,"#color[2]{#mu = 4.0}")
lat2.SetTextAngle(-27)
lat2.DrawLatex(-8,7.1,"#color[1]{#mu = 0.67}")
lat2.DrawLatex(-7,8.1,"#color[1]{#mu = 1.0}")
lat2.DrawLatex(-6,9.7,"#color[1]{#mu = 2.0}")
lat2.DrawLatex(-4,11.3,"#color[1]{#mu = 4.0}")
#lat2.DrawLatex(-12,6.2,"#color[2]{#mu = 0.5}")
#lat2.DrawLatex(-11,7.7,"#color[2]{#mu = 1.0}")
#lat2.DrawLatex(-10,10.,"#color[2]{#mu = 2.0}")
#lat2.DrawLatex(-8,12.4,"#color[2]{#mu = 4.0}")
lat2.SetTextAngle(-48)
lat2.DrawLatex(-2,-9.2,"#color[861]{#mu_{decay} = 0.67}")
lat2.DrawLatex(9.2,-9.4,"#color[861]{#mu_{decay} = 1.0}")
lat2.SetTextAngle(-65)
lat2.DrawLatex(11.5,5.5,"#color[861]{#mu_{decay} = 1.5}")





  
canv.Update()
canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/qqH_BSM_hzz_cWWMinuscB_vs_cHW.png")
canv.SaveAs("/eos/home-j/jlangfor/www/CMS/thesis/chapter7/scaling_functions/qqH_BSM_hzz_cWWMinuscB_vs_cHW.pdf")
