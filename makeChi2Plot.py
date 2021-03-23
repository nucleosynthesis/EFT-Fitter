import ROOT
import pickle
import numpy as np
from collections import OrderedDict as od

# For interpolation
from scipy import interpolate

# Input options
from optparse import OptionParser
def get_options():
  parser = OptionParser()
  parser.add_option('--poi', dest='poi', default='cG', help="Main poi to plot")
  parser.add_option('--param_set', dest='param_set', default='params.HEL', help="Definition of parameters (in module)")
  parser.add_option('--otherPOIs', dest='otherPOIs', default='', help="Comma separated list of profiled pois")
  parser.add_option('--inputPkl', dest='inputPkl', default='results.pkl', help="Input pkl file storing results")
  parser.add_option('--doProfiledPOIFrac', dest='doProfiledPOIFrac', default=False, action="store_true", help="Plot fractional profiled poi instead of the pull")
  parser.add_option('--doLinear', dest='doLinear', default=False, action="store_true", help="Add linear lines to plot")
  parser.add_option('--outputDir', dest='outputDir', default='.', help="Output plot directory")
  return parser.parse_args()
(opt,args) = get_options()

if opt.doLinear: modes = ['fixed','fixed_linear','profiled','profiled_linear']
else: modes = ['fixed','profiled']
# Load results
with open(opt.inputPkl,"rb") as fpkl: results = pickle.load(fpkl)

# Extract full list of pois
from importlib import import_module
pois = import_module(opt.param_set).pois

if opt.otherPOIs == "all":
  opoistr = ""
  for poi in pois.keys():
    if poi!=opt.poi: opoistr += "%s,"%poi
  opt.otherPOIs = opoistr[:-1]
  
# Function to extract poi best-fit and +-1/2sigma points
def extractValsV2( _p, _dchi2 ):

  # Best-fit
  i_bf = _dchi2.argmin()
  bf = _p[i_bf]

  # Add union of intervals stuff: if i+1 and i-1 are > i && dchi2 <= 4 (minimum)
  i_min, minimum = [], []
  up1, up2, down1, down2 = [], [], [], []
  dchi2_min = []
  for i in range(1,len(_p)-1):
    if(_dchi2[i] < 4. )&( _dchi2[i-1] > _dchi2[i])&( _dchi2[i+1] > _dchi2[i]): i_min.append(i)

  # Loop over minima
  for i in i_min: 
    minimum.append( _p[i] )
    dchi2_min.append(_dchi2[i])
    
    # Find confidence intervals for each minimum
    found_1upsigma, found_2upsigma = False, False
    if _dchi2[i] < 1.:
      for j in range(i,len(_p)):
        if not found_1upsigma:
          if _dchi2[j] >= 1.:
            j_1up = j
            found_1upsigma = True
    for j in range(i,len(_p)):
      if not found_2upsigma:
        if _dchi2[j] >= 4.:
          j_2up = j
          found_2upsigma = True
    if found_1upsigma: up1.append( _p[j_1up] )
    else: up1.append(None)
    if found_2upsigma: up2.append( _p[j_2up] ) 
    else: up2.append(None)

    found_1downsigma, found_2downsigma = False, False
    if _dchi2[i] < 1.:
      for j in range(i,0,-1):
        if not found_1downsigma:
          if _dchi2[j] >= 1.:
            j_1down = j
            found_1downsigma = True
    for j in range(i,0,-1):
      if not found_2downsigma:
        if _dchi2[j] >= 4.:
          j_2down = j
          found_2downsigma = True
    if found_1downsigma: down1.append( _p[j_1down] )
    else: down1.append(None)
    if found_2downsigma: down2.append( _p[j_2down] ) 
    else: down2.append(None)

  return bf, np.array(minimum), np.array(dchi2_min), np.array(up1), np.array(up2), np.array(down1), np.array(down2)


# Function to extract poi best-fit and +-1/2sigma points
def extractVals( _p, _dchi2 ):

  # Best-fit
  i_bf = _dchi2.argmin()
  bf = _p[i_bf]

  # + confidence intervals
  found_1upsigma, found_2upsigma = False, False
  for i in range(i_bf,len(_p)):
    if not found_1upsigma:
      if _dchi2[i] >= 1.:
        i_1up = i
        found_1upsigma = True
  
    if not found_2upsigma:
      if _dchi2[i] >= 4.:
        i_2up = i
        found_2upsigma = True
  
  if found_1upsigma: up1 = abs(_p[i_1up]-bf)
  else: up1 = None
  if found_2upsigma: up2 = abs(_p[i_2up]-bf)
  else: up2 = None
  
  found_1downsigma, found_2downsigma = False, False
  for i in range(i_bf,0,-1):
    if not found_1downsigma:
      if _dchi2[i] >= 1.:
        i_1down = i
        found_1downsigma = True
  
    if not found_2downsigma:
      if _dchi2[i] >= 4.:
        i_2down = i
        found_2downsigma = True

  if found_1downsigma: down1 = -1*abs(_p[i_1down]-bf)
  else: down1 = None
  if found_2downsigma: down2 = -1*abs(_p[i_2down]-bf)
  else: down2 = None

  return bf, up1, up2, down1, down2
      
  

# Do interpolations and make dchi2 graphs
grs = od()
for poi in pois:
  for mode in modes: #['fixed','profiled','fixed_linear','profiled_linear']:
    print(poi)
    #dchi2 = results[poi][mode]['dchi2'][1:]
    dchi2 = results[poi][mode]['dchi2']
    #p = results[poi][mode]['pvals'][1:]
    p = results[poi][mode]['pvals']
    if"linear" in mode: f = interpolate.interp1d(p,dchi2)
    else: f = interpolate.interp1d(p,dchi2,"cubic")
    pext = np.linspace( p.min(), p.max(), 10000 )
    dchi2ext = f(pext)
    results[poi][mode]['pvals_ext'] = pext
    results[poi][mode]['dchi2_ext'] = dchi2ext

    bf, up1, up2, down1, down2 = extractVals(pext,dchi2ext)
    results[poi][mode]['bestfit'] = bf
    results[poi][mode]['up01sigma'] = up1 
    results[poi][mode]['up02sigma'] = down1 
    results[poi][mode]['down01sigma'] = down1
    results[poi][mode]['down02sigma'] = down2

    bf, minimum, dchi2_min, up1, up2, down1, down2 = extractValsV2(pext,dchi2ext)
    results[poi][mode]['bestfitv2'] = bf
    results[poi][mode]['minimumv2'] = minimum
    results[poi][mode]['dchi2_minv2'] = dchi2_min
    results[poi][mode]['up01crossv2'] = up1
    results[poi][mode]['up02crossv2'] = up2
    results[poi][mode]['down01crossv2'] = down1
    results[poi][mode]['down02crossv2'] = down2

    if poi == opt.poi: 
      grs["%s_%s"%(poi,mode)] = ROOT.TGraph()
      grs["%s_%s_ext"%(poi,mode)] = ROOT.TGraph()
      for i in range(len(p)): grs["%s_%s"%(poi,mode)].SetPoint( grs["%s_%s"%(poi,mode)].GetN(), p[i], dchi2[i] )
      for i in range(len(pext)): grs["%s_%s_ext"%(poi,mode)].SetPoint( grs["%s_%s_ext"%(poi,mode)].GetN(), pext[i], dchi2ext[i] )

# Make profiled poi curves
mode = "profiled"
if opt.otherPOIs != '':
  for opoi in opt.otherPOIs.split(","):
    #p = results[opt.poi][mode]['pvals'][1:]
    p = results[opt.poi][mode]['pvals']
    pext = np.linspace( p.min(), p.max(), 10000 )
    #allpvals = results[opt.poi][mode]['allpvals'][1:]
    allpvals = results[opt.poi][mode]['allpvals']

    # Extract index of other POI and fill values
    i_op = pois.keys().index(opoi)
    op = []
    for j in range(len(allpvals)): op.append( allpvals[j][i_op] )
    op = np.array(op)

    # Variation with respect to minimum
    fop = (op-pois[opoi]['range'][0])/abs(pois[opoi]['range'][1]-pois[opoi]['range'][0])

    # Pull with respect to closest minimum
    obf = results[opoi][mode]['bestfit']
    oup1 = results[opoi][mode]['up01sigma']
    odown1 = results[opoi][mode]['down01sigma']
    pullop = []
    for i in range(len(op)):
      #if op[i]>obf: pullop.append((op[i]-obf)/abs(odown1))
      #else: pullop.append((op[i]-obf)/abs(oup1))

      iomin = abs(results[opoi][mode]['minimumv2']-op[i]).argmin()
      omin = results[opoi][mode]['minimumv2'][iomin]
      if results[opoi][mode]['up01crossv2'][iomin] == None: oup1 = 0.5*abs(results[opoi][mode]['up02crossv2'][iomin]-omin)
      else: oup1 = abs(results[opoi][mode]['up01crossv2'][iomin]-omin)
      if results[opoi][mode]['down01crossv2'][iomin] == None: odown1 = -0.5*abs(results[opoi][mode]['down02crossv2'][iomin]-omin)
      else: odown1 = -1*abs(results[opoi][mode]['down01crossv2'][iomin]-omin)

      if op[i]>omin: pullop.append((op[i]-omin)/abs(odown1))
      else: pullop.append((op[i]-omin)/abs(oup1))
    pullop = np.array(pullop)
    
    # Interpolate
    f_frac = interpolate.interp1d(p,fop,'cubic')
    fop_ext = f_frac(pext)
    f_pull = interpolate.interp1d(p,pullop,'cubic')
    pullop_ext = f_pull(pext)

    # Make graphs
    grs["%s_profiledFrac_%s"%(opt.poi,opoi)] = ROOT.TGraph()
    grs["%s_profiledFrac_ext_%s"%(opt.poi,opoi)] = ROOT.TGraph()
    grs["%s_profiledPull_%s"%(opt.poi,opoi)] = ROOT.TGraph()
    grs["%s_profiledPull_ext_%s"%(opt.poi,opoi)] = ROOT.TGraph()
    for i in range(len(p)):
      grs["%s_profiledFrac_%s"%(opt.poi,opoi)].SetPoint(grs["%s_profiledFrac_%s"%(opt.poi,opoi)].GetN(),p[i],fop[i])
      grs["%s_profiledPull_%s"%(opt.poi,opoi)].SetPoint(grs["%s_profiledPull_%s"%(opt.poi,opoi)].GetN(),p[i],pullop[i])
    for i in range(len(pext)):
      grs["%s_profiledFrac_ext_%s"%(opt.poi,opoi)].SetPoint(grs["%s_profiledFrac_ext_%s"%(opt.poi,opoi)].GetN(),pext[i],fop_ext[i])
      grs["%s_profiledPull_ext_%s"%(opt.poi,opoi)].SetPoint(grs["%s_profiledPull_ext_%s"%(opt.poi,opoi)].GetN(),pext[i],pullop_ext[i])

  
    
 
# Plot
styleMap = od()
styleMap['fixed_ext'] = {'LineColor':9,'LineWidth':2,'LineStyle':1,'MarkerSize':0}
styleMap['fixed'] = {'MarkerColor':9,'LineWidth':0,'MarkerSize':.75,'MarkerStyle':20}
styleMap['fixed_dummy'] = {'LineColor':9,'MarkerColor':9,'LineWidth':2,'LineStyle':1,'MarkerSize':.75,'MarkerStyle':20}

styleMap['fixed_linear_ext'] = {'LineColor':9,'LineWidth':1,'LineStyle':2,'MarkerSize':0}
styleMap['fixed_linear'] = {'MarkerColor':9,'LineWidth':0,'MarkerSize':.75,'MarkerStyle':24}
styleMap['fixed_linear_dummy'] = {'LineColor':9,'MarkerColor':9,'LineWidth':1,'LineStyle':2,'MarkerSize':.75,'MarkerStyle':24}

styleMap['profiled_ext'] = {'LineColor':ROOT.kBlack,'LineWidth':2,'MarkerSize':0}
styleMap['profiled'] = {'MarkerColor':ROOT.kBlack,'LineWidth':0,'MarkerSize':.75,'MarkerStyle':20}
styleMap['profiled_dummy'] = {'LineColor':ROOT.kBlack,'MarkerColor':ROOT.kBlack,'LineWidth':2,'MarkerSize':.75,'MarkerStyle':20}

styleMap['profiled_linear_ext'] = {'LineColor':ROOT.kBlack,'LineWidth':1,'LineStyle':2,'MarkerSize':0}
styleMap['profiled_linear'] = {'MarkerColor':ROOT.kBlack,'LineWidth':0,'MarkerSize':.75,'MarkerStyle':24}
styleMap['profiled_linear_dummy'] = {'LineColor':ROOT.kBlack,'MarkerColor':ROOT.kBlack,'LineWidth':1,'LineStyle':2,'MarkerSize':.75,'MarkerStyle':24}

styleMap['profiledPull_ext'] = {'LineWidth':2,'LineStyle':1,'MarkerSize':0}
styleMap['profiledPull'] = {'LineWidth':0,'MarkerSize':0,'MarkerStyle':24}
styleMap['profiledFrac_ext'] = {'LineWidth':2,'LineStyle':1,'MarkerSize':0}
styleMap['profiledFrac'] = {'LineWidth':0,'MarkerSize':0,'MarkerStyle':24}

colorMap = od()
colorMap['cG'] = {'LineColor':ROOT.kAzure+1,'MarkerColor':ROOT.kAzure+1}
colorMap['cA'] = {'LineColor':ROOT.kRed-4,'MarkerColor':ROOT.kRed-4}
colorMap['cHW'] = {'LineColor':ROOT.kGreen-3,'MarkerColor':ROOT.kGreen-3}
colorMap['cWWMinuscB'] = {'LineColor':ROOT.kOrange-3,'MarkerColor':ROOT.kOrange-3}
colorMap['cu'] = {'LineColor':ROOT.kMagenta-7,'MarkerColor':ROOT.kMagenta-7}
colorMap['cd'] = {'LineColor':ROOT.kViolet+6,'MarkerColor':ROOT.kViolet+6}
#colorMap['cd'] = {'LineColor':ROOT.kYellow+1,'MarkerColor':ROOT.kYellow+1}
colorMap['cl'] = {'LineColor':ROOT.kCyan-7,'MarkerColor':ROOT.kCyan-7}
#colorMap['cl'] = {'LineColor':ROOT.kViolet+6,'MarkerColor':ROOT.kViolet+6}

# POI str
import math
m = "%g"%math.log(1/pois[opt.poi]['multiplier'],10)
if m == '1': m = ''
if opt.poi == "cWWMinuscB":
  pstr = "(c_{WW} #minus c_{B}) x 10^{%s}"%m
else:
  pstr = "c_{%s} x 10^{%s}"%(opt.poi.split("c")[-1],m)


ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

if opt.otherPOIs == '':
  canv = ROOT.TCanvas()
  canv.SetBottomMargin(0.15)
  canv.SetTickx()
  canv.SetTicky()

else:
  canv = ROOT.TCanvas("canv","canv",600,600)
  pad1 = ROOT.TPad("pad1","pad1",0,0.25,1,1)
  pad1.SetTickx()
  pad1.SetTicky()
  pad1.SetBottomMargin(0.25)
  pad1.SetLeftMargin(0.12)
  pad1.Draw()
  pad2 = ROOT.TPad("pad2","pad2",0,0,1,0.35)
  pad2.SetTickx()
  pad2.SetTicky()
  pad2.SetTopMargin(0.15)
  pad2.SetBottomMargin(0.25)
  pad2.SetLeftMargin(0.12)
  pad2.Draw()
  padSizeRatio = 0.75/0.35
  pad1.cd()

h_axes = ROOT.TH1F("haxes","",100, results[opt.poi][mode]['pvals_ext'].min(), results[opt.poi][mode]['pvals_ext'].max() )
h_axes.SetMaximum(15.)
h_axes.SetMinimum(0.)
h_axes.SetTitle("")
if opt.otherPOIs == '':
  h_axes.GetXaxis().SetTitle(pstr)
  #h_axes.GetXaxis().SetTitleOffset(0.9)
  h_axes.GetXaxis().SetTitleSize(0.05)
  h_axes.GetXaxis().SetLabelSize(0.035)
else:
  h_axes.GetXaxis().SetLabelSize(0.)
  h_axes.GetXaxis().SetTitle(pstr)
  #h_axes.GetXaxis().SetTitleOffset(0.9)
  h_axes.GetXaxis().SetTitleSize(0.05)
  h_axes.GetXaxis().SetLabelSize(0.035)

h_axes.GetYaxis().SetTitle("#Delta#chi^{2}")
h_axes.GetYaxis().SetTitleSize(0.05)
h_axes.GetYaxis().SetTitleOffset(0.8)
h_axes.GetYaxis().SetLabelSize(0.035)
h_axes.GetYaxis().SetLabelOffset(0.007)
h_axes.GetYaxis().CenterTitle()
h_axes.Draw()


for mode in modes:
  for k,v in styleMap[mode].items(): getattr(grs["%s_%s"%(opt.poi,mode)],"Set%s"%k)(v)
  for k,v in styleMap["%s_ext"%mode].items(): getattr(grs["%s_%s_ext"%(opt.poi,mode)],"Set%s"%k)(v)
  grs["%s_%s"%(opt.poi,mode)].Draw("Same P")
  grs["%s_%s_ext"%(opt.poi,mode)].Draw("Same C")

# Lines
hlines = {}
yvals = [1,4]
for i in range(len(yvals)):
  yval = yvals[i]
  hlines['hline_%g'%i] = ROOT.TLine(results[opt.poi][mode]['pvals_ext'].min(),yval,results[opt.poi][mode]['pvals_ext'].max(),yval)
  hlines['hline_%g'%i].SetLineColorAlpha(ROOT.kRed,0.5)
  hlines['hline_%g'%i].SetLineStyle(2)
  hlines['hline_%g'%i].SetLineWidth(1)
  hlines['hline_%g'%i].Draw("SAME")

# Box with legend
x_range = results[opt.poi][mode]['pvals_ext'].max()-results[opt.poi][mode]['pvals_ext'].min()
box = ROOT.TBox(results[opt.poi][mode]['pvals_ext'].min()+0.01*x_range, h_axes.GetMaximum()*0.8, results[opt.poi][mode]['pvals_ext'].max()-0.01*x_range, h_axes.GetMaximum()*0.995)
box.SetFillStyle(1001)
box.SetFillColor(ROOT.kWhite)
box.Draw("Same")
yval = h_axes.GetMaximum()*0.8
hlines['hline_box'] = ROOT.TLine(results[opt.poi][mode]['pvals_ext'].min(),yval,results[opt.poi][mode]['pvals_ext'].max(),yval)
hlines['hline_box'].SetLineColorAlpha(ROOT.kBlack,0.5)
hlines['hline_box'].SetLineWidth(1)
hlines['hline_box'].Draw("SAME")

# Dummy graphs
grs_dummy = od()
for mode in modes:
  grs_dummy["%s_%s"%(opt.poi,mode)] = ROOT.TGraph()
  for k,v in styleMap["%s_dummy"%mode].items(): getattr(grs_dummy["%s_%s"%(opt.poi,mode)],"Set%s"%k)(v)

if opt.doLinear: 
  leg = ROOT.TLegend(0.15,0.78,0.85,0.89)
  leg.SetNColumns(2)
  leg.SetFillStyle(0)
  leg.SetLineColor(0)
  leg.SetTextSize(0.032)
  leg.AddEntry( grs_dummy["%s_profiled"%opt.poi], "Profiled", "LP") 
  leg.AddEntry( grs_dummy["%s_profiled_linear"%opt.poi], "Profiled (Lin. terms only)", "LP") 
  leg.AddEntry( grs_dummy["%s_fixed"%opt.poi], "Other c_{p} = 0", "LP") 
  leg.AddEntry( grs_dummy["%s_fixed_linear"%opt.poi], "Other c_{p} = 0 (Lin. terms only)", "LP") 
  leg.Draw("Same")
else:
  leg = ROOT.TLegend(0.15,0.78,0.55,0.89)
  leg.SetFillStyle(0)
  leg.SetLineColor(0)
  leg.SetTextSize(0.032)
  leg.AddEntry( grs_dummy["%s_profiled"%opt.poi], "Profiled", "LP") 
  leg.AddEntry( grs_dummy["%s_fixed"%opt.poi], "Other c_{p} = 0", "LP") 
  leg.Draw("Same")


# Text
lat0 = ROOT.TLatex()
lat0.SetTextFont(42)
lat0.SetTextAlign(31)
lat0.SetNDC()
lat0.SetTextSize(0.045)
lat0.DrawLatex(0.9,0.92,"35.9-137 fb^{-1} (13 TeV)")

lat1 = ROOT.TLatex()
lat1.SetTextFont(42)
lat1.SetTextAlign(12)
lat1.SetTextSize(0.035)
xpos = 0.05*(results[opt.poi]['profiled']['pvals_ext'].max()-results[opt.poi]['profiled']['pvals_ext'].min())+results[opt.poi]['profiled']['pvals_ext'].min()
xposinv = results[opt.poi]['profiled']['pvals_ext'].max()-0.05*(results[opt.poi]['profiled']['pvals_ext'].max()-results[opt.poi]['profiled']['pvals_ext'].min())
lat1.DrawLatex(xpos,1.,"#color[2]{#bf{1#sigma}}")
lat1.DrawLatex(xpos,4.,"#color[2]{#bf{2#sigma}}")
  
if opt.otherPOIs != '':
  pad2.cd()
  h_axes_ratio = ROOT.TH1F("haxes_ratio","",100,results[opt.poi][mode]['pvals_ext'].min(), results[opt.poi][mode]['pvals_ext'].max() )
  if opt.doProfiledPOIFrac:
    h_axes_ratio.SetMaximum(1.15)
    h_axes_ratio.SetMinimum(-0.15)
  else:
    h_axes_ratio.SetMaximum(2.9)
    h_axes_ratio.SetMinimum(-2.9)
  h_axes_ratio.SetTitle("")
  h_axes_ratio.GetXaxis().SetTitle(pstr)
  h_axes_ratio.GetXaxis().SetTitleSize(0.05*padSizeRatio)
  h_axes_ratio.GetXaxis().SetLabelSize(0.035*padSizeRatio)
  h_axes_ratio.GetXaxis().SetLabelOffset(0.007)
  h_axes_ratio.GetXaxis().SetTickLength(0.03*padSizeRatio)
  h_axes_ratio.GetYaxis().SetLabelSize(0.035*padSizeRatio)
  h_axes_ratio.GetYaxis().SetTitleSize(0.05*padSizeRatio)
  h_axes_ratio.GetYaxis().SetTitleOffset(0.8/padSizeRatio)
  h_axes_ratio.GetYaxis().SetLabelOffset(0.007)
  h_axes_ratio.GetYaxis().CenterTitle()
  if opt.doProfiledPOIFrac: h_axes_ratio.GetYaxis().SetTitle("#Delta(c)")
  else: h_axes_ratio.GetYaxis().SetTitle("(c-#hat{c})/#sigma_{c}")
  h_axes_ratio.SetLineWidth(0)
  h_axes_ratio.Draw()

  for opoi in opt.otherPOIs.split(","):
    if opt.doProfiledPOIFrac: mode = "profiledFrac"
    else: mode = "profiledPull"
    for k,v in styleMap[mode].items(): getattr(grs["%s_%s_%s"%(opt.poi,mode,opoi)],"Set%s"%k)(v)
    for k,v in styleMap["%s_ext"%mode].items(): getattr(grs["%s_%s_ext_%s"%(opt.poi,mode,opoi)],"Set%s"%k)(v)
    for k,v in colorMap[opoi].items(): 
      getattr(grs["%s_%s_%s"%(opt.poi,mode,opoi)],"Set%s"%k)(v)
      getattr(grs["%s_%s_ext_%s"%(opt.poi,mode,opoi)],"Set%s"%k)(v)

    grs["%s_%s_%s"%(opt.poi,mode,opoi)].Draw("Same P")
    grs["%s_%s_ext_%s"%(opt.poi,mode,opoi)].Draw("Same C")

  # Draw lines
  hlines_r = {}
  if opt.doProfiledPOIFrac: yvals = [0,0.5,1]
  else: yvals = [-2,-1,0,1,2]
  for i in range(len(yvals)):
    yval = yvals[i]
    hlines_r['hline_%g'%i] = ROOT.TLine(results[opt.poi]['profiled']['pvals_ext'].min(),yval,results[opt.poi]['profiled']['pvals_ext'].max(),yval)
    hlines_r['hline_%g'%i].SetLineColorAlpha(ROOT.kGray,0.5)
    hlines_r['hline_%g'%i].SetLineStyle(2)
    hlines_r['hline_%g'%i].SetLineWidth(1)
    hlines_r['hline_%g'%i].Draw("SAME")

  # Text
  lat2 = ROOT.TLatex()
  lat2.SetTextFont(42)
  lat2.SetTextAlign(12)
  lat2.SetTextSize(0.035*padSizeRatio)
  xpos = 0.05*(results[opt.poi]['profiled']['pvals_ext'].max()-results[opt.poi]['profiled']['pvals_ext'].min())+results[opt.poi]['profiled']['pvals_ext'].min()
  xposinv = results[opt.poi]['profiled']['pvals_ext'].max()-0.05*(results[opt.poi]['profiled']['pvals_ext'].max()-results[opt.poi]['profiled']['pvals_ext'].min())
  if opt.doProfiledPOIFrac: 
    lat2.DrawLatex(xpos,0.,"#color[15]{#bf{c_{min}}}")
    lat2.DrawLatex(xpos,1.0,"#color[15]{#bf{c_{max}}}")
  else:
    lat2.DrawLatex(xposinv,-2,"#color[15]{#bf{-2#sigma}}")
    lat2.DrawLatex(xposinv,-1,"#color[15]{#bf{-1#sigma}}")
    lat2.DrawLatex(xposinv,1,"#color[15]{#bf{1#sigma}}")
    lat2.DrawLatex(xposinv,2,"#color[15]{#bf{2#sigma}}")

  # Legend
  legs = {}
  nopois = len(opt.otherPOIs.split(","))
  legLength = (1.0-0.12)/(nopois)
  for iop, opoi in enumerate(opt.otherPOIs.split(",")):
    xpos = 0.12+iop*legLength
    legs[opoi] = ROOT.TLegend(xpos,0.86,xpos+legLength,1.0)
    legs[opoi].SetFillStyle(0)
    legs[opoi].SetLineColor(0)
    legs[opoi].SetTextSize(0.04*padSizeRatio)
    if opoi == "cWWMinuscB": pstr = "#scale[0.8]{(c_{WW} #minus c_{B})}"
    else: pstr = "c_{%s}"%opoi.split("c")[-1]
    legs[opoi].AddEntry( grs["%s_%s_ext_%s"%(opt.poi,mode,opoi)], pstr, "L")
    legs[opoi].Draw("Same")


canv.Update()
canv.SaveAs("%s/%s.png"%(opt.outputDir,opt.poi))
canv.SaveAs("%s/%s.pdf"%(opt.outputDir,opt.poi))
