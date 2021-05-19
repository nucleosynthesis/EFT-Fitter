from collections import OrderedDict as od
import json

params = [u'0J_PTH_0_10_Tag0', u'0J_PTH_0_10_Tag1', u'0J_PTH_0_10_Tag2', u'0J_PTH_GT10_Tag0', u'0J_PTH_GT10_Tag1', u'0J_PTH_GT10_Tag2', u'1J_PTH_0_60_Tag0', u'1J_PTH_0_60_Tag1', u'1J_PTH_0_60_Tag2', u'1J_PTH_120_200_Tag0', u'1J_PTH_120_200_Tag1', u'1J_PTH_120_200_Tag2', u'1J_PTH_60_120_Tag0', u'1J_PTH_60_120_Tag1', u'1J_PTH_60_120_Tag2', u'GE2J_PTH_0_60_Tag0', u'GE2J_PTH_0_60_Tag1', u'GE2J_PTH_0_60_Tag2', u'GE2J_PTH_120_200_Tag0', u'GE2J_PTH_120_200_Tag1', u'GE2J_PTH_120_200_Tag2', u'GE2J_PTH_60_120_Tag0', u'GE2J_PTH_60_120_Tag1', u'GE2J_PTH_60_120_Tag2', u'PTH_200_300_Tag0', u'PTH_200_300_Tag1', u'PTH_300_450_Tag0', u'PTH_300_450_Tag1', u'PTH_450_650_Tag0', u'PTH_GT650_Tag0', u'THQ_LEP', u'TTH_HAD_PTH_0_60_Tag0', u'TTH_HAD_PTH_0_60_Tag1', u'TTH_HAD_PTH_0_60_Tag2', u'TTH_HAD_PTH_120_200_Tag0', u'TTH_HAD_PTH_120_200_Tag1', u'TTH_HAD_PTH_120_200_Tag2', u'TTH_HAD_PTH_120_200_Tag3', u'TTH_HAD_PTH_200_300_Tag0', u'TTH_HAD_PTH_200_300_Tag1', u'TTH_HAD_PTH_200_300_Tag2', u'TTH_HAD_PTH_60_120_Tag0', u'TTH_HAD_PTH_60_120_Tag1', u'TTH_HAD_PTH_60_120_Tag2', u'TTH_HAD_PTH_GT300_Tag0', u'TTH_HAD_PTH_GT300_Tag1', u'TTH_LEP_PTH_0_60_Tag0', u'TTH_LEP_PTH_0_60_Tag1', u'TTH_LEP_PTH_0_60_Tag2', u'TTH_LEP_PTH_120_200_Tag0', u'TTH_LEP_PTH_120_200_Tag1', u'TTH_LEP_PTH_200_300_Tag0', u'TTH_LEP_PTH_60_120_Tag0', u'TTH_LEP_PTH_60_120_Tag1', u'TTH_LEP_PTH_60_120_Tag2', u'TTH_LEP_PTH_GT300_Tag0', u'VBFLIKEGGH_Tag0', u'VBFLIKEGGH_Tag1', u'VBFTOPO_BSM_Tag0', u'VBFTOPO_BSM_Tag1', u'VBFTOPO_JET3VETO_HIGHMJJ_Tag0', u'VBFTOPO_JET3VETO_HIGHMJJ_Tag1', u'VBFTOPO_JET3VETO_LOWMJJ_Tag0', u'VBFTOPO_JET3VETO_LOWMJJ_Tag1', u'VBFTOPO_JET3_HIGHMJJ_Tag0', u'VBFTOPO_JET3_HIGHMJJ_Tag1', u'VBFTOPO_JET3_LOWMJJ_Tag0', u'VBFTOPO_JET3_LOWMJJ_Tag1', u'VBFTOPO_VHHAD_Tag0', u'VBFTOPO_VHHAD_Tag1', u'VH_MET_Tag0', u'VH_MET_Tag1', u'VH_MET_Tag2', u'WH_LEP_PTV_0_75_Tag0', u'WH_LEP_PTV_0_75_Tag1', u'WH_LEP_PTV_75_150_Tag0', u'WH_LEP_PTV_75_150_Tag1', u'WH_LEP_PTV_GT150_Tag0', u'ZH_LEP_Tag0', u'ZH_LEP_Tag1']

# Json file storing purity matrix with mapping to STXS bins
with open("functions/SMEFT/purity_stage1p2_extended_pruned_fromWS_mH125p38.json","r") as fj: purity = json.load(fj)

# Json file storing Aj, Bjk coefficients
with open("functions/SMEFT/xs_coeffs.json","r") as fj: xs_coeffs = json.load(fj)

with open("functions/SMEFT/dec_coeffs.json","r") as fj: dec_coeffs = json.load(fj)

pois = {}

# Json file storing merged bin definitions
with open("functions/SMEFT/recofit_merge.json","r") as fj: merges = json.load(fj)
# Load STXS cross sections: reweight functions by SM cross section in merging
with open("functions/SMEFT/XS.json","r") as fj: XSMap = json.load(fj)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Directory to store EFT functions
class functionDirectory:
  def __init__(self,name): 
    self.name = name

  # For decay channel scaling
  def createDecayFunction(self,d):
    # First check if already created function
    if hasattr(self,"dec_%s"%d): return getattr(self,"dec_%s"%d)
    # Otherwise create function
    def dec_function(pois):
      eq = 1.
      # Add linear and square terms
      for poi,poiVal in pois.items():
        if "A_%s"%poi in dec_coeffs[d]: eq += dec_coeffs[d]["A_%s"%poi]*poiVal
        if "B_%s_2"%poi in dec_coeffs[d]: eq += dec_coeffs[d]["B_%s_2"%poi]*poiVal*poiVal
      # Add cross terms
      poiNames = pois.keys()
      for ipoi in poiNames:
        for jpoi in poiNames:
          if "B_%s_%s"%(ipoi,jpoi) in dec_coeffs[d]: 
            ipoiVal, jpoiVal = pois[ipoi], pois[jpoi]
            eq += dec_coeffs[d]["B_%s_%s"%(ipoi,jpoi)]*ipoiVal*jpoiVal
      return eq
    # Add decay scaling as attribute in case other inputs share same name
    setattr(self,"dec_%s"%d,dec_function)
    return dec_function

  # For production cross section bin scaling
  def createXSFunction(self,p):
    # First check if already created function
    if hasattr(self,"xs_%s"%p): return getattr(self,"xs_%s"%p)
    # Otherwise create function
    def xs_function(pois):
      eq_contributions = od()
      # Loop over contributions in recobin (from purity matrix) and build function
      # weight according to contribution to recobin
      for stxsbin,spurity in purity[p].items():
        eq_contributions[stxsbin] = 0
        
        # If merged bin: loop over nominal STXS bins and add equations with relative fraction
        if stxsbin in merges:
          mbins = merges[stxsbin]
          xs_tot = 0
          for mb in mbins:
            xs_tot += XSMap[mb]
            # Add linear and square terms
            for poi,poiVal in pois.items():            
              if "A_%s"%poi in xs_coeffs[mb]: eq_contributions[stxsbin] += XSMap[mb]*xs_coeffs[mb]["A_%s"%poi]*poiVal
              if "B_%s_2"%poi in xs_coeffs[mb]: eq_contributions[stxsbin] += XSMap[mb]*xs_coeffs[mb]["B_%s_2"%poi]*poiVal*poiVal
            # Add cross terms
            poiNames = pois.keys()
            for ipoi in poiNames:
              for jpoi in poiNames:
                if "B_%s_%s"%(ipoi,jpoi) in xs_coeffs[mb]: 
                  ipoiVal, jpoiVal = pois[ipoi], pois[jpoi]
                  eq_contributions[stxsbin] += XSMap[mb]*xs_coeffs[mb]["B_%s_%s"%(ipoi,jpoi)]*ipoiVal*jpoiVal
          # Divide through by total xs
          if xs_tot != 0: eq_contributions[stxsbin] = eq_contributions[stxsbin]/xs_tot
          else: 
            print(" --> [ERROR] Total cross section for merged bin is zero")
            sys.exit(1)

        # If not a merged STXS bin
        else:
          # Add linear and square terms
          for poi,poiVal in pois.items():
            if "A_%s"%poi in xs_coeffs[stxsbin]: eq_contributions[stxsbin] += xs_coeffs[stxsbin]["A_%s"%poi]*poiVal
            if "B_%s_2"%poi in xs_coeffs[stxsbin]: eq_contributions[stxsbin] += xs_coeffs[stxsbin]["B_%s_2"%poi]*poiVal*poiVal
          # Add cross terms
          poiNames = pois.keys()
          for ipoi in poiNames:
            for jpoi in poiNames:
              if "B_%s_%s"%(ipoi,jpoi) in xs_coeffs[stxsbin]: 
                ipoiVal, jpoiVal = pois[ipoi], pois[jpoi]
                eq_contributions[stxsbin] += xs_coeffs[stxsbin]["B_%s_%s"%(ipoi,jpoi)]*ipoiVal*jpoiVal

        # Weight contribution by purity in recobin
        eq_contributions[stxsbin] = eq_contributions[stxsbin]*spurity 

      # Sum contributions
      eq = 1.
      for v in eq_contributions.values(): eq += v
      # Return final equation
      return eq
    # Add xs scaling as attribute in case other inputs share same name
    setattr(self,"xs_%s"%p,xs_function)
    return xs_function

  def EFTFunction(self,XS,GAMMA,TOT):
    def eft_function(pois):
      return XS(pois)*(GAMMA(pois)/TOT(pois))
    return eft_function

  def addFunction(self,name):
    # Extract production and decay channel
    recobin = name
    dec = "gamgam"
    # Create production scaling
    xs = self.createXSFunction(recobin)
    # Create decay scaling
    gamma = self.createDecayFunction(dec)
    tot = self.createDecayFunction('tot')
    # Create total scaling
    x = self.EFTFunction(xs,gamma,tot)
    setattr(self,name,x)

  def getfunction(self,name):
    return getattr(self,name)

# Create function directory, with scaling functions for each input parameter
myfuncs = functionDirectory("SMEFT_recofit")
for p in params: 
  myfuncs.addFunction(p)

functions = od()
for p in params: 
  functions[p] = myfuncs.getfunction(p)
