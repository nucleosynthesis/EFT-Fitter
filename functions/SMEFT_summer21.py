from collections import OrderedDict as od
import json
import numpy as np

params = ['mu_XS_simVBF_lowmjj_highpthjj_BR_ZZ', 'mu_XS_VBF_highMJJ_PTHGT200_BR_gamgam', 'mu_XS_ggH_PTH_GT200_BR_ZZ', 'mu_XS_ttH_PTH_120_200_BR_gamgam', 'mu_XS_ggH_0J_PTH_0_10_BR_gamgam', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_0_60_BR_tautau', 'mu_XS_ttH_PTH_0_60_BR_ZZ', 'mu_XS_VBF_MJJ_120_350_BR_tautau', 'mu_XS_WH_lep_PTV_GT150_BR_gamgam', 'mu_XS_ttH_PTH_200_300_BR_bb', 'mu_XS_WH_lep_PTV_GT150_BR_WW', 'mu_XS_totZH_lep_PTV_75_150_BR_bb', 'mu_XS_WH_lep_PTV_75_150_BR_gamgam', 'mu_XS_ttH_PTH_0_60_BR_tautau', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_0_60_BR_ZZ', 'mu_XS_ttH_PTH_120_200_BR_bb', 'mu_XS_ggH_1J_PTH_0_60_BR_tautau', 'mu_XS_ttH_PTH_GT300_BR_WW', 'mu_XS_ggH_1J_PTH_60_120_BR_tautau', 'mu_XS_WH_lep_PTV_75_150_BR_WW', 'mu_XS_WH_lep_PTV_150_250_0J_BR_tautau', 'mu_XS_totZH_lep_PTV_0_75_BR_tautau', 'mu_XS_simVBF_lowmjj_BR_WW', 'mu_XS_simVBF_lowmjj_highpthjj_BR_tautau', 'mu_XS_VBF_MJJ_60_120_BR_WW', 'mu_XS_totZH_lep_PTV_0_75_BR_gamgam', 'mu_XS_ttH_PTH_0_60_BR_gamgam', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_60_120_BR_gamgam', 'mu_XS_ggH_PTH_GT450_BR_gamgam', 'mu_XS_ttH_PTH_60_120_BR_bb', 'mu_XS_VBF_rest_BR_ZZ', 'mu_XS_ggH_0J_PTH_0_10_BR_ZZ', 'mu_XS_WH_lep_PTV_75_150_BR_tautau', 'mu_XS_totZH_lep_PTV_GT150_BR_WW', 'mu_XS_ttH_PTH_200_300_BR_ZZ', 'mu_XS_VBF_MJJ_60_120_BR_ZZ', 'mu_XS_ggH_1J_PTH_60_120_BR_WW', 'mu_XS_totZH_lep_PTV_75_150_BR_WW', 'mu_XS_ggH_PTH_200_300_BR_tautau', 'mu_XS_ttH_PTH_120_200_BR_ZZ', 'mu_XS_totZH_lep_PTV_75_150_BR_gamgam', 'mu_XS_VBF_MJJ_60_120_BR_tautau', 'mu_XS_simVBF_lowmjj_lowpthjj_BR_ZZ', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_60_120_BR_tautau', 'mu_XS_WH_lep_PTV_75_150_BR_bb', 'mu_XS_totZH_lep_PTV_GT75_BR_tautau', 'mu_XS_simVBF_highmjj_BR_WW', 'mu_XS_simVBF_lowmjj_lowpthjj_BR_gamgam', 'mu_XS_ttH_PTH_60_120_BR_WW', 'mu_XS_WH_lep_PTV_GT75_BR_ZZ', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_60_120_BR_WW', 'mu_XS_ggH_1J_PTH_0_60_BR_WW', 'mu_XS_simVBF_highmjj_highpthjj_BR_tautau', 'mu_XS_VBF_MJJ_0_60_BR_tautau', 'mu_XS_ggH_1J_PTH_120_200_BR_ZZ', 'mu_XS_totZH_lep_PTV_150_250_GE1J_BR_bb', 'mu_XS_VBF_MJJ_60_120_BR_gamgam', 'mu_XS_ggH_0J_PTH_GT10_BR_ZZ', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_120_200_BR_ZZ', 'mu_XS_totZH_lep_PTV_0_75_BR_WW', 'mu_XS_ggH_0J_PTH_0_10_BR_tautau', 'mu_XS_simVBF_lowmjj_highpthjj_BR_gamgam', 'mu_XS_WH_lep_PTV_0_75_BR_tautau', 'mu_XS_VBF_MJJ_120_350_BR_WW', 'mu_XS_ggH_1J_PTH_120_200_BR_WW', 'mu_XS_VBF_highMJJ_PTHGT200_BR_WW', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_60_120_BR_ZZ', 'mu_XS_ttH_PTH_GT300_BR_ZZ', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_120_200_BR_gamgam', 'mu_XS_simVBF_highmjj_highpthjj_BR_gamgam', 'mu_XS_ggH_PTH_200_300_BR_gamgam', 'mu_XS_simVBF_highmjj_highpthjj_BR_ZZ', 'mu_XS_simVBF_highmjj_lowpthjj_BR_tautau', 'mu_XS_ggH_0J_PTH_0_10_BR_WW', 'mu_XS_ttH_PTH_60_120_BR_gamgam', 'mu_XS_simVBF_lowmjj_lowpthjj_BR_tautau', 'mu_XS_WH_lep_PTV_0_75_BR_ZZ', 'mu_XS_ggH_1J_PTH_0_60_BR_ZZ', 'mu_XS_ggH_1J_PTH_60_120_BR_gamgam', 'mu_XS_ttH_PTH_60_120_BR_tautau', 'mu_XS_ttH_PTH_200_300_BR_gamgam', 'mu_XS_VBF_highMJJ_PTHGT200_BR_tautau', 'mu_XS_ggH_1J_PTH_120_200_BR_tautau', 'mu_XS_ggH_PTH_200_300_BR_WW', 'mu_XS_WH_lep_PTV_150_250_0J_BR_bb', 'mu_XS_ttH_PTH_200_300_BR_tautau', 'mu_XS_ggH_PTH_300_450_BR_gamgam', 'mu_XS_WH_lep_PTV_GT250_BR_bb', 'mu_XS_ttH_PTH_GT300_BR_gamgam', 'mu_XS_totZH_lep_PTV_150_250_0J_BR_bb', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_120_200_BR_tautau', 'mu_XS_ggH_PTH_GT300_BR_WW', 'mu_XS_ggH_1J_PTH_0_60_BR_gamgam', 'mu_XS_WH_lep_PTV_0_75_BR_gamgam', 'mu_XS_ttH_PTH_60_120_BR_ZZ', 'mu_XS_ttH_PTH_120_200_BR_tautau', 'mu_XS_VBF_MJJ_0_60_BR_WW', 'mu_XS_WH_lep_PTV_150_250_GE1J_BR_bb', 'mu_XS_ttH_PTH_0_60_BR_bb', 'mu_XS_ggH_0J_PTH_GT10_BR_WW', 'mu_XS_ttH_PTH_GT300_BR_bb', 'mu_XS_WH_lep_PTV_0_75_BR_WW', 'mu_XS_WH_lep_PTV_GT250_BR_tautau', 'mu_XS_totZH_lep_PTV_GT150_BR_gamgam', 'mu_XS_ggH_1J_PTH_60_120_BR_ZZ', 'mu_XS_VBF_LT2J_BR_WW', 'mu_XS_totZH_lep_PTV_GT250_BR_bb', 'mu_XS_ttH_PTH_GT300_BR_tautau', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_120_200_BR_WW', 'mu_XS_WH_lep_PTV_150_250_GE1J_BR_tautau', 'mu_XS_ttH_PTH_120_200_BR_WW', 'mu_XS_ggH_0J_PTH_GT10_BR_tautau', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_0_60_BR_gamgam', 'mu_XS_ggH_1J_PTH_120_200_BR_gamgam', 'mu_XS_simVBF_highmjj_lowpthjj_BR_ZZ', 'mu_XS_ggH_PTH_300_450_BR_tautau', 'mu_XS_ggH_GE2J_MJJ_0_350_PTH_0_60_BR_WW', 'mu_XS_simVBF_highmjj_lowpthjj_BR_gamgam', 'mu_XS_ttH_PTH_0_60_BR_WW', 'mu_XS_ggH_PTH_GT450_BR_tautau', 'mu_XS_ggH_0J_PTH_GT10_BR_gamgam', 'mu_XS_ttH_PTH_200_300_BR_WW']

# Json file storing Aj, Bjk coefficients
with open("functions/SMEFT/xs_coeffs.json","r") as fj: xs_coeffs = json.load(fj)

with open("functions/SMEFT/dec_coeffs.json","r") as fj: dec_coeffs = json.load(fj)

#pois = {}

# Json file storing merged bin definitions
with open("functions/SMEFT/summer21_prefit_asimov_merge.json","r") as fj: merges = json.load(fj)
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
    def dec_function(pois,rpois=None,rmatrix=None):
      # Nominal poi names
      cnames = pois.keys()
      Nc = len(cnames)
      # Create row vector of A terms
      A = np.zeros(Nc)
      for i,ic in enumerate(cnames): 
        if "A_%s"%ic in dec_coeffs[d]: A[i] = dec_coeffs[d]["A_%s"%ic]
      # Create matrix of B terms (squared on diagonal, cross terms on 1/2 off-diagonal)
      B = np.zeros((Nc,Nc))
      for i,ic in enumerate(cnames):
        for j,jc in enumerate(cnames):
          # Squared terms
          if j==i: 
            if "B_%s_2"%ic in dec_coeffs[d]: B[i][j] = dec_coeffs[d]["B_%s_2"%ic]
          # Cross terms
          elif j>i:
            if "B_%s_%s"%(ic,jc) in dec_coeffs[d]: B[i][j] = dec_coeffs[d]["B_%s_%s"%(ic,jc)]
            elif "B_%s_%s"%(jc,ic) in dec_coeffs[d]: B[i][j] = dec_coeffs[d]["B_%s_%s"%(jc,ic)]
      # Matrix calculation
      # Nominal poi basis
      if( rpois is None ):
        cvec = np.array(list(pois.values()))
        cvecT = cvec.T
        mu = 1+A.dot(cvec)+cvecT.dot(B.dot(cvec))
      # Use rotated basis
      else:
        rcvec = np.array(list(rpois.values()))
        rcvecT = rcvec.T
        R = rmatrix
        RT = R.T
        mu = 1+A.dot(R.dot(rcvec))+rcvecT.dot(RT.dot(B.dot(R.dot(rcvec))))
      return mu
    # Add decay scaling as attribute in case other inputs share same name
    setattr(self,"dec_%s"%d,dec_function)
    return dec_function

  # For production cross section bin scaling
  def createXSFunction(self,p):
    # First check if already created function
    if hasattr(self,"xs_%s"%p): return getattr(self,"xs_%s"%p)
    # Otherwise create function
    def xs_function(pois,rpois=None,rmatrix=None):
      cvec = np.array(list(pois.values()))
      cvecT = cvec.T
      cnames = pois.keys()
      Nc = len(cnames)
      # If merged bin: loop over bins (nominal STXS bins) and add equations with relative fraction
      if p in merges:
        mbins = merges[p]
        mu = 0
        xs_tot = 0
        for mb in mbins:
          xs_tot += XSMap[mb]
          # Create row vector of A terms
          A = np.zeros(Nc)
          for i,ic in enumerate(cnames):
            if "A_%s"%ic in xs_coeffs[mb]: A[i] = xs_coeffs[mb]["A_%s"%ic]
          # Create matrix of B terms (squared on diagonal, cross terms on 1/2 off-diagonal)
          B = np.zeros((Nc,Nc))
          for i,ic in enumerate(cnames):
            for j,jc in enumerate(cnames):
              # Squared terms
              if j==i:
                if "B_%s_2"%ic in xs_coeffs[mb]: B[i][j] = xs_coeffs[mb]["B_%s_2"%ic]
              # Cross terms
              elif j>i:
                if "B_%s_%s"%(ic,jc) in xs_coeffs[mb]: B[i][j] = xs_coeffs[mb]["B_%s_%s"%(ic,jc)]
                elif "B_%s_%s"%(jc,ic) in xs_coeffs[mb]: B[i][j] = xs_coeffs[mb]["B_%s_%s"%(jc,ic)] 
          # Calculate mu term for merged bin, weighted by xs
          # Nominal poi basis
          if( rpois is None ):
            cvec = np.array(list(pois.values()))
            cvecT = cvec.T
            mu += XSMap[mb]*(A.dot(cvec)+cvecT.dot(B.dot(cvec)))
          # Rotated basis
          else:
            rcvec = np.array(list(rpois.values()))
            rcvecT = rcvec.T
            R = rmatrix
            RT = R.T
            mu += XSMap[mb]*(A.dot(R.dot(rcvec))+rcvecT.dot(RT.dot(B.dot(R.dot(rcvec)))))
        # Divide through by total xs and add 1 (SM)
        if xs_tot != 0: mu = 1+(mu/xs_tot)
        else: 
          print(" --> [ERROR] Total cross section for merged bin is zero")
          sys.exit(1)

      # If not a merged STXS bin
      else:
        # Create row vector of A terms
        A = np.zeros(Nc)
        for i,ic in enumerate(cnames):
          if "A_%s"%ic in xs_coeffs[p]: A[i] = xs_coeffs[p]["A_%s"%ic]
        # Create matrix of B terms (squared on diagonal, cross terms on 1/2 off-diagonal)
        B = np.zeros((Nc,Nc))
        for i,ic in enumerate(cnames):
          for j,jc in enumerate(cnames):
            # Squared terms
            if j==i:
              if "B_%s_2"%ic in xs_coeffs[p]: B[i][j] = xs_coeffs[p]["B_%s_2"%ic]
            # Cross terms
            elif j>i:
              if "B_%s_%s"%(ic,jc) in xs_coeffs[p]: B[i][j] = xs_coeffs[p]["B_%s_%s"%(ic,jc)]
              elif "B_%s_%s"%(jc,ic) in xs_coeffs[p]: B[i][j] = xs_coeffs[p]["B_%s_%s"%(jc,ic)] 
        # Matrix calculation
        # Nominal poi basis
        if( rpois is None ):
          cvec = np.array(list(pois.values()))
          cvecT = cvec.T
          mu = 1+A.dot(cvec)+cvecT.dot(B.dot(cvec))
        else:
          rcvec = np.array(list(rpois.values()))
          rcvecT = rcvec.T
          R = rmatrix
          RT = R.T
          mu = 1+A.dot(R.dot(rcvec))+rcvecT.dot(RT.dot(B.dot(R.dot(rcvec))))

      # Return scaling
      return mu

    # Add xs scaling as attribute in case other inputs share same name
    setattr(self,"xs_%s"%p,xs_function)
    return xs_function

  def EFTFunction(self,XS,GAMMA,TOT):
    def eft_function(pois,RPOIS=None,RMATRIX=None):
      return XS(pois,rpois=RPOIS,rmatrix=RMATRIX)*(GAMMA(pois,rpois=RPOIS,rmatrix=RMATRIX)/TOT(pois,rpois=RPOIS,rmatrix=RMATRIX))
    return eft_function

  def addFunction(self,name):
    # Extract production and decay channel
    prod = name.split("_BR_")[0].split("_XS_")[-1]
    dec = name.split("_BR_")[-1]
    # Create production scaling
    xs = self.createXSFunction(prod)
    # Create decay scaling
    gamma = self.createDecayFunction(dec)
    tot = self.createDecayFunction('tot')
    # Create total scaling
    x = self.EFTFunction(xs,gamma,tot)
    setattr(self,name,x)

  def getfunction(self,name):
    return getattr(self,name)

# Create function directory, with scaling functions for each input parameter
myfuncs = functionDirectory("SMEFT_summer21")
for p in params: 
  myfuncs.addFunction(p)

functions = od()
for p in params: 
  functions[p] = myfuncs.getfunction(p)
