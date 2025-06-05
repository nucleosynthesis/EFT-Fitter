from collections import OrderedDict as od
import sys

def Gamma_hgg(pois):
   k = pois["kappa_gamma"]
   return k*k

def Gamma_gluglu(pois):
   k = pois["kappa_glu"]
   return k*k

def Gamma_tot(pois): 
   kc = pois["kappa_c"]
   kb = pois["kappa_b"]
   kZ = pois["kappa_Z"]
   kW = pois["kappa_W"]
   ktau = pois["kappa_tau"]
   kmu = pois["kappa_mu"]
   ks  = pois["kappa_s"]
   kglu = Gamma_gluglu(pois)
   kgam = Gamma_hgg(pois)
   kZg  = pois['kappa_Zg']
   sumtotals = 0.58+0.22+0.08+0.06+0.027+0.029+0.0023+0.0016+0.00025+0.00022
   kH2 = 0.58*kb*kb + 0.22*kW*kW + 0.08*kglu*kglu +\
    0.06*ktau*ktau +0.027*kZ*kZ +0.029*kc*kc +0.0023*kgam*kgam \
    +0.0016*kZg*kZg + 0.00025*ks*ks+0.00022*kmu*kmu
   return kH2/sumtotals

def xs_ZH(pois):
    kZ = pois["kappa_Z"]
    return kZ*kZ

def xsBR_hZZ(pois):
    kZ = pois["kappa_Z"]
    return xs_ZH(pois)*kZ*kZ/Gamma_tot(pois)

def xsBR_hWW(pois):
    kW = pois["kappa_W"]
    return xs_ZH(pois)*kW*kW/Gamma_tot(pois)

def xsBR_hbb(pois):
    kb = pois["kappa_b"]
    return xs_ZH(pois)*kb*kb/Gamma_tot(pois)

def xsBR_hcc(pois):
    kc = pois["kappa_c"]
    return xs_ZH(pois)*kc*kc/Gamma_tot(pois)

def xsBR_hss(pois):
    ks = pois["kappa_s"]
    return xs_ZH(pois)*ks*ks/Gamma_tot(pois)

def xsBR_hgluglu(pois):
    kglu = Gamma_gluglu(pois)
    return xs_ZH(pois)*kglu/Gamma_tot(pois)

def xsBR_hgamgam(pois):
    kgam = Gamma_hgg(pois)
    return xs_ZH(pois)*kgam/Gamma_tot(pois)

def xsBR_hZgam(pois):
    kZg = pois["kappa_Zg"]
    return xs_ZH(pois)*kZg*kZg/Gamma_tot(pois)

def xsBR_hmumu(pois):
    kmu = pois["kappa_mu"]
    return xs_ZH(pois)*kmu*kmu/Gamma_tot(pois)

def xsBR_htautau(pois):
    ktau = pois["kappa_tau"]
    return xs_ZH(pois)*ktau*ktau/Gamma_tot(pois)

# Define gradient functions 
def dGamma_tot(pois,param):
    sumtotals = 0.58+0.22+0.08+0.06+0.027+0.029+0.0023+0.0016+0.00025+0.00022

    if param == "kappa_b": return 2*0.58*pois["kappa_b"]/sumtotals
    elif param == "kappa_W": return 2*0.22*pois["kappa_W"]/sumtotals
    elif param == "kappa_glu": return 2*0.08*pois["kappa_glu"]/sumtotals
    elif param == "kappa_tau": return 2*0.06*pois["kappa_tau"]/sumtotals
    elif param == "kappa_Z": return 2*0.027*pois["kappa_Z"]/sumtotals
    elif param == "kappa_c": return 2*0.029*pois["kappa_c"]/sumtotals
    elif param == "kappa_gamma": return 2*0.0023*pois["kappa_gamma"]/sumtotals
    elif param == "kappa_Zg": return 2*0.0016*pois["kappa_Zg"]/sumtotals
    elif param == "kappa_s": return 2*0.00025*pois["kappa_s"]/sumtotals
    elif param == "kappa_mu": return 2*0.00022*pois["kappa_mu"]/sumtotals
    else: return 0.0


def dxs_ZH(pois,param):
    kZ = pois["kappa_Z"]
    if param=="kappa_Z": return 2*kZ
    else: return 0.0 


def dXSBRgeneral(param,pois):
    k = pois[param]
    p1 = dxs_ZH(pois,param)*k*k/Gamma_tot(pois)
    p2 = xs_ZH(pois)*(2*k/Gamma_tot(pois))
    p3 = xs_ZH(pois)*k*k*(1./Gamma_tot(pois)**2)*dGamma_tot(pois,param)
    return p1 + p2 + p3


def dxsBR_hZZ(pois,param):
    if param == "kappa_Z":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hWW(pois,param):
    if param == "kappa_W":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hbb(pois,param):
    if param == "kappa_b":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hcc(pois,param):
    if param == "kappa_c":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hss(pois,param):
    if param == "kappa_s":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hgluglu(pois,param):
    if param == "kappa_glu":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hgamgam(pois,param):
    if param == "kappa_gamma":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hZgam(pois,param):
    if param == "kappa_Zg":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_hmumu(pois,param):
    if param == "kappa_mu":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

def dxsBR_htautau(pois,param):
    if param == "kappa_tau":
        return dXSBRgeneral(param, pois)
    else: 
        return 0.0

functions = od()
functions["ZH_had"]=xs_ZH
functions["ZH_lep"]=xs_ZH
functions["ZH_hbb"]=xsBR_hbb
functions["ZH_hcc"]=xsBR_hcc
functions["ZH_hss"]=xsBR_hss
functions["ZH_hgluglu"]=xsBR_hgluglu
functions["ZH_hgamgam"]=xsBR_hgamgam
functions["ZH_tautau"]=xsBR_htautau
functions["ZH_mumu"]=xsBR_hmumu
functions["ZH_hZZ"]=xsBR_hZZ
functions["ZH_hWW"]=xsBR_hWW
functions["ZH_hZgam"]=xsBR_hZgam


grad_functions = od()
grad_functions["ZH_had"]=dxs_ZH
grad_functions["ZH_lep"]=dxs_ZH
grad_functions["ZH_hbb"]=dxsBR_hbb
grad_functions["ZH_hcc"]=dxsBR_hcc
grad_functions["ZH_hss"]=dxsBR_hss
grad_functions["ZH_hgluglu"]=dxsBR_hgluglu
grad_functions["ZH_hgamgam"]=dxsBR_hgamgam
grad_functions["ZH_tautau"]=dxsBR_htautau
grad_functions["ZH_mumu"]=dxsBR_hmumu
grad_functions["ZH_hZZ"]=dxsBR_hZZ
grad_functions["ZH_hWW"]=dxsBR_hWW
grad_functions["ZH_hZgam"]=dxsBR_hZgam
