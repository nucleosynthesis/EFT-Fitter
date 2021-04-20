from collections import OrderedDict as od
import sys
# LHC-HXSWG kappa functions, taken from HIG-17-031 
def Gamma_hgg(pois):
   kV = pois["kappa_V"]
   kF = pois["kappa_F"]
   return 1.59*kV*kV + 0.07*kF*kF - 0.67*kV*kF

def Gamma_gluglu(pois):
   kF = pois["kappa_F"]
   return 1.04*kF*kF + 0.002*kF*kF - 0.038*kF*kF

def Gamma_tot(pois): 
   kV = pois["kappa_V"]
   kF = pois["kappa_F"]
   kglu2 = Gamma_gluglu(pois)
   kgam2 = Gamma_hgg(pois)
   kZg2  = kgam2
   sumtotals = 0.58+0.22+0.08+0.06+0.026+0.029+0.0023+0.0015+0.00025+0.00022
   kH2 = 0.58*kF*kF + 0.22*kV*kV + 0.08*kglu2 +0.06*kF*kF +0.026*kV*kV +0.029*kF*kF +0.0023*kgam2 +0.0015*kZg2 + 0.00025*kF*kF+0.00022*kF*kF
   return kH2/sumtotals

def BR_hgg(pois):
   return Gamma_hgg(pois)/Gamma_tot(pois) 

def ggZH(pois):
   kV = pois["kappa_V"]
   kF = pois["kappa_F"]
   return (2.46*kV*kV+0.47*kF*kF-1.94*kF*kV)/(2.46+0.47-1.94)

def ZH(pois):
   kV = pois["kappa_V"]
   return kV*kV

def WH(pois):
   kV = pois["kappa_V"]
   return kV*kV

def qqH(pois):
   kV = pois["kappa_V"]
   return kV*kV

def ggH(pois):
   kF = pois["kappa_F"]
   return kF*kF

def ttH(pois):
   kF = pois["kappa_F"]
   return kF*kF

def tHq(pois):
   kV = pois["kappa_V"]
   kF = pois["kappa_F"]
   return (2.63*kF*kF + 3.58*kV*kV - 5.21*kF*kV)/(2.63+3.58-5.21)

def tHW(pois):
   kV = pois["kappa_V"]
   kF = pois["kappa_F"]
   return 2.91*kF*kF + 2.40*kV*kV - 4.22*kF*kV/(2.91+2.40-4.22)

# functions I need map the yields of the contributing processes onto the reco yields

# Leptonic categories
def THQ_LEP(pois):
    cttH = 1.+4.+7.+8.+9.
    ctHq = 24.
    ctHW = 6. 
    cWH  = 7.+10.+9.
    cZH  = 5.
    cggZH = 3. 
    cqqH = 4.
    sumc = cttH+ctHq+ctHW+cWH+cggZH+cZH+cqqH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggZH*ggZH(pois)+cWH*WH(pois)+cZH*ZH(pois)+cqqH*qqH(pois))
    return s*BR_hgg(pois)


def TTH_LEP_PTH_GT300_Tag0(pois):
    cttH = 7.+62.
    ctHq = 3.
    ctHW = 9. 
    cWH  = 1+2.+10.
    cZH  = 2.
    cggZH = 3. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW+cWH+cggZH+cZH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggZH*ggZH(pois)+cWH*WH(pois)+cZH*ZH(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_200_300_Tag0(pois):
    cttH = 1.+86.+1.
    ctHq = 3.
    ctHW = 5. 
    cWH  = 1.+1.
    cZH  = 1.
    cggZH = 1. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW+cWH+cggZH+cZH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggZH*ggZH(pois)+cWH*WH(pois)+cZH*ZH(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_120_200_Tag1(pois):
    cttH = 1.+78.+1.
    ctHq = 3.
    ctHW = 2. 
    cWH  = 2.+2.+5.
    cZH  = 2.
    cggZH = 1. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW+cWH+cggZH+cZH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggZH*ggZH(pois)+cWH*WH(pois)+cZH*ZH(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_120_200_Tag0(pois):
    cttH = 1.+90.+2.
    ctHq = 2.
    ctHW = 2. 
    cWH  = 1.+1.
    cZH  = 0.
    cggZH = 0. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW+cWH+cggZH+cZH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggZH*ggZH(pois)+cWH*WH(pois)+cZH*ZH(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_60_120_Tag2(pois):
    cttH = 1.+91.+2.
    ctHq = 3.
    ctHW = 2. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois))
    return s*BR_hgg(pois)


def TTH_LEP_PTH_60_120_Tag1(pois):
    cttH = 1.+91.+3.+1.
    ctHq = 2.
    ctHW = 1. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_60_120_Tag0(pois):
    cttH = 95.+2.
    ctHq = 1.
    ctHW = 1. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_0_60_Tag2(pois):
    cttH = 88.+2.
    ctHq = 3.
    ctHW = 1. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_0_60_Tag1(pois):
    cttH = 94.+3.+1.
    ctHq = 2.
    ctHW = 1. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois))
    return s*BR_hgg(pois)

def TTH_LEP_PTH_0_60_Tag0(pois):
    cttH = 94.+4.
    ctHq = 1.
    ctHW = 0. 
    #cggH = 0.
    sumc = cttH+ctHq+ctHW
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois))
    return s*BR_hgg(pois)

# Hadronic categories
def TTH_HAD_PTH_GT300_Tag1(pois):
    cggH = 20.+4.+2.
    cttH = 46.
    ctHq = 11.
    ctHW = 7.
    cqqH = 7.
    sumc = cttH+ctHq+ctHW+cggH+cqqH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois)+cqqH*qqH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_GT300_Tag0(pois):
    cggH = 6.+1.+1.
    cttH = 1.+74.
    ctHq = 8.
    ctHW = 7. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_200_300_Tag2(pois):
    cggH = 19.
    cttH = 57.+1.
    ctHq = 11.
    ctHW = 4. 
    cqqH = 9.
    sumc = cttH+ctHq+ctHW+cggH+cqqH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois)+cqqH*qqH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_200_300_Tag1(pois):
    cggH = 9.
    cttH = 2.+75.+1.
    ctHq = 7.
    ctHW = 3. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_200_300_Tag0(pois):
    cggH = 0.
    cttH = 2.+90.
    ctHq = 4.
    ctHW = 3. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_120_200_Tag3(pois):
    cggH = 12.+3.
    cttH = 1.+62.+1.
    ctHq = 8.
    ctHW = 2. 
    cqqH = 9.
    sumc = cttH+ctHq+ctHW+cggH+cqqH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois)+cqqH*qqH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_120_200_Tag2(pois):
    cggH = 7.+3.
    cttH = 1.+74.+1.
    ctHq = 6.
    ctHW = 2. 
    cqqH = 5.
    sumc = cttH+ctHq+ctHW+cggH+cqqH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois)+cqqH*qqH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_120_200_Tag1(pois):
    cggH = 3.+1.
    cttH = 2.+83.+1.
    ctHq = 4.
    ctHW = 2. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)


def TTH_HAD_PTH_120_200_Tag0(pois):
    cggH = 1.+1.
    cttH = 1.+91.+2.
    ctHq = 2.
    ctHW = 1. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_60_120_Tag2(pois):
    cggH = 1.+1.
    cttH = 1.+89.+2.
    ctHq = 3.
    ctHW = 1. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    #print("inside TTH_HAD_PTH_60_120_Tag2",s,ttH(pois),tHq(pois),tHW(pois),ggH(pois),BR_hgg(pois))

    return s*BR_hgg(pois)

def TTH_HAD_PTH_60_120_Tag1(pois):
    cggH = 5.
    cttH = 1.+91.
    ctHq = 2.
    ctHW = 1. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_60_120_Tag0(pois):
    cggH = 1.
    cttH = 1.+93.+4.
    ctHq = 1.
    ctHW = 1. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_0_60_Tag2(pois):
    cggH = 1.
    cttH = 90.+2.+1
    ctHq = 3.
    ctHW = 1. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

def TTH_HAD_PTH_0_60_Tag1(pois):
    cggH = 0.
    cttH = 93.+3.
    ctHq = 2.
    ctHW = 1. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)


def TTH_HAD_PTH_0_60_Tag0(pois):
    cggH = 2.
    cttH = 94.+2.
    ctHq = 1.
    ctHW = 0. 
    sumc = cttH+ctHq+ctHW+cggH
    s = (1./sumc)*(cttH*ttH(pois)+ctHq*tHq(pois)+ctHW*tHW(pois)+cggH*ggH(pois))
    return s*BR_hgg(pois)

functions = od()
functions["TTH_HAD_PTH_60_120_Tag2"]=TTH_HAD_PTH_60_120_Tag2
functions["THQ_LEP"]=THQ_LEP
functions["TTH_HAD_PTH_GT300_Tag0"]=TTH_HAD_PTH_GT300_Tag0
functions["TTH_HAD_PTH_0_60_Tag0"]=TTH_HAD_PTH_0_60_Tag0
functions["TTH_HAD_PTH_GT300_Tag1"]=TTH_HAD_PTH_GT300_Tag1 
functions["TTH_HAD_PTH_0_60_Tag1"]=TTH_HAD_PTH_0_60_Tag1   
functions["TTH_LEP_PTH_0_60_Tag0"]=TTH_LEP_PTH_0_60_Tag0
functions["TTH_HAD_PTH_0_60_Tag2"]=TTH_HAD_PTH_0_60_Tag2    
functions["TTH_LEP_PTH_0_60_Tag1"]=TTH_LEP_PTH_0_60_Tag1 
functions["TTH_HAD_PTH_120_200_Tag0"]=TTH_HAD_PTH_120_200_Tag0 
functions["TTH_LEP_PTH_0_60_Tag2"]=TTH_LEP_PTH_0_60_Tag2 
functions["TTH_HAD_PTH_120_200_Tag1"]=TTH_HAD_PTH_120_200_Tag1 
functions["TTH_LEP_PTH_120_200_Tag0"]=TTH_LEP_PTH_120_200_Tag0 
functions["TTH_HAD_PTH_120_200_Tag2"]=TTH_HAD_PTH_120_200_Tag2 
functions["TTH_LEP_PTH_120_200_Tag1"]=TTH_LEP_PTH_120_200_Tag1 
functions["TTH_HAD_PTH_120_200_Tag3"]=TTH_HAD_PTH_120_200_Tag3 
functions["TTH_LEP_PTH_200_300_Tag0"]=TTH_LEP_PTH_200_300_Tag0 
functions["TTH_HAD_PTH_200_300_Tag0"]=TTH_HAD_PTH_200_300_Tag0 
functions["TTH_LEP_PTH_60_120_Tag0"]=TTH_LEP_PTH_60_120_Tag0 
functions["TTH_HAD_PTH_200_300_Tag1"]=TTH_HAD_PTH_200_300_Tag1 
functions["TTH_LEP_PTH_60_120_Tag1"]=TTH_LEP_PTH_60_120_Tag1 
functions["TTH_HAD_PTH_200_300_Tag2"]=TTH_HAD_PTH_200_300_Tag2 
functions["TTH_LEP_PTH_60_120_Tag2"]=TTH_LEP_PTH_60_120_Tag2 
functions["TTH_HAD_PTH_60_120_Tag0"]=TTH_HAD_PTH_60_120_Tag0  
functions["TTH_LEP_PTH_GT300_Tag0"]=TTH_LEP_PTH_GT300_Tag0 
functions["TTH_HAD_PTH_60_120_Tag1"]=TTH_HAD_PTH_60_120_Tag1 
