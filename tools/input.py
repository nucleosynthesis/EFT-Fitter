import numpy as np
import array
from scipy import linalg

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract all required info for input measurement
class INPUT:

  def __init__(self,inputMeasurement,functions,doAsimov):

    # Name
    self.name = inputMeasurement['name']

    # X0
    self.X0 = []
    self.XList = []
    
    if (len(inputMeasurement['X'].keys()))>1:
      self.type = "measurement"
      for x, vals in inputMeasurement['X'].items(): 
        if doAsimov: self.X0.append(1.)
        else: self.X0.append(float(vals['bestfit']))
        self.XList.append(x)
      self.X0 = np.asarray(self.X0)
    else:
      self.type = "spline"
      self.X0 = inputMeasurement['X']["likelihood"]['likelihood']
      self.XList = self.X0.getParameters()
      inputMeasurement['X']['likelihood']["merged"] = False
      for k in self.XList  : 
        inputMeasurement['X'][k] = {"likelihood":1} # could use this better
        inputMeasurement['X'][k]["merged"] = False

    # Scaling functions
    self.ProdScaling = []
    self.DecayScaling = []
    for x, vals in inputMeasurement["X"].items():
      if x == "likelihood": continue
      # first check if the function is already provided 
      if x in functions.keys(): 
       continue
      # Then appeal to prod*decay
      print("Having to rely on scaling function logic!")
      # Extract production and decay strings
      prod = "_".join(x.split("_")[:-1])
      dec = x.split("_")[-1]

      # For merged STXS bins
      if vals['merged']:
        terms = od()
        for stxs_bin,frac in vals['STXS_fractions'].items():
          _terms = extractTerms(functions[stxs_bin],multiplier=frac)
          for k,v in _terms.items():
            if k in terms: terms[k] = terms[k]+v
            else: terms[k] = v
        self.ProdScaling.append(terms)
        # Add function to list
        functions[prod] = termsToFunction(terms)

      # For unmerged STXS bin/process
      else: 
        self.ProdScaling.append( extractTerms(functions[prod]) )

      # Add list of partial width and total width scaling
      self.DecayScaling.append( [extractTerms(functions[dec]),extractTerms(functions['tot'])] )

    if self.type=="spline": 
      nbins = len( self.XList )
      self.nbins = nbins
      return

    # Error matrix
    corr = []
    print(inputMeasurement)
    for ix in inputMeasurement['X']:
      for jx in inputMeasurement['X']:
        if (ix,jx) in inputMeasurement['rho'].keys(): p = inputMeasurement['rho'][(ix,jx)]
        elif (jx,ix) in inputMeasurement['rho'].keys(): p = inputMeasurement['rho'][(jx,ix)]
        else:
          if ix == jx: p = 1.
          else: 
            p = 0.
            print(" --> [WARNING] No correlation info given for (%s,%s). Assuming = 0"%(ix,jx))
        corr.append(p)
    corr = array.array('d',corr)

    # Extract symmetrized errors
    # FIXME: add something more appropriate for observed!
    err = []
    for x,vals in inputMeasurement['X'].items():
      if doAsimov: err.append( 0.5*(vals['Up01SigmaExp']+vals['Down01SigmaExp']) )
      else: err.append( 0.5*(vals['Up01Sigma']+vals['Down01Sigma']) )
    
    # Do some squarification and inverting
    nbins = len( inputMeasurement['X'].keys() )
    corr_sq = [ corr[i:i+nbins] for i in range(0,len(corr),nbins)]
    cov_sq = [ [ corr_sq[i][j]*err[i]*err[j] for i in range(nbins)] for j in range(nbins) ]
    self.err_mat = np.array(cov_sq)
    self.Vinv = linalg.inv(self.err_mat)
    self.nbins = nbins
