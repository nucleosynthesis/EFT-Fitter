# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np
verbose=False
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chi2: function to minimize
#  * Takes as input array of POIS to fit, if list is empty then will not reset
def GetChi2(P,args=[]):#FIT):
  #print(args)
  FIT=args
  chi2 = 0.

  # Set POIS for minimizer
  if len(P) != 0:
    FIT.setPOIS({ipoi:P[ip] for ip,ipoi in enumerate(FIT.PToFitList)})
    #for ip,ipoi in enumerate(FIT.PToFitList): 
    #  FIT.setPOIS({ipoi:P[ip]})

  init_i = len(FIT.PToFitList)
  nuisance_parameter_vals = [P[i] for i  in range(init_i,len(P))]
  # Set nuisance parameters -> need to have the FIT object keep track of which one is which, but once its 
  # decided it will be fixed and can be ignored 
  #print("nuisance param vals = ",nuisance_parameter_vals)
  if FIT.has_uncerts: 
    #for iN in range(len(nuisance_parameter_vals)) : 
    FIT.setNuisances({iN:nuisance_parameter_vals[iN] for iN in range(len(nuisance_parameter_vals))})
    nuisance_parameters = np.asarray(nuisance_parameter_vals)
    nuisance_parametersT = nuisance_parameters.T
    chi2 += nuisance_parametersT.dot(FIT.tHVinv.dot(nuisance_parameters))
  # Calculate chi2 terms
  for _input in FIT.INPUTS:
    X0 = _input.X0
    # first thing should return 0 if there's no systematics. 
    #X = np.asarray( [ (1+FIT.evaluateTHUncertainty(_input.XList[i]))*FIT.evaluateScalingFunctions(_input.ProdScaling[i])*(FIT.evaluateScalingFunctions(_input.DecayScaling[i][0])/FIT.evaluateScalingFunctions(_input.DecayScaling[i][1])) for i in range(_input.nbins)] )
    X = np.asarray( [ (1+FIT.evaluateTHUncertainty(_input.XList[i]))*FIT.evaluateScalingFunctions(_input.XList[i]) for i in range(_input.nbins)] )
    if _input.type == "spline":
          #chi2 += sum([ (X[j]-1)/1 for j in range(len(X)) ])
          chi2 += X0.evaluate({"%s"%_input.XList[j]:X[j] for j in range(len(X))})
    else: 
          xarr = X-X0
          xarrT = xarr.T
          chi2 += xarrT.dot( _input.Vinv.dot( xarr ) )

  if verbose:
    print(" ----------------------------------------------------------------------------")
    print(" --> Scaling functions: Linear" if FIT.linearOnly else " --> Scaling functions: Quadratic")
    print("")
    print(" --> Parameters of interest:")
    for ip,ipoi in enumerate(FIT.PList): 
     if ipoi in FIT.PToFitList: print("   * %s = %.3f"%(ipoi,FIT.P0[ip]))
     else: print("   * %s = %.3f (fixed)"%(ipoi,FIT.P0[ip]))

    for _input in FIT.INPUTS:
      print("\n --> Inputs: %s"%_input.name)
      if _input.type == "spline" : continue
      X0 = _input.X0
      #xvals_list = [ FIT.evaluateScalingFunctions(_input.ProdScaling[i])*(FIT.evaluateScalingFunctions(_input.DecayScaling[i][0])/FIT.evaluateScalingFunctions(_input.DecayScaling[i][1])) for i in range(_input.nbins)]
      xvals_list = [ FIT.evaluateScalingFunctions(_input.XList[i]) for i in range(_input.nbins)]
      X = np.asarray(xvals_list)
      for ix in range(_input.nbins):
        print("   * %-50s : X0 = %.6f, X(p) = %.6f"%(_input.XList[ix],X0[ix],X[ix]))
      print("   * Vij = ")
      for i in range(len(_input.Vinv[0])):
        vstr =  "           ("
        for j in range(len(_input.Vinv[0])): vstr += " %-5.2f "%_input.Vinv[i][j]
        vstr += ")"
        print(vstr)

    print("\n --> Result: chi2 = %.8f"%chi2)
      
  return chi2