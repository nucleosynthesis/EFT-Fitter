from dataclasses import dataclass
import numpy as np

@dataclass
class dataresult: 
	pvals: np.array 
	chi2 : np.array 
	allpvals: np.array 
	allparams: list 
	allpredictions: list 