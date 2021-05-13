# Define functions
def extractTerms(_func,multiplier=1):
  terms = od()
  terms['const'] = 0
  f = re.sub(" ","",_func)
  f = re.sub("-","+-",f)
  for t in f.split("+"):
    c = t.split("*")
    if len(c) == 1: terms['const'] = float(c[0])*multiplier
    elif len(c) == 2: terms['A_%s'%c[1]] = float(c[0])*multiplier
    else:
      if c[1] == c[2]: terms['B_%s'%c[1]] = float(c[0])*multiplier
      else: terms['B_%s_%s'%(c[1],c[2])] = float(c[0])*multiplier
  return terms

# Terms to function
def termsToFunction(_terms):
  f = ""
  for k,v in _terms.items():
    if k == "const": f += "%.4f"%v
    elif k[0] == "A": f += "+%.4f*%s"%(v,k.split("_")[-1])
    elif k[0] == "B": 
      if len(k.split("_"))==2: f += "+%.4f*%s*%s"%(v,k.split("_")[-1],k.split("_")[-1])
      else: f += "+%.4f*%s*%s"%(v,k.split("_")[1],k.split("_")[2])
  f = re.sub("\+-","-",f)
  return f

# Function for printing matrix
def printMatrix(_mat):
  for i in range(len(_mat[0])):
    vstr =  "("
    for j in range(len(_mat[0])): vstr += " %-5.2f "%_mat[i][j]
    vstr += ")"
    print(vstr)