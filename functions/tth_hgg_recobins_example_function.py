from collections import OrderedDict as od

params = [
"TTH_HAD_PTH_60_120_Tag2"  
, "THQ_LEP"                  
, "TTH_HAD_PTH_GT300_Tag0"  
, "TTH_HAD_PTH_0_60_Tag0"    
, "TTH_HAD_PTH_GT300_Tag1" 
, "TTH_HAD_PTH_0_60_Tag1"    
, "TTH_LEP_PTH_0_60_Tag0" 
, "TTH_HAD_PTH_0_60_Tag2"    
, "TTH_LEP_PTH_0_60_Tag1" 
, "TTH_HAD_PTH_120_200_Tag0" 
, "TTH_LEP_PTH_0_60_Tag2" 
, "TTH_HAD_PTH_120_200_Tag1" 
, "TTH_LEP_PTH_120_200_Tag0" 
, "TTH_HAD_PTH_120_200_Tag2" 
, "TTH_LEP_PTH_120_200_Tag1" 
, "TTH_HAD_PTH_120_200_Tag3" 
, "TTH_LEP_PTH_200_300_Tag0" 
, "TTH_HAD_PTH_200_300_Tag0" 
, "TTH_LEP_PTH_60_120_Tag0" 
, "TTH_HAD_PTH_200_300_Tag1" 
, "TTH_LEP_PTH_60_120_Tag1" 
, "TTH_HAD_PTH_200_300_Tag2" 
, "TTH_LEP_PTH_60_120_Tag2" 
, "TTH_HAD_PTH_60_120_Tag0"  
, "TTH_LEP_PTH_GT300_Tag0" 
, "TTH_HAD_PTH_60_120_Tag1"
]

# nothing but a super weird way to have a function that returns the value 
# given in the dictionary provided, corresponding to a particular name
class functionDirectory:
  def __init__(self,name): 
    self.name = name
  def createFunction(self,p):
    def a_function(pois={p:1.0}):
     lookup = str(p)
     return pois[lookup] 
    return a_function
  def addfunction(self,name):
    x = self.createFunction(name)
    setattr(self,name,x)
  def getfunction(self,name):
    return getattr(self,name)

myfuncs = functionDirectory("tth_hgg_functions")

for p in params: 
  myfuncs.addfunction(p)

functions = od()
for p in params: 
  functions[p] = myfuncs.getfunction(p)

#print(functions)
#print(functions["TTH_HAD_PTH_60_120_Tag1"]({"TTH_HAD_PTH_60_120_Tag1":3.8}))
#print(functions["TTH_LEP_PTH_200_300_Tag0"]({"TTH_LEP_PTH_60_120_Tag2":9.6,"TTH_LEP_PTH_200_300_Tag0":999.9 }))


