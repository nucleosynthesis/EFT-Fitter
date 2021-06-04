from collections import OrderedDict as od

name = "HIG19015_chi2"

def findCrossings(pts,c):
  # start from the left and find the point where we cross c
  # pts are poi,chi2
  xings =[]
  for i,pt in enumerate(pts):     
    if i==len(pts)-1: break
    pt2 = pts[i+1]
    #print(pt,pt2,c)
    #print(pt[1] < c and pt2[1] > c)
    #print(pt[1] > c and pt2[1] < c)
    if (pt[1] < c and pt2[1] > c) or (pt[1] > c and pt2[1] < c): 
      m = (pt2[1]-pt[1])/(pt2[0]-pt[0])
      b = pt[1] - m*pt[0]
      x = (c-b)/m 
      xings.append(x)
  return xings

def findMin(pts):
  pts2 = [[pt[1],pt[0]] for pt in pts]
  pts2.sort()
  return pts2[0][1]

def findMinnLL(pts):
  pts2 = [[pt[1],pt[0]] for pt in pts]
  pts2.sort()
  return pts2[0][0]

params = [
'0J_PTH_0_10_Tag0', '0J_PTH_0_10_Tag1', '0J_PTH_0_10_Tag2', '0J_PTH_GT10_Tag0', '0J_PTH_GT10_Tag1', '0J_PTH_GT10_Tag2', '1J_PTH_0_60_Tag0', '1J_PTH_0_60_Tag1', '1J_PTH_0_60_Tag2', '1J_PTH_120_200_Tag0', '1J_PTH_120_200_Tag1', '1J_PTH_120_200_Tag2', '1J_PTH_60_120_Tag0', '1J_PTH_60_120_Tag1', '1J_PTH_60_120_Tag2', 'GE2J_PTH_0_60_Tag0', 'GE2J_PTH_0_60_Tag1', 'GE2J_PTH_0_60_Tag2', 'GE2J_PTH_120_200_Tag0', 'GE2J_PTH_120_200_Tag1', 'GE2J_PTH_120_200_Tag2', 'GE2J_PTH_60_120_Tag0', 'GE2J_PTH_60_120_Tag1', 'GE2J_PTH_60_120_Tag2', 'PTH_200_300_Tag0', 'PTH_200_300_Tag1', 'PTH_300_450_Tag0', 'PTH_300_450_Tag1', 'PTH_450_650_Tag0', 'PTH_GT650_Tag0', 'THQ_LEP', 'TTH_HAD_PTH_0_60_Tag0', 'TTH_HAD_PTH_0_60_Tag1', 'TTH_HAD_PTH_0_60_Tag2', 'TTH_HAD_PTH_120_200_Tag0', 'TTH_HAD_PTH_120_200_Tag1', 'TTH_HAD_PTH_120_200_Tag2', 'TTH_HAD_PTH_120_200_Tag3', 'TTH_HAD_PTH_200_300_Tag0', 'TTH_HAD_PTH_200_300_Tag1', 'TTH_HAD_PTH_200_300_Tag2', 'TTH_HAD_PTH_60_120_Tag0', 'TTH_HAD_PTH_60_120_Tag1', 'TTH_HAD_PTH_60_120_Tag2', 'TTH_HAD_PTH_GT300_Tag0', 'TTH_HAD_PTH_GT300_Tag1', 'TTH_LEP_PTH_0_60_Tag0', 'TTH_LEP_PTH_0_60_Tag1', 'TTH_LEP_PTH_0_60_Tag2', 'TTH_LEP_PTH_120_200_Tag0', 'TTH_LEP_PTH_120_200_Tag1', 'TTH_LEP_PTH_200_300_Tag0', 'TTH_LEP_PTH_60_120_Tag0', 'TTH_LEP_PTH_60_120_Tag1', 'TTH_LEP_PTH_60_120_Tag2', 'TTH_LEP_PTH_GT300_Tag0', 'VBFLIKEGGH_Tag0', 'VBFLIKEGGH_Tag1', 'VBFTOPO_BSM_Tag0', 'VBFTOPO_BSM_Tag1', 'VBFTOPO_JET3VETO_HIGHMJJ_Tag0', 'VBFTOPO_JET3VETO_HIGHMJJ_Tag1', 'VBFTOPO_JET3VETO_LOWMJJ_Tag0', 'VBFTOPO_JET3VETO_LOWMJJ_Tag1', 'VBFTOPO_JET3_HIGHMJJ_Tag0', 'VBFTOPO_JET3_HIGHMJJ_Tag1', 'VBFTOPO_JET3_LOWMJJ_Tag0', 'VBFTOPO_JET3_LOWMJJ_Tag1', 'VBFTOPO_VHHAD_Tag0', 'VBFTOPO_VHHAD_Tag1', 'VH_MET_Tag0', 'VH_MET_Tag1', 'VH_MET_Tag2', 'WH_LEP_PTV_0_75_Tag0', 'WH_LEP_PTV_0_75_Tag1', 'WH_LEP_PTV_75_150_Tag0', 'WH_LEP_PTV_75_150_Tag1', 'WH_LEP_PTV_GT150_Tag0', 'ZH_LEP_Tag0', 'ZH_LEP_Tag1'
]

# Bestfit + uncertainties: if -sigma goes below zero then use expected
X = od()
for p in params:
  in_p="inputs/HIG_19_015/root/output_ic_crab_full_wotheory_RECO_"+p+"_combine.txt"
  fi = open(in_p,"r")
  points  = []
  for i,line in enumerate(fi.readlines()): 
   if i==0: par  = line.split()[0]
   else:
    #if float(line.split()[1]) > MAXC2-0.001: continue
    # we assume they should be from 0 so reset why not
    points.append([float(line.split()[0]),float(line.split()[1])])
 
  points.sort()
  mini = findMin(points)
  minNLL = findMinnLL(points)
  # need to shift by minimum value 
  points = [[p[0],p[1]-minNLL] for p in points]
  crossings1 = findCrossings(points,1)
  if len(crossings1)>1: 
   p1s  =crossings1[1]-mini
   m1s=mini-crossings1[0]
  else : 
   p1s=crossings1[0]-mini
   m1s  = mini

  X[p]        = {"bestfit":mini, "Up01Sigma":p1s,  "Down01Sigma":m1s, "Up01SigmaExp":0., "Down01SigmaExp":0., "merged":False}

