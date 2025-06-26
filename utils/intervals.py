MAXC2 = 10

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

def mymin(l):
  if len(l)==0: return 0
  else: return min(l)