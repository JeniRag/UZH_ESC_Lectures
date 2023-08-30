# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 14:10:47 2021

@author: Sujeni

ESC202-exercise_L1
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Rectangle 

"""classes"""
class particle:
	def __init__(self,r): #r=location, self referencing to the class itself
		self.r=r
        
class cell:
    def __init__(self, rLow, rHigh, lower, upper):
        self.rLow = rLow
        self.rHigh = rHigh
        
        self.iLower = lower
        self.iUpper = upper
        
        self.pLower = None
        self.pUpper = None

        self.Node=True
        self.Leaf=False
   
"""functions"""

def partition( A, i, j, v, dim):
    ##Variant 1
    #"""
    l=i
    for k in range(i,j+1): 
        if A[k].r[dim]<v:
            A[k],A[l]=A[l],A[k]
            l+=1
    s=l
    return s
    #"""
    """
    #Variant2
    s=0
    l=j
    jnew=j
    for k in range(i,j):
        if A[k].r[dim]<v:
            for l in range(jnew,i,-1):
                if k<l and A[l].r[dim]>v:
                    A[k],A[l]=A[l],A[k]
                    jnew=l
                if k>=l:
                    s=k
                    break
    return s #start of the upper partition
    """
"""functions"""

def split_variables(A,root,dim):
    v = (root.rLow[dim]+root.rHigh[dim])/2 #split value
    s = partition(A, root.iLower, root.iUpper,  v, dim )
    return v, s

def cLower(A,root,dim,v,s,ax):
    rLow=np.copy(root.rLow)
    rHigh=np.copy(root.rHigh) 
    rHigh[dim]=v
    cLow = cell(rLow, rHigh, root.iLower, s-1)
    
    root.pLower=cLow
    return cLow

def cHigher(A,root,dim,v,s,ax):
    rHigh=np.copy(root.rHigh)     #make Copy
    rLow=np.copy(root.rLow)#lowest edge
    rLow[dim]=v #highest edge
    cHigh=cell(rLow,rHigh, s, root.iUpper)
   
    root.pUpper=cHigh #linkage to parental cell
    return cHigh

def treebuild(A,root,dim,tree_list,ax):
   
    v,s = split_variables(A,root,dim)
   
    print("s="+str(s)+" ,v="+str(v))
    
    plottree(root,ax)
    # May have two parts: lower.. s-1 and s... upper
    tree_list.append(root)
        #if there is a lower part:
    if s > root.iLower:

        cLow=cLower(A,root,dim,v,s,ax)
        #if there are more than 8 particles in cell:
        if len(A[root.iLower: s-1])>=value:
            
            treebuild(A, cLow, 1-dim,tree_list,ax)  
          
            plottree(cLow,ax)
            tree_list.append(cLow)
        else:
            root.Leaf=True
            root.Node=False
            plottree(cLow,ax)
            tree_list.append(cLow)
            
    #else if there is an upper part
    if s <= root.iUpper:
        
        cHigh=cHigher(A,root,dim,v,s,ax)
       
        #if there are more than 8 particles in cell:
       
        if len(A[s:root.iUpper])>=value:
            treebuild(A,cHigh, 1- dim,tree_list,ax)
            plottree(cHigh,ax)
            tree_list.append(cHigh)
        else:
            root.Leaf=True
            root.Node=False
            plottree(cHigh,ax)
            tree_list.append(cHigh)
            
#grafical representation of tree 
def plottree(root,ax):
    #draw a rectangle specified by rLow and rHigh
    # matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
   
    ax.add_patch( Rectangle(root.rLow, root.rHigh[0]-root.rLow[0], root.rHigh[1]-root.rLow[1],angle=0.0, fc='none',ec ='b', 
                        lw = 1))
                            
def plot_particle(array,ax):
    for i in array:
        ax.scatter(i.r[0], i.r[1], c='b',s=10)

"""array with particle coordinates"""
#Create array A with particles
n=200#number of particles
ntest1=20

array=np.zeros(n, dtype=particle)
for i in range(0, n):
    element=np.array([random.random(),random.random()])
    array[i]=particle(element)

#Test partition function with known solution
test_array1=np.zeros(ntest1,dtype=particle)
for i in range(0,ntest1):
    element=np.array([(i+1)/ntest1,(i+1)/ntest1])
    test_array1[i]=particle(element)
    
"""test"""
#array=test_array1

"""
for p in array:
   print(p.r)
"""
value=1
   
rLow=np.array([0.0,0.0])
rHigh=np.array([1.0,1.0])

lower = 0
upper = len(array)-1

root = cell(rLow, rHigh, lower, upper)

dim = 0

#p1=partition(array,0,len(array)-1,0.2,1)
#p0=partition(array,0,len(array)-1,0.2,0)

#print(p0)
#print(p1)

fig = plt.figure() 
ax = fig.add_subplot(111)
ax.set_title("exercise 1") 

tree_list=[]

treebuild(array,root, dim,tree_list,ax)
plot_particle(array,ax)

plt.savefig("Tree.png", dpi=300)   
 