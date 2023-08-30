# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 14:20:49 2021

@author: Sujeni

ESC202-exercise 2
"""
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.patches import Rectangle 
from timeit import default_timer as timer

"""classes"""

class particle:
    def __init__(self, r):
        self.r=r
        self.count=0
        self.color="mediumturquoise"
          
class cell:
    def __init__(self, rLow, rHigh, lower, upper):
        self.rLow = rLow
        self.rHigh = rHigh
        
        self.iLower = lower
        self.iUpper = upper
        
        self.pLower = None
        self.pUpper = None

        self.Node=True
        self.isLeaf=False
        
        self.c= rLow+(rHigh-rLow)/2 #center
        self.b=self.c-rLow#abs(rHigh-self.c) #diagonal from center to edge
        
        self.dist2=0
        self.used=False
   
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
            root.isLeaf=True
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
            root.isLeaf=True
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
        ax.scatter(i.r[0], i.r[1], c=i.color,s=10)
        
""""""
def euc_dist(r,a): #euclidean distance
    d2=(r[0]-a[0])**2+(r[1]-a[1])**2 #squared
    return d2

def mean_position(A,root):
    sum=np.zeros(2)
    length=0
    for a in A[root.iLower: root.iUpper+1]:
        sum+=a.r
        length+=1
   
    root.c=sum/length
        
def d2(A,root, p):
    
    mean_position(A,root)
    c=root.c #center of cell
    #b=abs(root.iUpper-c)
    b=root.b
    d2=0 # 
    for d in range(2):
        t=abs(c[d]-p[d])-b[d]
        if t>0:
            d2+=t**2
            
    root.dist2=d2
    return d2

def Ballwalk (A,c,p,rmax):
    
    r2max=rmax**2
    p.count=0
 
    if c.isLeaf==True:
        
        for a in A[c.iLower: c.iUpper+1]:      
            if euc_dist(p.r, a.r)<r2max:           
                p.count+=1
                #a.color="yellowgreen"

    else:
        
        if c.pLower!= None:
            if d2(A,c.pLower,p.r)<r2max:
               
               p.count+=Ballwalk(A,c.pLower,p, rmax)
              
        if c.pUpper != None:
            if d2(A,c.pUpper,p.r)<r2max: ##right
               p.count+= Ballwalk(A,c.pUpper,p,rmax)
                  
    return p.count

def plot_circle(p,ax,rmax):
    circle=plt.Circle((p.r[0],p.r[1]), rmax, alpha=0.2)
    ax.add_patch(circle)
"""    
def method2(c,p,L,d):
    d2=0
    b=c.b
    for d in range(0,d2):
        t=c[d]-p.r[d]
        if t<0:
            t1=t+L
        else:
            t1=t-L
        t=abs(t)-b[d]
        t1=abs(t1)-b[d]
        if t1<t:
            t=t1
        if t>0:
            d2+=t*t
    return d2
"""
def final_plots(N,ax,rmax):
    
    array=np.zeros(N, dtype=particle)
    for i in range(0, N):
        element=np.array([random.random(),random.random()])
        array[i]=particle(element)
    
    point=np.array([0.6,0.5])
    array[0]=particle(point)
    p1=array[0]
    #p1=array[int(random.random()*N)]
    
    
    rLow=np.array([0.0,0.0])
    rHigh=np.array([1.0,1.0])
    
    lower = 0
    upper = len(array)-1
    
    root = cell(rLow, rHigh, lower, upper)
    
    dim = 0

    tree_list=[]

    treebuild(array,root, dim,tree_list,ax)
    
    start=timer()
    for p in array:
        count=Ballwalk (array,root,p,rmax)
        p.count=count
        #print(str(p.r)+str(count))
    end=timer()
    time=round(end-start,4)
    plt.show() 
    p1.color="purple"

    print("N=" +str(N)+", time="+ str(time) +", particle: "+ str(np.round(p1.r,2))+", count= "+str(count))
    plot_circle(p1,ax,rmax)
    plot_particle(array,ax)
    ax.set_title("Exercise 2") 
    ax.text(0, -0.04, "N=" +str(N)+", time="+ str(time) + ", particle: "+ str(np.round(p1.r,2))+", count= "+str(count))

""""""
value=8
N_list=[50,100,200,400,500]
n=len(N_list)

fig=[]
ax=[]

for x in range(n):
    fig.append("fig{0}".format(x))
    ax.append("ax{0}".format(x))

    fig[x]=plt.figure()
    ax[x]=fig[x].add_subplot(111)
    ax[x].set_title("exercise 2")

for i in range(n):
    final_plots(N_list[i],ax[i], 0.2)

plt.savefig("./Particle.png", dpi=300)

"""
The time increases, if the particle number is increased (if I count the neighbors for all particles).
I used method 1 as the intersection test, because
I don't need to apply a condition for rmax like in method 2.
"""