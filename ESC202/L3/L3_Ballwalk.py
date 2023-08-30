# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 14:00:41 2021

@author: Sujeni

ESC202- exericse 3
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib as mpl
import matplotlib.cm as cm
import random
from matplotlib.patches import Rectangle 
from timeit import default_timer as timer
from heapq import *

"""classes"""

class particle:
    def __init__(self, r):
        self.r=r
        self.count=0
        self.color="mediumturquoise"
        self.rho=0.0
        self.m=1 
        self.rho=0 #initial
    
          
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

class PRIQ:
    def __init__(self, n):
        self.heap=[]
        sentinel=(0,None)
        
        for i in range(n):
            heappush(self.heap, sentinel)
        
    def show(self):
        return list(self.heap)
    
    def dist2(self): #heap.dist2()
       
        if self.heap[0][0]==0:
            d2=float('inf')
        else:
            d2=self.heap[0][0] ####heappop?
            d2=1/d2
        return d2

    def replace(self, dist,p):
        if dist==0:
            heapreplace(self.heap, (float('inf'), p))
        else:
            heapreplace(self.heap, ((1.0/dist), p))
        
        
"""treebuild functions"""

def partition( A, i, j, v, dim):
  
    l=i
    for k in range(i,j+1): 
        if A[k].r[dim]<v:
            A[k],A[l]=A[l],A[k]
            l+=1
    s=l
    return s
  
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


""""""
def euc_dist(r,a): #euclidean distance
    d2=(r[0]-a[0])**2+(r[1]-a[1])**2 #squared
    return d2
        
def d2(A,root, p):

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

""" density"""
def mass(queue):
    heap=queue.show()
    ma=0
    for i in range(len(heap)):
        ma+=heap[i][1].m
    return ma

def density(p,queue):
    r=np.sqrt(queue.dist2())
    
    m_sum=mass(queue) 
    rho=m_sum/(2*np.pi*r**2)
    p.rho=rho
    
    return rho
"""neighboursearch"""
def Ballwalk(A,c,p,queue,rOffset):
    pi=p.r+rOffset
    if c.isLeaf==True:
         
        for a in A[c.iLower: c.iUpper+1]: 
            dist2=euc_dist(pi,a.r)
                 
            if dist2<queue.dist2():  
                queue.replace(dist2, a)
                a.color="orange"
                
    else:
       
        if c.pUpper !=None and c.pLower!=None: 
            #calculate distances once and use results later in order to save some time
            d2_pLower = d2(A,c.pLower,p.r)
            d2_pUpper = d2(A,c.pUpper,p.r)
            if d2_pLower<d2_pUpper: 
                #here we only look at the cell if it is close enough
                if d2_pLower < queue.dist2():
                   Ballwalk(A,c.pLower,p, queue, rOffset)
                if d2_pUpper < queue.dist2():
                   Ballwalk(A,c.pUpper,p,queue,rOffset)
        
            else:
                #same here as before
                if d2_pUpper < queue.dist2():
                  Ballwalk(A,c.pUpper,p, queue,rOffset)
                if d2_pLower < queue.dist2():
                  Ballwalk(A,c.pLower,p, queue,rOffset)
                   
        elif c.pUpper !=None and c.pLower==None:
          
           Ballwalk(A,c.pUpper,p, queue, rOffset)
            
        elif c.pLower !=None and c.pUpper==None:
           
           Ballwalk(A,c.pLower,p,queue, rOffset)

    return queue 

def Ballwalk_periodic(A,c,p,queue,period):
    for y in [0.0, -period[1], period[1]]:
        for x in [0.0, -period[0], period[0]]:
            rOffset=np.array([x,y])
            queue=Ballwalk(A,c,p,queue,rOffset)
    return queue
            
def extract_dist(p,array,heap):
    heap=heap.show()
    d2_list=[]
    for i in range(len(heap)):
        d2_list.append(1/heap[i][0])
    return d2_list

"""plot functions"""
#grafical representation of tree 
def plottree(root,ax):
    #draw a rectangle specified by rLow and rHigh
    # matplotlib.patches.Rectangle(xy, width, height, angle=0.0, **kwargs)
    ax.add_patch( Rectangle(root.rLow, root.rHigh[0]-root.rLow[0], root.rHigh[1]-root.rLow[1],angle=0.0, fc='none',ec ='black', 
                        lw = 1))
                            
def plot_particle(array,ax):
    
    for i in array:
        ax.scatter(i.r[0], i.r[1], c=i.color,s=20)
        
def densities_color(root,array,n,ax):
    queue=PRIQ(n)
    rho_list= np.zeros(len(array))
    for p in array:
        heap=Ballwalk(array,root,p,queue)
        dens=density(p,heap)
        rho_list[i] = dens
    norm=matplotlib.colors.Normalize(dmin=min(rho_list), dmax=max(rho_list))
    
    for p in array:
        p.color=cmap(norm(p.rho))    
                                                                
def plot_circle(p,ax,rmax):
    circle=plt.Circle((p.r[0],p.r[1]), rmax, alpha=0.2)
    ax.add_patch(circle)
    
def final_plot(root,array,n,period,ax,fig):
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    rho_list= np.zeros(len(array))

    for i in range(len(array)):
        x[i] = array[i].r[0]
        y[i] = array[i].r[1]
        queue=PRIQ(n)
        p=array[i]
        heap=Ballwalk_periodic(array,root,p,queue,period)
        density(p,heap)
        rho_list[i] = p.rho
    
 
    norm=matplotlib.colors.Normalize(vmin=min(rho_list), vmax=max(rho_list))
   
    cmap = cm.hot
    ax.scatter(x, y, s = 25, cmap = cmap, c = rho_list, alpha = 0.5)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="density")
                     
    #sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    #for i in range(len(array)):
    #    array[i].color= sm.to_rgba(rho_list[i])
       
        #cmap(norm(p.rho))
    #plot_particle(array,ax)    
    #fig.colorbar(sm)    
        
"""Test Code"""
# O(N**2) Test Code
# k = Number of nearest neighbors

def test_code(p,A,k):
    p.color="black"
    
    NN = []
    d2NN = []
    for i in range(k):
        d2min = float('inf')
       #for q in range(len(A)):
        for q in A:
            if p != q and q not in NN: 
                #d2 = p.dist2(q)
                d2=euc_dist(p.r,q.r)
                if d2 < d2min:
                    d2min = d2
                    qmin = q
        NN.append(qmin)
       
        d2NN.append(d2min)
    for i in NN:
        i.color="red"        
                
    return NN, d2NN
    # Here NN and d2NN lists for particle p are filled.
    # Compare them with the lists you got from the recursive algorithm

   
"""Parameters & Settings"""
N=1000
period=np.array([1,1])
array=np.zeros(N, dtype=particle)
for i in range(0, N):
    element=np.array([random.random(),random.random()])
    array[i]=particle(element)

    #p1=array[int(random.random()*N)]
    
value=8

cmap = plt.cm.rainbow

p=array[30]
p.color="black"

rLow=np.array([0.0,0.0])
rHigh=np.array([1.0,1.0])

lower = 0
upper = len(array)-1

root = cell(rLow, rHigh, lower, upper)

dim = 0
rmax=0
#final_plots(100,array,p, ax1,rmax)
tree_list=[]

#treebuild(array,root, dim,tree_list,ax)

""" Testing """

fig1=plt.figure()
ax1=fig1.add_subplot(111)

fig2=plt.figure()
ax2=fig2.add_subplot(111)

treebuild(array,root,dim,tree_list,ax1)
treebuild(array,root,dim,tree_list,ax2)

NN, d2NN=test_code(p,array,32)
plot_particle(array,ax2)
ax2.set_title("exercise 3- testcode")
for i in array:
    i.color="mediumturquoise"

queue=PRIQ(32)
heap1=Ballwalk_periodic(array,root,p,queue, period)
d2_list=extract_dist(p,array,heap1)
p.color="black"
plot_particle(array,ax1)
ax1.set_title("exercise 3- Ballwalk boundary")

"""density plot """
fig3=plt.figure()
ax3=fig3.add_subplot(111)
treebuild(array,root,dim,tree_list,ax3)

final_plot(root,array,32,period,ax3,fig3)
ax3.set_title("exercise 3- density" )

fig1.savefig("Ballwalk_boundary.png")
fig2.savefig("Test_Code.png")
fig3.savefig("Density_Plot.png")

