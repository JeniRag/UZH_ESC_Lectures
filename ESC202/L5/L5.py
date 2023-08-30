# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 14:01:11 2021

@author: Sujeni

ESC202 - Exercise 5
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
import copy
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

"""classes"""

class particle:
    def __init__(self, r):
        self.r=r
        self.count=0
        self.color="mediumturquoise"
        
        self.rho=0.0
        self.m=1.0 
        self.rho=0.0 #density
        self.a=np.zeros_like(r, dtype='float') #acceleration
        self.v=np.ones_like(r,dtype='float') #velocity
        self.vpred=np.zeros_like(r,dtype='float')
        
        self.edot=np.zeros_like(r, dtype='float')
        self.e= 10.0 #internal Energy
        self.epred=0.0
        self.h=0.0
        self.c= 0.0
        
        self.rb=r*0 #coordinates of neigbours
             
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
    m_sum=mass(queue) 
    rho=m_sum/(2*np.pi*queue.dist2())
    p.rho=rho
    
    return rho

def density_SPH(p,queue): #Monaghan Kernel
    h=np.sqrt(queue.dist2())
    heap=queue.show()
    dens=0
    
    for i in range(len(heap)):
        r_diff=p.r-heap[i][1].r
        r=np.sqrt(r_diff[0]**2+r_diff[1]**2)
        
        dens+=heap[i][1].m*W(r,h)
        
    p.rho=dens
    #print(dW(r,h))
    #print(W(r,h))
    return dens

def W(r,h):
    d=2
    sigma=40/(np.pi*7)
    div=(r/h)
    factor=sigma/h**d
    w=0
    if div>=0 and div<0.5:
        w=factor*(6*div**3-6*div**2+1)
    elif div>=0.5 and div<=1:
        w=factor*(2*(1-(div))**3)
    elif div>1:
        w=0
    return w

"""neighboursearch"""
def Ballwalk(A,c,p,queue,rOffset):
    
    pir=p.r+rOffset
    
    if c.isLeaf==True:
         
        for a in A[c.iLower: c.iUpper+1]: 
            dist2=euc_dist(pir,a.r)  
            if dist2<queue.dist2():
                a_new=copy.deepcopy(a)
                a_new.r-=rOffset
                queue.replace(dist2, a_new)
                
                
    else:
       
        if c.pUpper !=None and c.pLower!=None: 
            #calculate distances once and use results later in order to save some time
            d2_pLower = d2(A,c.pLower,pir)
            d2_pUpper = d2(A,c.pUpper,pir)
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
    d2_list.append(1/heap[i][0])
    return d2_list

"""
def dW(r,h):
    d=2
    sigma=40/(np.pi*7)
    div=(r/h)
    factor=sigma/h**d
    w=0
    if div>=0 and div<0.5:
        w=factor*(18*r**2/(h**3))-12*r/(h**2)
    elif div>=0.5 and div<=1:
        w=factor*1/h*((1-r/h)**2*3) ##########calculate derivative
    elif div>1:
        w=0
    return w
"""

"""exercise 5"""
def dW(r,h): ###ergibt negatives ergebnis
   
  
    d=2
    sigma=40/(np.pi*7)
    div=(r/h)
    factor=6*sigma/h**(d+1)
    w=0
    if  div>=0 and div<0.5:
        w=factor*3*(div)**2-2*div
    elif div>=0.5 and div<=1.0:
        w=factor*-1*(1-div)**2##########calculate derivative
    elif div>1:
        w=0.0
    return w

def Drift1(array, period, dt):
    for p in array:

        p.r+=p.v*dt
        p.r %=period
        #####p.r %= period
        p.vpred= p.v+p.a*dt
        p.epred=p.e +p.edot*dt     ####################edot ??
        
def Drift2(array, period, dt):
    for p in array:
        p.r+=p.v*dt
        p.r %=period

def Kick(array, dt):
    for p in array:
        p.v+=p.a*dt
        
        p.e=p.edot*dt

def gradient(ra,rb, h):
    r=ra-rb
    x=r[0]
    y=r[1]
    r_abs=np.sqrt(x**2+y**2) ############dot product?
    
    return dW(r_abs,h)*r/(r_abs+1/10000) ########r_abs 0???#######################################
    #difference vector r_diff

 
def dvdt(p,queue): #acceleration
    sum=0.0
    gamma=2
    h=np.sqrt(queue.dist2())
    heap=queue.show()
    
    for i in range(len(heap)):
        if heap[i][1] != p:
        
            b=heap[i][1]
            rb=b.r
            sum+=-b.m*(p.c**2/(p.rho*gamma)+b.c**2/(b.rho*gamma))*gradient(p.r,rb,h) #######W and Pi
    return sum

def dedt(p, queue):
    h=np.sqrt(queue.dist2())
    sum=0
    gamma=2
    heap=queue.show()
    for i in range(len(heap)):
        if heap[i][1] != p:
            b=heap[i][1]
    
            sum+=b.m*(p.v-b.v)
            rb=b.r
        
    return p.c**2/(gamma* p.rho)*sum*gradient(p.r,rb,h)

def NN_density(root,array,n,period,density_func,ax,fig):
    #x = np.zeros(len(array))
    #y = np.zeros(len(array))
    queue_list=np.zeros_like((array))
    rho_list= np.zeros(len(array))
    
    for i in range(len(array)):
        #x[i] = array[i].r[0]
        #y[i] = array[i].r[1]
        queue=PRIQ(n)
        p=array[i]
        heap=Ballwalk_periodic(array,root,p,queue,period)
        queue_list[i]=heap
        density_func(p,heap)
        rho_list[i] = p.rho
    return queue_list, rho_list
        
def Calcsound(array):
    gamma=2
    for p in array:
        p.c=np.sqrt(gamma*(gamma-1)*p.epred)

def Calcforce(array, root, dim, tree_list, n, period, density_func, dt, ax, fig):
    treebuild(array,root,dim,tree_list,ax)
    queue_list, rho_list=NN_density(root,array,n,period,density_func,ax,fig)
    Calcsound(array)
    SPH_Force(array, queue_list)
    
    return queue_list, rho_list

def SPH_Force(array,queue_list):
    for i in range(len(array)):
        array[i].a=dvdt(array[i], queue_list[i])
        array[i].edot=dedt(array[i], queue_list[i])

def SPH(array, root, dim, tree_list, n, period, density_func, dt, nstep,ax, fig):
    #Drift1(array, queue, dt)
    queue_list, rho_list=Calcforce(array, root, dim, tree_list, n, period, density_func, dt, ax, fig)
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    X=np.zeros(shape=(nstep, len(array)))
    Y=np.zeros(shape=(nstep, len(array)))
    R=np.zeros(shape=(nstep, len(array)))
    
    for j in range(0, nstep):
        Drift1(array,period, dt/2)
        queue_list, rho_list=Calcforce(array, root, dim, tree_list, n, period, density_func, dt, ax, fig)
        Kick(array,dt)
        Drift2(array,period,dt/2)
    #plot_particle(array,ax)


        for i in range(len(array)):
            x[i] = array[i].r[0]
            y[i] = array[i].r[1]
        X[j]=x
        Y[j]=y
        R[j]=rho_list
    return X, Y, R
        
    
    
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
        
def densities_color(root,array,n,density_func,ax):
    queue=PRIQ(n)
    rho_list= np.zeros(len(array))
    for p in array:
        heap=Ballwalk(array,root,p,queue)
        dens=density_func(p,heap)
        rho_list.append(dens)
    norm=matplotlib.colors.Normalize(dmin=min(rho_list), dmax=max(rho_list))
    
    for p in array:
        p.color=cmap(norm(p.rho))    
                                                                
def plot_circle(p,ax,rmax):
    circle=plt.Circle((p.r[0],p.r[1]), rmax, alpha=0.2)
    ax.add_patch(circle)
    
def final_plot(root,array,n,period,density_func,ax,fig):
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    rho_list= np.zeros(len(array))

    for i in range(len(array)):
        x[i] = array[i].r[0]
        y[i] = array[i].r[1]
        queue=PRIQ(n)
        p=array[i]
        heap=Ballwalk_periodic(array,root,p,queue,period)
        density_func(p,heap)
        rho_list[i] = p.rho

    norm=matplotlib.colors.Normalize(vmin=min(rho_list), vmax=max(rho_list))
   
    cmap = cm.hot
    ax.scatter(x, y, s = 25, cmap = cmap, c = rho_list, alpha = 0.5)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label="density")
                     
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
N=100
period=np.array([1,1])
array=np.zeros(N, dtype=particle)
for i in range(0, N):
    element=np.array([random.random(),random.random()])
    array[i]=particle(element)

    #p1=array[int(random.random()*N)]
    
value=8
nn=32 #neares neighbour
nstep=30
#nEnd=1


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
"""
p=array[100]
p.color="black"

fig1=plt.figure()
ax1=fig1.add_subplot(111)

fig2=plt.figure()
ax2=fig2.add_subplot(111)

treebuild(array,root,dim,tree_list,ax1)
treebuild(array,root,dim,tree_list,ax2)

NN, d2NN=test_code(p,array,nn)
plot_particle(array,ax2)
ax2.set_title("Exercise 3- Testcode")
for i in array:
    i.color="mediumturquoise"

queue=PRIQ(nn)
#print(queue.show())
heap1=Ballwalk_periodic(array,root,p,queue, period)
heap2=heap1.show()

for i in range(len(heap2)):
    heap2[i][1].color="orange"
    
#print(heap1.show())
#d2_list=extract_dist(p,array,heap1)
p.color="black"
plot_particle(array,ax1)
ax1.set_title("Exercise 3- Ballwalk Boundary")
"""

"""density plot """
"""
fig3=plt.figure()
ax3=fig3.add_subplot(111)
treebuild(array,root,dim,tree_list,ax3)

final_plot(root,array,nn,period,density,ax3,fig3)
ax3.set_title("Exercise 3- Density" )
"""
fig4=plt.figure()
ax4=fig4.add_subplot(111)
treebuild(array,root,dim,tree_list,ax4)

final_plot(root,array,nn,period,density_SPH,ax4,fig4)
ax4.set_title("Exercise 4- Monaghan Kernel Density" )

                                                      

"""ex 5 plot"""

cmap = cm.hot

fig5=plt.figure()
ax5=fig5.add_subplot(111)
ax5.set_title("Exercise 5- SPH")
#Calcforce(array, root, dim, tree_list, 32, period, density_SPH, 0.1, ax5, fig5)
X,Y,R=SPH(array, root, dim, tree_list, 32, period, density_SPH, 0.01,nstep, ax5, fig5)

def update(frame,plot,ax):
 
        plot[0].remove()
        plot[0]=ax.scatter(X[frame],Y[frame],s=25, cmap=cmap, c=R[frame], alpha=0.5)
        
        
p1=[ax5.scatter(X[0], Y[0], s=25, cmap=cmap, c=R[0], alpha=0.5)]
norm=matplotlib.colors.Normalize(vmin=np.min(R), vmax=np.max(R))

anim = animation.FuncAnimation(fig5, update, fargs=(p1,ax5),
                               frames=len(R)-1, interval=500, blit=False)   #greater interval is slower

#ax.scatter(x, y, s = 25, cmap = cmap, c = rho_list, alpha = 0.5)
fig5.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax5, label="density")

anim.save('Exercise5.mp4', fps=10, extra_args=['-vcodec', 'libx264'])

plt.show()  
              