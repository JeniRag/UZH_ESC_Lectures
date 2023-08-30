# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 14:27:57 2020

@author: Sujeni
"""
"""
Week 9-Exercise 8
"""
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import math


def CupperplateVolt(U,R,ax):
    #Cupper plate Volt
    L,J=U.shape
    U[int(round(L/2)),int(np.ceil(J*0.25)):int(round(J*0.75)):1]=1000
    R[int(round(L/2)),int(np.ceil(J*0.25)):int(round(J*0.75)):1]=0
    y=L/2
    x1=J/4
    x2=J*3/4
    ax.plot([x1,x2], [y,y])
    #boundary conditions
    for l in range(0,L):
        for j in range (0,J):
            if j==0 or j==J-1:
                U[l,j]=0
                R[l,j]=0
            if l==0 or l==L-1:
                U[l,j]=0
                R[l,j]=0
                
def Matrixes(U,omega,ax):
    L=Grid_y
    J=Grid_x
    M=np.zeros_like(U)
    P=np.zeros_like(U)
    Unew=np.zeros_like(U)
    
    
    return ( M, P, R, Unew)

#Weight
W=np.array([[0,1,0], [1,-4,1],[ 0,1,0]])

#create checked condition matrix
def condition(U):
    C=np.zeros_like(U, dtype=bool)
    C[::2, ::2]=True
    C[1::2, 1::2]=True
    return C
   
def one_step(U,Unew, W,R,M,C,P):
    ndimage.convolve(U,W, output=P, mode = "constant", cval=0)
    np.multiply(R,P,out=M)
    Unew[C]=U[C]+M[C] #update black points
    
    ndimage.convolve(Unew,W, output=P, mode = "constant", cval=0)
    np.multiply(R,P,out=M)
    Unew[~C]=U[~C]+M[~C] #update red points
    
    return Unew

def SOR(U,J,L,W,omega,threshold,ax): 
    output=[]
    M, P,R, Unew= Matrixes(U,omega,ax)
    C=condition(U)
    
    delt=100 #just an initial value
    N=0 
    average=[np.average(U)]
    
    while (delt)>threshold: #was while np.linalg.norm(delt)>threshold: before
        Unew=one_step(U,Unew,W,R,M,C,P)
        average.append(np.average(Unew))
        #delt=average[len(average)-1]-average[len(average)-2]
        delt = np.max(np.abs(M))
        U=Unew
        N+=1
        
    output.append(U)
    output.append(N)
    return output


"""exercise 8"""

e=-1.6*10**(-19) #charge of electron in C
Me=9.11*10**(-31) #Mass of electron in kg

def index(x,y):
    x1=np.array(x/dx, dtype=int) #dx=0.01
    y1=np.array(y/dy, dtype=int)

    return (x1,y1)

def acc(x,y):
    global U
    
    j, l=index(x,y)
    J, L=U.shape
    
    j=np.maximum(j, np.zeros_like(j))
    j=np.minimum(j, (J-2)*np.ones_like(j))
    l=np.maximum(l, np.zeros_like(l))
    l=np.minimum(l, (L-2)*np.ones_like(l))
   
    
    t=(x-j*dx)/dx
    u=(y-l*dy)/dy
    
    phi1=U[l,j] #Unew=the calculated potential from plot1
    phi2=U[l, j+1]
    phi3=U[l+1, j]
    phi4=U[l+1, j+1]
    
    ax=(e/Me)*(-1/dx)*((1-u)*(phi2-phi1)+u*(phi4-phi3))
    ay=(e/Me)*(-1/dy)*((1-t)*(phi3-phi1)+t*(phi4-phi2))

    return (ax,ay)


def LeapFrog(t,acc,h, vector): #input vector: (x,y,vx,vy)
    x,y=vector[0], vector[1]
    vx,vy=vector[2], vector[3]
    ax,ay=acc(x,y)
    
    x_half= x+h/2*vx
    vx_new=vx+h*ax
    x_new=x_half+h/2*vx_new
    
    y_half= y+h/2*vy
    vy_new=vy+h*ay
    y_new=y_half+h/2*vy_new
    
    return (x_new,y_new,vx_new,vy_new)

def odeSolver(t,vector0, dfdt, h, nSteps, solverStepFunc):
   
    x_list=[vector0[0]]
    y_list=[vector0[1]]
    #vx_list=[vector0[2]]
    #vy_list=[vector0[3]]
    vector=vector0 #initial vector (x,y,vx,vy)
    t_list=[t]
    
    for n in range(nSteps):
        
        x_new, y_new,  vx_new, vy_new = solverStepFunc(t,acc, h, vector) 
        
        x_list.append(x_new) 
        y_list.append(y_new)
            
            #vx_list.append(vx_new) 
            #vy_list.append(vy_new)
            
        vector=np.array([x_new,y_new, vx_new, vy_new])
       
        
        t+=h
        t_list.append(t)
    return (x_list,y_list, t_list) #returns x,y

def generateElectron(n):
    
    y=np.linspace(0.6*Box_width, 0.9*Box_heigth,n)
    x=np.zeros_like(y)
    phi=np.linspace(-math.pi/2, math.pi/2, n)
    vx = 1e6*np.ones_like(x)
    vx=np.ones_like(y)
    vy = np.ones_like(y)
    vx = 1e6 * np.cos(phi)*vx
    vy = 1e6 * np.sin(phi)*vy
    return np.array([x, y, vx, vy])


def placePlate(p1, p2, R, U, potential, ax):
    x1, y1 = p1
    x2, y2 = p2
    if x2 < x1:
        t = x2
        x2 = x1
        x1 = t
        t = y2
        y2 = y1
        y1 = t
    a = (y2 - y1) / (x2 - x1)
    j1 = int(x1 / dx)
    j2 = int(x2 / dx)
    l1 = int(y1 / dy)
    l2 = int(y2 / dy)

    n = max(j2-j1+1, l2-l1+1)
    for i in range(n+1):
        x = x1 + i*(x2 - x1)/n
        y = y1 + a*(x - x1)
        j = int(x / dx)
        l = int(y / dy)
        R[l, j] = 0
        U[l, j] = potential
    ax.plot([x1, x2], [y1, y2], color="black")
    
def Box(Box_width,Box_heigth, axs):
    x1,y1=Box_width, Box_heigth/Box_heigth*1/4
    x2, y2= Box_width, Box_heigth/Box_heigth*3/4
    axs.plot([x1,y1],[x2,y2], color='green')
        
        
"""Parameters"""

N=200
omega=2/(1+np.pi/N)
Grid_x, Grid_y=N,N
threshold=0.001
rows, cols = 2, 2 # SR: overall 4 elements
fig, axs = plt.subplots(rows, cols, sharex=False, sharey=False)
plot1=axs[0][0]
plot2=axs[0][1]
plot3=axs[1][0]
plot4=axs[1][1]
"""Box Settings"""
Box_width, Box_heigth= 0.01, 0.01 #0,01
dx=Box_width/(Grid_x-1)
dy=Box_heigth/(Grid_y-1)

"""plot Potential"""
U=np.zeros(shape=(Grid_x,Grid_y))
R=np.ones_like(U)*omega/4
w=2/(1+np.pi/N)

#CupperplateVolt(U,R,axs[0])
placePlate([0.002,0.002],[0.004, 0.001],R,U,100,plot2)
placePlate([0,0.004],[0.001, 0.005],R,U,-1,plot2)
placePlate([0.006,0],[0.009, 0.00],R,U,-5,plot2)
placePlate([0.0098,0.004],[0.0099, 0.0041],R,U,230,plot2)
placePlate([0.009, 0.001], [0.0099, 0.001], R, U, 40, plot2)
#placePlate([0.005,0.001],[0.008, 0.001],R,U,-10,axs[1])
#placePlate([0.009,0.004],[0.0091, 0.006], R, U, 100, axs[1])
start=timer()


U,N=SOR(U,Grid_x,Grid_y,W,omega,threshold,plot1)


end=timer()
time=round(end-start,2)
plot1.contour(U)

levels=(np.arange(0, 1000, 100))
plot1.contour(U, levels=levels)
plot1.set_title("Contourplot")
plot1.text(5, 1, "iterations:"+str(N)+", time: "+str(time) )

"""opening and detectors"""
#plot2.plot([0.01,0.001],[0.01, 0.004], color='green', linewidth=2)
#Box(Box_width, Box_heigth, axs[1])

"""electrons"""
n_electron=40
vector0=generateElectron(n_electron)

tEnd=10*10**(-9)
h=10*10**(-12)
N2=round(tEnd/h)

x,y, t= odeSolver(0,vector0, acc, h, N2, LeapFrog)
"""histogram"""
histogramm_N=np.zeros_like(t,dtype='int')

t_detector=np.zeros(shape=( n_electron))
N_total=0

x_array=np.array(x)
y_array=np.array(y)

for i in range(0, len(t)-1):
    N=0
    for j in range(0, n_electron-1): #electrons which hit the box shouldn't be counted when they hit the box. So I keep them out of the box for the detector calculation.
        if x_array[i][j]<0:
            x_array[i:n_electron-1,j]=-0.001
        if x_array[i][j]>Box_width:
            x_array[i:n_electron-1,j]=Box_width
        if y_array[i][j]<0:
            y_array[i:n_electron,j]=-0.001
        if y_array[i][j]>Box_heigth:
            y_array[i:n_electron-1, j]=Box_heigth
        
       
for i in range(0, len(t)-1):
    N=0
    for j in range(0, n_electron-1):
        if x_array[i][j]>0.00999 and x_array[i][j]<0.01 and y_array[i][j]>0.001 and y_array[i][j]<0.004: # when it hits the detector
            N+=1
            t_detector[j]=t[i]
            
    histogramm_N[i]=N
N_total=sum(histogramm_N)
rate=round(N_total/n_electron*100  ,1)
    
plot3.vlines(t, ymin=np.zeros_like(histogramm_N), ymax=histogramm_N)
plot3.set_xlabel("time [s]")
plot3.set_ylabel("electrons hit detector")
plot3.text( 0,0, "N total:"+str(N_total) +", hitting rate:"+str(rate)+"%")



plot2.plot(x,y)
plot2.set_xlim(0,Box_width)
plot2.set_ylim(0,Box_heigth)
plt.savefig("./Electron_Movement.png", dpi=300)

