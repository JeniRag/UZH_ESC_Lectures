# -*- coding: utf-8 -*-
"""
Created on Sat Nov  7 10:16:46 2020

@author: Sujeni
"""
"""
Week 8-Exercise 7
"""
import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from timeit import default_timer as timer

omega=np.array([2/(1+np.pi/101),1,1.5,2])

L=101#rows
J=101#columns
threshold=0.01

def Matrixes(L,J, omega):
    U=np.zeros(shape=(L,J))
    R=np.ones_like(U)*omega/4
    M=np.zeros_like(U)
    P=np.zeros_like(U)
    Unew=np.zeros_like(U)
    
    #Cupper plate Volt
    U[int(round(L/2)),int(np.ceil(J*0.25)):int(round(J*0.75)):1]=1000
    R[int(round(L/2)),int(np.ceil(J*0.25)):int(round(J*0.75)):1]=0
    
    #boundary conditions
    for l in range(0,L):
        for j in range (0,J):
            if j==0 or j==J-1:
                U[l,j]=1
                R[l,j]=0
            if l==0 or l==L-1:
                U[l,j]=1
                R[l,j]=0
    return np.array([ M, P, R, U, Unew])

#Weight
W=np.array([[0,1,0], [1,-4,1],[ 0,1,0]])

#create checked condition matrix
def condition(L,J):
    
    C=np.zeros(shape=(L,J), dtype=bool)
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

def SOR(J,L,W,omega,threshold): 
    output=[]
    M, P,R = Matrixes(L,J,omega)[0], Matrixes(L,J,omega)[1], Matrixes(L,J,omega)[2]
    U, Unew= Matrixes(L,J,omega)[3],Matrixes(L,J,omega)[4]
    C=condition(J,L)
    
    delt=100 #just an initial value
    N=0 
    average=[np.average(U)]
    
    while np.linalg.norm(delt)>threshold:
        Unew=one_step(U,Unew,W,R,M,C,P)
        average.append(np.average(Unew))
        delt=average[len(average)-1]-average[len(average)-2]
        U=Unew
        N+=1
        
    output.append(U)
    output.append(N)
    return output

rows, cols = 1, 4 # SR: overall 4 elements
fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True)
axs.flatten()
for i, w in enumerate(omega): #plot for different omega values
    start=timer()
    plot=axs[i]
    Unew=SOR(J,L,W,w,threshold)[0]
    N=SOR(J,L,W,w,threshold)[1]
    end=timer()
    time=round(end-start,2)
    plot.contour(Unew)
    axs[i].set_title("omega="+str(w))
    axs[i].text(5, 1, "iterations:"+str(N)+", time: "+str(time) )
    
    
plt.savefig("./Exercise7.png")
