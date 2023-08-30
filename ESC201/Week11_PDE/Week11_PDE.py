# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 14:14:12 2020

@author: Sujeni
"""
"""
Week 11-Exercise10
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage



"""Rolling, For-Loop, Convolve"""
def Rolling(rho):
    rho_next=np.roll(rho,-1) #following data at position j
    rho_previous=np.roll(rho,1) #previous data at actual position
    return rho_next, rho_previous

def For_Loop(rho):
    rho_next=np.zeros_like(rho)
    rho_previous=np.zeros_like(rho)
    
    for j in range(len(rho)):
        rho_previous[j]=rho[j-1] #rho[-1] gives already the last element
        if j==len(rho)-1:
            rho_next[j]=rho[0]
        else:
            rho_next[j]=rho[j+1]
    
    return rho_next, rho_previous

def convolve(rho):
   
    rho_next=ndimage.convolve(rho, np.array([1,0,0]), mode='wrap')
    rho_previous=ndimage.convolve(rho, np.array([0,0,1]), mode='wrap')
    
    return rho_next, rho_previous

"""Calculation functions"""
def Lax(rho,method):
    rho_plus_one, rho_minus_one=method(rho)
    rho_new=1/2*(rho_plus_one+rho_minus_one)-0.5*C*(rho_plus_one-rho_minus_one) #Vorzeichen geändert
    rho=rho_new
    return rho

def Upwind_scheme(rho, method):
    rho_next, rho_previous=method(rho)
    rho_new=rho-C*(rho-rho_previous)
    rho=rho_new
    return rho

def Lax_Wendorf(rho, method):
    rho_next, rho_previous=method(rho)
    rho_new=0.5*C*(1+C)*rho_previous+(1-C**2)*rho-0.5*C*(1-C)*rho_next
    rho=rho_new
    return rho
"""Iteration"""
def timestep1(t, rho, func, method):
    for i in range(t):
        rho=func(rho, method)
    return rho

"""Plot function"""
def plot(axs, n, time_step, timestep,func, method, rho):
    index=np.zeros_like(rho)
    for i in range(len(rho)):
        index[i]=i
    for j in range(0,n+1):
        if j % time_step==0:
            rho1=timestep(j,rho,func,method)
            axs.plot(index, rho1)
    axs.set_title(func.__name__)
    
"""Parameters"""
C=0.5
N = 200

"""Initialization"""
rho = np.zeros(N)
for j in range(int(N/2)-10, int(N/2)+10):
  rho[j] = 1
  
"""Plotting"""   
rows, cols = 3,1# SR: overall 4 elements

fig1, axs1 = plt.subplots(rows, cols, sharex=False, sharey=False)
fig1.suptitle("Rolling")
func=[Lax, Upwind_scheme, Lax_Wendorf]
for i in range(0,len(func)):
    plot(axs1[i], 400, 50, timestep1, func[i], Rolling, rho)
    
fig2,axs2 =plt.subplots(rows, cols, sharex=False, sharey=False)
fig2.suptitle("For Loop")
for i in range(0,len(func)):
    plot(axs2[i], 400, 50, timestep1, func[i], For_Loop, rho)

fig3,axs3 =plt.subplots(rows, cols, sharex=False, sharey=False)
fig3.suptitle("Convolve")
for i in range(0,len(func)):
    plot(axs3[i], 400, 50, timestep1, func[i], convolve, rho)
    
fig1.savefig("Rolling.png")
fig2.savefig("For_Loop.png")
fig3.savefig("Convolve.png")

