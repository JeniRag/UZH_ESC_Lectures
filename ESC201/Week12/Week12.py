# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:30:57 2020

@author: Sujeni
"""
"""
Week 12-Exercise11
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import ndimage

"""Parameters"""
dx=1
dy=1
dt=1

sigma_x=3
sigma_y=3

x0=sigma_x*2
y0=sigma_y*2

N=50
A=10

Ca=0.4
Cb=0.4

a=Ca*dx/dt
b=Cb*dy/dt

N_time=200

"""Initialization"""
rho = np.zeros(shape=(N,N))
for x in range(N):
    for y in range(N):
        rho[x][y]= A*np.exp(-((x-x0)**2/(2*sigma_x**2)+(y-y0)**2/(2*sigma_y**2)))
        
"""Function definitions"""
def CTU(rho): #Cornor transport upwind method
    rho_left=np.roll(rho, 1, axis=0) #rho(j-1,l)
    
    rho_star=(1-Ca)*rho+Ca*rho_left
    rho_new=(1-Cb)*rho_star+Cb*np.roll(rho_star, 1, axis=1)
    return rho_new

def CIR(rho):
    rho_left=np.roll(rho,1,axis=0)
    rho_down=np.roll(rho, 1, axis=1)
    
    rho_new=rho-Ca*(rho-rho_left)-Cb*(rho-rho_down)

    return rho_new

def plot(rho, func, axs):
    
    for i in range(N_time):
        rho_new=func(rho)
        if i%50==0:
            axs.contour(rho_new)
        rho=rho_new
    axs.set_title(func.__name__)
 

"""plotting"""
rows, cols= 1,2
fig1, axs1 = plt.subplots(rows, cols, sharex=True, sharey=True)

plot1=axs1[0]
plot2=axs1[1]

plot1.set_xlim(0,N)
plot1.set_ylim(0,N)

rho_CTU=plot(rho,CTU, plot1)
rho_CIR=plot(rho,CIR, plot2)