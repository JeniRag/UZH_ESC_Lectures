# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 15:17:06 2020

@author: Sujeni
"""
"""
Week2-Exercise 1- Orbit
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math


fig, ax = plt.subplots()
xdata, ydata = [], []
ln, = plt.plot([], [], 'ro', animated=True)

#Parameters
a=1 # in AU 
e =0.5 
#functions
def f(E,M): #Kepler equation
    return E - e*np.sin(E) - M

T=pow(a,3/2)

def fder(E): #first derivative of Keplar equation
    return 1 - e*np.cos(E)

def newton(f, fder, M, Estart):
    """Newtons algorithm"""
   
    E=Estart
    
    correction=100
    
    while correction>0.1:
        correction=-f(E,M)/fder(E)
        Enew=E+correction
        E=Enew
       
    return E
      
def init():
    plt.cla()
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.scatter(0,0, color='y')
    return ln,

def update(frame): 
    M = 2*np.pi*(frame/T)
    E = newton(f, fder, M, M)
    x =a*np.cos(E)-a*e
    y =a*math.sqrt(1-e*e)*np.sin(E)
    xdata.append(x)
    ydata.append(y)
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, frames=np.linspace(0, 100,50),
                    init_func=init, blit=True, interval=50)
plt.show()

