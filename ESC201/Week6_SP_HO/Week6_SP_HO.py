# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 20:51:04 2020

@author: Sujeni
"""
"""
ESC-201-Exercise 5
"""
import numpy as np
import matplotlib.pyplot as plt
import math


e= 6 #epsiolon
n=6

H0=np.linspace(1,6,n)

q0=-0.3
q01=-0.2 #angle


nSteps=1000
h=0.1

def SimplePendulum(H,q):
    p=math.sqrt((H+e*np.cos(q))*2) 
    return p

def HarmonicOscillator(H,q):
    p=math.sqrt(2*H-q**2) 
    return p
    
def dSP(q): #derivative of q in SimplePendulum
    dUdq=e*np.sin(q)
    return dUdq

def dHO(q): #derivative of q in HarmonicOscillator
    dUdq=q
    return dUdq

def LeapFrog(p,q,h,df): 
    qhst=q+1/2*h*p #qnhst=q halfstep
    pnew= p+h*(-df(qhst))
    qnew=qhst+1/2*h*pnew
    
    return np.array([pnew,qnew])

def odeSolver(H0,q0, func, df, h, nSteps, solverStepFunc):
   
    qn=q0
    pn=func(H0,q0)
    plist = [pn]
    qlist= [qn]
   
    for n in range(nSteps):
        
        pnew = solverStepFunc(pn, qn, h, df)[0]
        qnew = solverStepFunc(pn, qn, h, df)[1]
                
        plist.append(pnew)
        qlist.append(qnew)
        
        pn=pnew
        qn=qnew
        
    return (np.array(plist), np.array(qlist))

fig1, ax1= plt.subplots()

ax1.set_title("Simple Pendulum")
ax1.set_xlabel("q")
ax1.set_ylabel("p")

fig2, ax2= plt.subplots()
ax2.set_title("Harmonic Oscillator")
ax2.set_xlabel("q")
ax2.set_ylabel("p")

for H in H0:
    p1, q1 = odeSolver(H,q01, SimplePendulum, dSP, h, nSteps, LeapFrog)
    p2, q2 = odeSolver(H,q0, HarmonicOscillator, dHO, h, nSteps, LeapFrog)
    
    ax1.scatter(q1,p1, s=0.1, label=H)
   
    ax2.scatter(q2, p2, s= 3, label=H)
    
ax1.legend()
ax2.legend()

fig1.savefig("./Simple_Pendulum.png")
fig2.savefig("./Harmonic_Oscillator.png")


    