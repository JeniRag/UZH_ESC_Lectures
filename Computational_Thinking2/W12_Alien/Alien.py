# -*- coding: utf-8 -*-
"""
Created on Thu May 12 15:54:37 2022

@author: Sujeni
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy
"""Simulation of Neutrion radiation"""

A=np.array([-5.0, -3.0]) #Alien position
E=np.array([0.0,0.0]) #Earth position
N=15
N_neut=5000
Nsample=10


#%%
def radiation(Nsteps,alien, earth, detectors):
    b =  alien[0]-earth[0]
    
    for i in range(Nsteps):
        alpha=np.random.uniform(-np.pi/2, np.pi/2)
        r = b/np.cos(alpha)
        l=r*np.sin(alpha)
        a = alien[1]-earth[1]
        
        detector_y= int(alien[1] +  l )
        x=np.tan(alpha)*b+a
        # x =int(x)
        # if 0<= x <= len(detectors)-1:
        #     detectors[x] +=1
        
        if 0 <= detector_y <= len(detectors)-1:
            detectors[ detector_y] += 1
    
    return detectors

def Prior(detectors):
    return detectors/np.sum(detectors)


A=np.array([5, 50.0]) #Alien position #x,y
E=np.array([0.0,0]) #Earth position
x_label=np.arange(-50,50)
detectors=np.zeros(100)
detectors=radiation(70000,A,E, detectors)
plt.plot(x_label,detectors)
plt.xlabel("detectors")
plt.ylabel("count")


#%%
prior=Prior(detectors)  
plt.plot(prior)  
  
  
#%%
def prior(x): #zähler
    
    return

def evidence(x,data):#x float, data 1Darray (detectors)
    return data[int(x)]/np.sum(data)

def integral(x, a,b):
    return ((-b*np.log(np.cos(np.pi/2))+a) -(-b*np.log(np.cos(-np.pi/2))+a))*x

def proposal_distribution(x):
    return np.randon.normal(x, np.full(2,1))

def manual_log_like_normal(mu,s,data):
    return np.sum( (-np.log(mu * np.sqrt(2*np.pi)) - ((data-mu)**2) / (2*s**2)))


def acceptance(proposed, current):
    A=np.random.normal(proposed,np.full(2,1))/np.random.normal(proposed,np.full(2,1))
    p=min(1,A  )
    u=np.random.random()
    
    if p>1:
        accepted.append(current)
    
    elif u<=p:
        accepted.append(proposed)
        current=proposed
    else:
        accepted.append(current)

def MCHC(N, prior, likelihood):
    
    d=detectors.max()
    data_average=np.abs(detectors-d).argmin()-len(detectors)
    
    print(data_average)
    
    accepted=[]
    rejected=[]
    d=2
    p0=np.random.uniform(-100,100,d) #proposed a,b
    current=p0
    sigma=1
    
    for i in range(N):
        proposed = np.random.normal(current,sigma)
        # (likelihood(data_average, proposed)*prior(data_average, proposed))/( likelihood(data_average, current) * prior(data_average, proposed))
        A=np.sum(np.random.normal(proposed,sigma)/np.random.normal(current,sigma)) 
        
        p=min(1,A )
        u=np.random.uniform(0,1)
        
        if u<=p:
            #accepted.append(np.array(proposed, dtype=int))
            accepted.append(proposed)
            current=proposed
        else:
            #accepted.append(np.array(current, dtype=int))
            accepted.append(current)
            
    return np.array(accepted)
            

a=MCHC(1000, evidence,manual_log_like_normal)
# plt.plot(a,"o")

#acc,count= np.unique(a, return_counts=True)
plt.figure()
plt.hist2d(a[:,0], a[:,1])

a=A[0]-E[0]
b=A[1]-E[1]

print(a,", ",b)
#%%%
from sympy import Symbol, sin, cos, tan, sec, csc, cot
from sympy.integrals.trigonometry import trigintegrate
from sympy.abc import x
  
# Using trigintegrate() method
gfg = trigintegrate(tan(x), x)
  
print(gfg)
        

#%%
random.seed(10)
P_thetha=1/np.pi

def Likelihood(xi,A):
    return
    
def P_xi(i, detectors, count_list):
    return count_list[i]/( np.sum(count_list))


def Q(y,x):

    return 1/np.sqrt(2*np.pi)*np.exp( -0.5*(y[0]-x[0])**2 + (-0.5)*(y[1]-x[1])**2 )

def pi(theta):
    sx=1
    sy=1
    my=0
    mx=0
    return 1/np.sqrt(2*np.pi*sx*sy)* np.exp( (theta[0]-mx)**2/sx + (theta[1]-my)**2/sy)

def pi_theta_A(): # theta| location of Aliens
    return 1/(2*np.pi)

def Like(j, A):
    
    return pi_theta_A()/pi_xi(j)

def MetropolisHasting(target, dtectors):
  

    current=np.zeros(2)
    current[0]=np.random.uniform(-N,E[0])
    current[1]=np.random.uniform(-N,N)
    accepted=[]
    
    for j in range(len(detectors)):
        for i in range(10):
            proposed=current+np.random.normal(0,1, 2)
            A=target(proposed)/target(current) *Q(current, proposed)/Q(proposed, current)
            p=min(1,A  )
            u=np.random.random()
            if p>1:
                accepted.append(current)
            
            elif u<=p:
                accepted.append(proposed)
                current=proposed
            else:
                accepted.append(current)
        
    return np.array(accepted)
    
thetas=np.random.uniform(-N, N , size=(Nsample,2) )
l=[]

    
e=MetropolisHasting(pi,detectors)

plt.figure()
plt.plot(e[:,0], e[:,1], "o")
plt.figure()
plt.hist(e[:,0], bins=100)
plt.hist(e[:,1], bins=100)

plt.hist2d(e[:,0], e[:,1])

