# -*- coding: utf-8 -*-
"""
Created on Wed May 11 18:14:33 2022

@author: Sujeni
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import math

N=100
case1= [50,50]
case2=[65,35]
case3=[80,20]
cases=[case1, case2, case3]
dx=0.001
x=np.linspace(0,1,int(1/dx))


from scipy.integrate import quad


def factorial(n,x):
    z=math.factorial(n)
    n=math.factorial(x)*math.factorial(n-x)
    return z/n

class Model():
    def __init__(self,alpha, beta):
        self.alpha=alpha
        self.beta=beta
        
        
    def beta_dist(self): #P(tehta)
        return scipy.stats.beta.pdf(x, self.alpha, self.beta)

    def prob(self, case):
        #factorial(N, case[0])
        L= x**case[0]*(1-x)**case[1] #Likelihood
        PDM=L*self.beta_dist() #theta^N*(1-theta)^T*B(theta,alpha,beta)
        return PDM, np.trapz(PDM,dx=dx)
    
        
M1=Model(1,1)
M2=Model(30,30)


M_list=[M1,M2]
p1=[]
p2=[]
p3=[]

Lc1=[]
Lc2=[]
Lc3=[]

def plot_Lc(Lc_list):
    figure = plt.figure()
    for i,L in enumerate(Lc_list):
        plt.subplot(1,2,1)
        plt.title("Model 1")
        plt.plot(L[0], label=("case " + str(i+1)))
        
        plt.legend()
        plt.subplot(1,2,2)
        plt.title("Model 2")
        plt.plot(L[1], label=("case "+ str(i+1)))
    
    plt.legend()
    
    
for M in M_list:
        
    L1,I1=M.prob(case1)
    p1.append(I1)
    Lc1.append(L1)
    
    L2,I2=M.prob(case2)
    p2.append(I2)
    Lc2.append(L2)
    
    L3,I3=M.prob(case3)
    p3.append(I3)
    Lc3.append(L3)


def BF(p):
    return p[0]/p[1]

ratios=[BF(p1), BF(p2), BF(p3)]
for i,r in enumerate(ratios):
    print("case ",i+1, "H,T=", cases[i][0], "," , cases[i][1], " -->BF = ", np.round(r,3))

plt.figure()
plt.title("beta distribution")
plt.plot(x, M1.beta_dist(), label="M1 ")
plt.plot(x, M2.beta_dist(), label="M2 ")
plt.legend()


# plt.figure()
# plt.title("Likelihood")
# plt.plot(x, M1.t1, label="M1")
# plt.plot(x, M2.t1, label="M2")
# plt.legend()

plot_Lc([Lc1, Lc2])