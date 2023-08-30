# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:18:15 2020

@author: Sujeni
"""
"""
ESC201
Week3_Exercise2_Feigenbaum
"""
import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer


n_x0=300 #number of starting values
n_a=800 #number of a values
n_xn=1000 #number of iteration

start=timer()

def LogisticEquation(xn,a,n):
    for i in range(0,n):
        xnew=a*xn*(1-xn)
        xn=xnew
    return xn

a_List=np.linspace(0,4,n_x0)
x0_List=np.linspace(0,1,n_a)

A,X=np.meshgrid(a_List,x0_List)
y=LogisticEquation(X,A,n_xn)

fig, ax=plt.subplots()
ax.set_xlim([0,4])   
ax.set_ylim([0,1])
ax.set_xlabel("a")
ax.set_ylabel("x (" +str(n_xn)+" iterations)")
ax.set_title("Exercise2: Feigenbaum")


for i in range (len(A)):
    
    ax.scatter(A[i],y[i],s=0.05, c='black')

#plt.show()
plt.savefig("./Feigenbaum.png", dpi=300)

end=timer()
time=round(end-start,2)
print("Duration:" +str(time)+ "s")