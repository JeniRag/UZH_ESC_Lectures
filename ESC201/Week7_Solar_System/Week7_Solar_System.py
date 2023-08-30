# -*- coding: utf-8 -*-
"""
Created on Sat Oct 31 10:33:20 2020

@author: Sujeni
"""
import pandas as pd
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


def readPlanets(filename, N=-1):
    # Loading text files with Pandas
    df = pd.read_csv(filename, sep=',', header=None,
                     names=['name', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'])

    # Data is now in a Pandas dataframe
    # print(df)

    name = np.array(df.loc[:, 'name'])
    m = np.array(df.loc[:, 'm'])
    r = np.array(df.loc[:, 'x':'z'])
    v = np.array(df.loc[:, 'vx':'vz'])

    if N > 0:
        name = name[0:N-1]
        m = m[0:N-1]
        r = r[0:N-1]
        v = v[0:N-1]

    return (name, r, v, m)
name, r, v, m = readPlanets("SolSystData.dat")


h = 4

Msolar=1.99*10**30

k1=6.6*10**(-11)*Msolar#k^2


fin=100 #end time


r12=np.zeros((9,3), dtype=float)
Fx=np.zeros((9,9), dtype=float)
Fy=np.zeros((9,9), dtype=float)
Fz=np.zeros((9,9), dtype=float)
#F=np.zeros((9,9),dtype=float)
ax=np.zeros((9), dtype=float)
ay=np.zeros((9), dtype=float)
az=np.zeros((9), dtype=float)
a=np.zeros((9,3), dtype=float)


vnew=np.empty((9,3),dtype=float)
rnew=np.empty((9,3), dtype=float)


def acc(r12):
    for i in range(0,9):
        for j in range(i+1,9):
            
            dr=r12[j]-r12[i]
            drx, dry, drz =dr[0],dr[1], dr[2]
            
            
            drn=dr[0]**2+dr[1]**2+dr[2]**2
            drn=np.sqrt(drn)
            
            Fx[i][j]=k1*m[i]*m[j]*(drx) /(drn)**3
            Fy[i][j]=k1*m[i]*m[j]*(dry) /(drn)**3
            Fz[i][j]=k1*m[i]*m[j]*(drz) /(drn)**3
            
            
            ax[i]+=Fx[i][j]/m[i] #Fij
            ax[j]-=Fx[i][j]/m[j] #Fji=-Fij
           
            ay[i]+=Fy[i][j]/m[i] #Fij
            ay[j]-=Fy[i][j]/m[j] #Fji=-Fij
            
            az[i]+=Fz[i][j]/m[i] #Fij
            az[j]-=Fz[i][j]/m[j] #Fji=-Fij
           
            #F[j]=k1*m[i]*m[j]*(r12[i]-r12[j])/((drn)**3)
            #a[j]-=F[j]/m[j]
        a[i]=np.array([ax[i],ay[i],az[i]])
    return a

rlist=[r]
# 2. for time in 0 to end time:
for N in range(0,fin,h): 
    
    # 3. First drift to give r1/2
    r12=r+h/2*v #r1/2
    #4. Calculate forces and accels at postion r1/2
    a=acc(r12)
                
#    5. Kick: v1 = v0 + ... using accels from step 4.
    vnew=v+h*acc(r12)
#    6. Second drift to give r1
    rnew=r12+h/2*vnew
    #print("rnew", N, rnew)
    
#    7. Store r1 for graphics
    rlist.append(rnew)

    r=rnew
    v=vnew
    
       
#       8. Plot all planets
innerSS=np.array([0,1,2,3,4])
#innerSS=np.array([0,1])
outerSS=np.array([0,5,6,7,8])
rows, cols = 2, 2 # SR: overall 4 elements
fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True) 
	
fig.tight_layout() 
plot1 = axs[0][0]
plot2 = axs[1][0]
plot3 = axs[0][1]
plot4 = axs[1][1]

def draw(title, System,plotxy,plotxz):
    for i in System:
        rx=[]
        ry=[]
        rz=[] 
        for j in range(0, len(rlist)):
        
            rx.append(rlist[j][i][0])
            ry.append(rlist[j][i][1])
            rz.append(rlist[j][i][2])
            
        plotxy.plot(rx,ry, label=name[i])
        plotxy.set_title(str(title)+ " XY")
       
        plotxy.legend()
        
        plotxz.plot(rx,rz, label=name[i])
        plotxz.set_title(str(title)+ " XZ")
        
        plotxz.legend()
       
draw("Inner Solar System", innerSS,plot1,plot2)
draw("Outer Solar System", outerSS,plot3,plot4)
