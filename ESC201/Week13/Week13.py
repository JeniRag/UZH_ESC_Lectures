# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:30:28 2020

@author: Sujeni
"""
""" Week 13 - Exercise 12 """
import numpy as np
import matplotlib.pyplot as plt

"""Parameters"""
N=100
e1=10e-5
Dmax=0.3
dx=1
dt=1/2
N_time=1000
gamma=3
e2=1
"""Initial Conditions"""
def Schock_Tube(N):
    rho=np.ones(N)
    U=np.zeros(shape=(N,3))
   
    u=np.zeros_like(rho)
    
    for i in range (int(N/2), N):
        rho[i]=4
        
    E=np.zeros_like(rho)
    
    E=1/2*rho*u**2+rho*e1
    U[:,0]=rho
    U[:,1]=u
    U[:,2]=E
    
    return U

def Blast_Wave(N):
    rho=np.ones(N)
    U=np.zeros(shape=(N,3))
 
    u=np.zeros_like(rho)
    
    E=np.zeros_like(rho)
      
    P=e1*rho*2
    P[int(N/2)]=1*rho[int(N/2)]*2
  
    E=1/2*rho*u**2+1/2*P
    U[:,0]=rho
    U[:,1]=rho*u
    U[:,2]=E

    return U

"""Methods"""
def Flux(U):   
    P=U[:,2]*2-U[:,1]*U[:,1]/U[:,0]
    
    F=np.zeros_like(U)
    
    F[:,0]=U[:,1]
    F[:,1]=U[:,1]*U[:,1]/U[:,0]+P
    F[:,2]=U[:,1]/U[:,0]*(U[:,2]+P) #u*(E+P)
    
    return F

def Flux_news(Flux,U):
    U_left=np.roll(U, 1, axis=1)
    F_new=1/2*(Flux(U)+ Flux(U_left))-1/2 *Dmax(U-U_left)
    return F_new

def LAX(U): #Method A
    U_left=np.roll(U,1,axis=0)
    U_right=np.roll(U,-1,axis=0)
   
    F=Flux(U)
    F_left=np.roll(F,1,axis=0)
    F_right=np.roll(F,-1,axis=0)
    
    U_new=1/2*(U_right+U_left)-dt/(2*dx)*(F_right-F_left)
    return U_new
    
def Method_B(U):
    P=U[:,2]*2-U[:,1]*U[:,1]/U[:,0]
   
    Cs=((gamma*np.abs(P)/U[:,0]))**0.5 
    Cs_left=np.roll(Cs,1)
    Cs_right=np.roll(Cs,-1)
    
    U_left=np.roll(U,1, axis=0)
    U_right=np.roll(U,-1, axis=0)
    
    Dmax_left=np.max(np.maximum(np.abs(U_left[:,1]/U_left[:,0])+Cs_left,np.abs(U[:,1]/U[:,0])+Cs)) #|ui-1|+Csi-1, |ui|+Csi
    
    Dmax_right=np.max(np.maximum(np.abs(U_right[:,1]/U_right[:,0])+Cs_right, np.abs(U[:,1]/U[:,0])+Cs))

    F_half_left=1/2*(Flux(U)+Flux(U_left))-1/2*Dmax_left*(U-U_left)
    F_half_right=1/2*(Flux(U)+Flux(U_right))-1/2*Dmax_right*(U_right-U)
    
    U_new=U-dt/dx*(F_half_right-F_half_left)
    
    return U_new

def Method_C(U):
    P=U[:,2]*2-U[:,1]*U[:,1]/U[:,0]
    Cs=(gamma*P/U[:,0])**0.5
    Cs_left=np.roll(Cs,1)
    Cs_right=np.roll(Cs,-1)

    U_left=np.roll(U,1, axis=0)
    U_right=np.roll(U,-1, axis=0)
    
    Dmax_left=np.max(np.maximum(np.abs(U_left[:,1]/U_left[:,0])+Cs_left,np.abs(U[:,1]/U[:,0])+Cs)) #|ui-1|+Csi-1, |ui|+Csi
    
    Dmax_right=np.max(np.maximum(np.abs(U_right[:,1]/U_right[:,0])+Cs_right, np.abs(U[:,1]/U[:,0])+Cs))
    
    F_half_left=1/2*(Flux(U)+Flux(U_left))-1/2*Dmax_left*(U-U_left)
    F_half_right=1/2*(Flux(U)+Flux(U_right))-1/2*Dmax_right*(U_right-U)
    
    ##intermediate
    U_int=U-dt/(2*dx)*(F_half_right-F_half_left)
    U_int_left=np.roll(U,1, axis=0)
    U_int_right=np.roll(U,-1, axis=0)
    
    P_int=U_int[:,2]*2-U_int[:,1]*U_int[:,1]/U_int[:,0]
    Cs_int=(gamma*P_int/U_int[:,0])**0.5
    Cs_int_left=np.roll(Cs,1)
    Cs_int_right=np.roll(Cs,-1)
    
    Dmax_int_left=np.max(np.maximum(np.abs(U_int_left[:,1]/U_left[:,0])+Cs_int_left,np.abs(U_int[:,1]/U_int[:,0])+Cs_int)) #|ui-1|+Csi-1, |ui|+Csi
    
    Dmax_int_right=np.max(np.maximum(np.abs(U_int_right[:,1]/U_int_right[:,0])+Cs_int_right, np.abs(U_int[:,1]/U_int[:,0])+Cs))
    
    F_half_left_int=1/2*(Flux(U_int)+Flux(U_int_left))-1/2*Dmax_int_left*(U_int-U_int_left)
    F_half_right_int=1/2*(Flux(U_int)+Flux(U_int_right))-1/2*Dmax_int_right*(U_int_right-U_int)
    
    U_new=U-dt/dx*(F_half_right_int-F_half_left_int)

    return U_new

"""Plotting"""
def plot_rho(initial, method,N,axs):
    N_index=np.zeros(N)
    for i in range(0,N):
        N_index[i]=i
   
    U=initial(N)
    axs.plot(N_index, U[:,0]) #plot tho
   
    for i in range(N_time):
        U_new=method(U)
       
        if i%100==0:
            axs.plot(N_index, U_new[:,0]) #rho
            
       # if i==N_time-1:
        #    axs.plot(N_index, U_new[:,0])
        U=U_new
     
    axs.set_title(method.__name__)
    
def plot_thE(initial, method,N,axs): #thermal Energy
    N_index=np.zeros(N)
    for i in range(0,N):
        N_index[i]=i
   
    U=initial(N)
    e=(U[:,2]-0.5*U[:,1]*U[:,1]/U[:,0])/U[:,0]
    axs.plot(N_index, e) #plot tho
   
    for i in range(N_time):
        U_new=method(U)
        e_new=(U_new[:,2]-0.5*U_new[:,1]*U_new[:,1]/U_new[:,0])/U_new[:,0]
        if i%100==0:
            axs.plot(N_index, e_new) #rho
            
       # if i==N_time-1:
        #    axs.plot(N_index, U_new[:,0])
        U=U_new
     
    axs.set_title(method.__name__)
 
method=[LAX, Method_B, Method_C]

rows, cols= 2,3

fig1, axs1 = plt.subplots(rows, cols, sharex=True, sharey=False)
fig1.suptitle("Schock Tube density & thermal Energy")

fig2, axs2 = plt.subplots(rows, cols, sharex=True, sharey=False)
fig2.suptitle("Blast Wave density & thermal Energy")

for i, m in enumerate(method):
    plot_rho(Schock_Tube, m, N, axs1[0][i]) #density Plot
    plot_thE(Schock_Tube, m, N, axs1[1][i]) #thermal Energy plot
    plot_rho(Blast_Wave, m, N, axs2[0][i])
    plot_thE(Blast_Wave, m, N, axs2[1][i])
    


