# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:19:29 2021

@author: Sujeni 

ESC202- Exercise 6 - 2D Ising model
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import math

"""Parameters"""
T_start=4
T_end=0.1
T_factor=0.98
 ##########
kB=1
J=1
Tc=J/(kB*0.4469)
Nx=50
Ny=50

class particle:
    def __init__(self, r,s):
        self.r=r
        self.spin=s
        self.H=None
        
        if self.spin==1:
            self.color='black'
        else:
            self.color='white'
    def set_color(self):
        if self.spin==1:
            self.color='black'
        elif self.spin==-1:
            self.color='white'
        
            
def create_particles(Nx, Ny):
    N=Nx*Ny
    particle_array=np.zeros(shape=(Nx,Ny), dtype=particle)
    for i in range(0,Nx):
        for j in range(0,Ny):
            r=np.array([i,j])
            s=random.randrange(-1,2,2) #random number -1 or 1
            particle_array[i][j]=particle(r,s)
           # particle_array=np.append(particle_array, p)
    return particle_array

def Energy(array):
    
    for i in range(0,len(array)):
        for j in range(0,len(array[0])):
            sum=0.0
            up=(j-1)%Ny
            down=(j+1)%Ny
            
            right=(i+1)%Nx
            left=(i-1)%Nx
            
            sum+=array[i][up].spin+array[i][down].spin
            sum+=array[right][j].spin+array[left][j].spin
            array[i][j].H=2*J*array[i][j].spin*sum

def Mag(spin_list):
    
    m=spin_list.sum()
    #print(m)
    return m/(Nx*Ny)

def extract_spin(array):
    spin_list=np.zeros(shape=(Nx,Ny),dtype='int')
    for i in range(0,Nx):
        for j in range(0,Ny):
            spin_list[i][j]=array[i][j].spin
    return spin_list

def extract_color(array):
    color_list=[]
    for i in range(0,Nx):
        for j in range(0,Ny):
            
           color_list.append(array[i][j].color)
    return color_list
    
def extract_coordinate(array):
    X=[]
    Y=[]
    for i in range(0,Nx):
        for j in range(0,Ny):
           X.append(array[i][j].r[0])
           Y.append(array[i][j].r[1])
    return X,Y
    
def ising_model(Nx,Ny, ax):
    array=create_particles(Nx,Ny)
    X,Y=extract_coordinate(array)
    
    
    #plot_init(array,ax)
    
    T=T_start
   
    spin_list=extract_spin(array)
    m=Mag(spin_list)
    m_list=[m]
    T_list=[T]
    c=[extract_color(array)]
    count=0
    
    while T>T_end:
        count+=1
        beta=1.0/(kB*T)
        Energy(array)
        r=random.random()
        #print(r)
        for i in range(0,Nx):
            for j in range(0,Ny):
                p=array[i][j]
                
                if p.H<0:
                    p.spin *= -1
                #r=random.random()
                
                elif r<math.exp(-beta*p.H) and math.exp(-beta*p.H)<1:
                        p.spin *=-1
                
                p.set_color()
        
        T_list.append(T)
        spin_list=extract_spin(array)
        m=Mag(spin_list)
        m_list.append(m)
        color_list=extract_color(array)
        c.append(color_list)
        
        #print("T="+str(T)+","+" <m>="+str(m))
        T*=T_factor
    #for p in array[1]:
    #    print(p.H)
    plot_init(array,ax)
    print(count)
    return X,Y, c, T_list, m_list          
                
def plot_init(array,ax):
    
    N=Nx*Ny
    coordinates=np.zeros(N, dtype=particle)
    for i in range(0,Nx):
        for j in range(0,Ny):
            ax.scatter(array[i][j].r[0],array[i][j].r[1],marker='s', color=array[i][j].color)
    ax.set_title("Ising-Model: Nx, Ny = "+str(Nx)+"," + str(Ny))

def plot_magnetisation(T_list, m_list, ax):
    ax.set_title('Magnetization')
    ax.set_xlabel('T')
    ax.set_ylabel('<m>')
    ax.plot(T_list, m_list ,'o-')
    
def update(frame,plot,ax):
 
        plot[0].remove()
        plot[0]=ax.scatter(X,Y, marker='s', color=color_list[frame])
        #print(array_list[frame])

    
fig=plt.figure()
ax=fig.add_subplot(111)

fig2=plt.figure()
ax2=fig2.add_subplot(111)

fig3=plt.figure()
ax3=fig3.add_subplot(111)
ax3.set_title("Ising-Model: Nx, Ny = "+str(Nx)+"," + str(Ny))

X,Y,color_list, T_list, m_list=ising_model(Nx,Ny,ax)

p1=[ax.scatter(X,Y, marker='s', color=color_list[0])]

anim = animation.FuncAnimation(fig3, update, fargs=(p1,ax3),
                               frames=len(color_list)-1, interval=100, repeat=False) #, blit=False)

plot_magnetisation(T_list, m_list, ax2)
anim.save('Exercise6_Ising-model.mp4', fps=2,dpi=200, extra_args=['-vcodec', 'libx264'])
plt.show()  
          