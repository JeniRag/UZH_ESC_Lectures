# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:29:20 2022

@author: Sujeni

W4 Game of Life
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import time
import matplotlib.animation as FuncAnimation
import matplotlib.animation as animation

from scipy import signal

from scipy import misc

Lx=1000
Ly=1000

d={  0: "dead", 1: "alive", 2: "sick"} 

class Node():
    def __init__(self, boundary,level):
        self.level=level
        self.b=boundary
        self.sub=[]
        
    # def add_sub(self):
    #     while self.level <10:
    #         Lx=self.b[0][0] -self.b[1][0]
    #         Ly=self.b[1][0] -self.b[1][1]
            
    #         sn1=Node(np.array([[ self.b[0][0], self.b[0][1]], [self.b[0][0]+ Lx/2, self.b[0][1]+Ly/2]]),self.level+1)
    #         sn2= Node(np.array([[ self.b[0][0] +Lx/2, self.b[0][1] +Ly/2], [self.b[0][0]+ Lx/2, self.b[0][1]+Ly]]),self.level+1)
    #         sn3=Node(np.array([[ self.b[0][0] +Lx/2, self.b[0][1] ], [self.b[0][0]+ Lx, self.b[0][1]+Ly/2]]),self.level+1)
    #         sn4=Node(np.array([[ self.b[0][0] +Lx/2, self.b[0][1]+Ly/2 ], [self.b[1][0]+ Lx, self.b[1][1]]]),self.level+1)
    #         # self.sub=[self.sn1, self.sn2, self.sn3, self.sn4]
    #         self.sub.append(sn1)
    #         self.sub.append(sn2)
    #         self.sub.append(sn3)
    #         self.sub.append(sn4)
            
       

def sub_nodes(n1):
    while n1.level<10:
        Lx=n1.b[0][0] -n1.b[1][0]
        Ly=n1.b[1][0] -n1.b[1][1]
        
        n11 =Node(np.array([[ n1.b[0][0], n1.b[0][1]], [n1.b[0][0]+ Lx/2, n1.b[0][1]+Ly/2]]), n1.level+1)
        n12= Node(np.array([[ n1.b[0][0] +Lx/2, n1.b[0][1] +Ly/2], [n1.b[0][0]+ Lx/2, n1.b[0][1]+Ly]]),n1.level+1)
        n13=Node(np.array([[ n1.b[0][0] +Lx/2, n1.b[0][1] ], [n1.b[0][0]+ Lx, n1.b[0][1]+Ly/2]]),n1.level+1)
        n14=Node(np.array([[ n1.b[0][0] +Lx/2, n1.b[0][1]+Ly/2 ], [n1.b[1][0]+ Lx, n1.b[1][1]]]),n1.level+1)
        
        n1.sub.append(n11)
        n1.sub.append(n12)
        n1.sub.append(n13)
        n1.sub.append(n14)
        
        for ns in n1.sub:
            sub_nodes(ns)
    
    
    
    

n1=Node(np.array([[0,0], [Lx, Ly]]), 0)  
sub_nodes(n1)    
        
        
    
def create_field( Lx, Ly):
    x=np.arange(Lx)
    y=np.arange(Ly)
    return np.zeros((Lx,Ly))
     

def update(f):
   
    kernel=np.array([[1,1,1], [1,0,1], [1,1,1]])  
    a=signal.convolve2d( f, kernel, boundary="wrap", mode="same")
    
    f2=np.copy(f)

    f2[(a==3)]=1 #if 3 nb-> becomes alive
    f2[( a>3)| (a<=1)]=0
    
    return f2
          
def random_cell(f): #populate with alive cells at random positions
    f[random.randint(0,Ly-1)][random.randint(0, Lx-1)] =1 #how to avoid double -> while Loop
 
def drop_horizontal(f, x=random.randint(0,Ly-1) ,y=random.randint(0,Lx-1)):
    while (x<1 or x>len(f) or y > np.shape(f)[1]-2):
        x=random.randint(0,Ly-1)
        y=random.randint(0, Lx-1)
    
    f[x][y:y+3]=1
    
def drop_p1(f,  x=random.randint(0,Ly-1) ,y=random.randint(0,Lx-1)):
    while(x< 3 or x>np.shape(f)[1]-3 or y<3 or y>np.shape(f)[1]-3):
        x=random.randint(0,Ly-1)
        y=random.randint(0, Lx-1)
    
    f[x][y]=1
    f[x-2,y-1:y+2]=1
    f[x+2, y-1: y+2]=1
    f[x-1:x+1,y-2]=1
    f[x-1:x+1,y+2]=1
    
def drop_circle(f,  x=random.randint(0,Ly-1) ,y=random.randint(0,Lx-1)):
    while(x< 3 or x>np.shape(f)[1]-4 or y<3 or y>np.shape(f)[1]-4):
        x=random.randint(0,Ly-1)
        y=random.randint(0, Lx-1)
    
    f[x][y]=1
    f[x-2,y-1:y+2]=1 #top
    f[x+2, y-2: y+3]=1 #bottom
    f[x-2:x+2,y-2]=1 #left
    f[x-2:x+2,y+2]=1 #right
        
def drop_glider(f,  x=random.randint(0,Ly-1) ,y=random.randint(0,Lx-1)):
    while (x<3 or y>np.shape(f)[1]-3):
         x=random.randint(0,Ly-1) 
         y=random.randint(0,Lx-1)
    f[x, y:y+3]=1
    f[x-1, y+2]=1
    f[x-2, y+1]=1
    
       
     
def random_cell_sick(f):
    f[random.randint(0,Ly-1)][random.randint(0, Lx-1)] =2

def game(Lx,Ly, it=100, sick=False, random=True, showPlot=True):
    g=create_field(Lx,Ly)
    f_list=[]
    
    if random==True:
        for i in range(Lx*Ly//2):
            random_cell(g)
         
        if sick==True:
            for j in range(Lx*Ly//4):
                  random_cell_sick(g)
        
        drop_horizontal(g)
        drop_circle(g)
        drop_p1(g)
     

    else:
        drop_p1(g, Lx//2, Ly//2)
        drop_glider(g, 4,4)
        drop_circle(g, Lx-6, Ly-6)

    if showPlot==True:
        plt.ion()
        im=plt.imshow(g, vmin=0, vmax=3)
            
        while True:  
            g=update(g)
    
            im.set_data(g)
            plt.draw()
            plt.pause(0.1) #s
            
    # """to see how many iterations in 1 min"""
    
    
    else: 
        f_list.append(g)
        
        start=time.time()
        iteration=0
        while time.time()-start<60:
            g=update(g)
            iteration+=1
        
        end=time.time()-start
        print("time: ", np.round(end,5), "s")
        f_list.append(g)
        
        print("iterations: ", iteration)
    
    return f_list

def Plot(f_list,s1=0, s=10):
    for i in range(s1,s):
        plt.figure()
        plt.imshow(f_list[i], vmin=0, vmax=3)
        
 
"""main"""

# f_list=game(Lx,Ly,100,sick=True, random=True, showPlot=True)
# plt.ioff()

# Plot(f_list,0,2)


"""1 In 1 min: 1383 iterations for 1000x1000"""