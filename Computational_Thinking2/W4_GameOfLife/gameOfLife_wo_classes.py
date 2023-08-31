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

Lx=20
Ly=20

d={  0: "dead", 1: "alive", 2: "sick", 100: "ignore border"} 


def create_field( Lx, Ly):
    x=np.arange(Lx)
    y=np.arange(Ly)
    return np.zeros((Lx,Ly))
     
def NB(f):
    # print(f)
    # print(" ")
    border=100
    bottom= np.roll(f, -1, axis=0) #######how to handle border
    
    # bottom[-1, :]=border
    
   
    top=np.roll(f ,1, axis=0)
    # top[0,:]=border
    
    
    right=np.roll(f, -1, axis=1)
    # right[:,-1]=border
    
    left=np.roll(f ,1, axis=1)
    # left[:,0]=border
    
    drb =np.roll(bottom, -1, axis=1) #diagonal right bottom
    # drb[:,-1]=border
   
    
    dlb= np.roll(bottom, 1, axis=1) #diagonal left bottom
    # dlb[:,0]=border
    
    
    drt= np.roll(top, -1, axis=1) #diagonal right top
    # drt[:,-1]=border

    
    # print(drt)
    dlt=np.roll(top, 1, axis=1)
    # dlt[:,0]=border
   
    
     # print(bottom)
     # print(top)
     # print(right)
     # print(left)
     # print(drb)
    # print(dlb)
     # print(dlt)
     
    dead=np.zeros_like(f)
    alive=np.zeros_like(f)
    sick=np.zeros_like(f)
    
    nb=[bottom, top, right, left, drb, dlb, drt, dlt]
    for n in nb:
        alive[(n==1) | (n==2) | (n==3) |(n==4)] +=1
        dead[(n==0)]  +=1  
        sick[(n==2) | (n==3) |(n==4) ] +=1
        
    return alive, dead,sick
    
def update(f):
   
    a,d,s=NB(f)
    f2=np.copy(f)
    """check if it is not executed after"""
    f2[(a==3)]=1 #if 3 nb-> becomes alive
    f2[((f==1)&(s>0))]=2 #if at least one neighbour of an alive cell is sick, it becomes sick  "
    f2[(f==2)] =3
    f2[(f==3)] =4
    
    f2[((a>3) |(a<=1) | (f==4))]=0 #more than 3 alive nb or less than 1 nb or cell is infected> cell dies
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

def game(Lx,Ly, it=100, sick=False, random=True):
    g=create_field(Lx,Ly)
    f_list=[]
    
    if random==True:
        for i in range(Lx*Ly//2):
            random_cell(g)
         
        if sick==True:
            for j in range(3):
                  random_cell_sick(g)
        
        drop_horizontal(g)
        drop_circle(g)
        drop_p1(g)
     
    else:
         drop_p1(g, Lx//2, Ly//2)
         drop_glider(g, 4,4)
         if sick==True:
             g[6][6]=2

    start=time.time()
    
    f_list.append(np.copy(g))
    for i in range(it):
        
        g=update(g)
        f_list.append(np.copy(g))
    
    end=time.time()-start
    print("time: ", np.round(end,5), "s")

    return f_list

def Plot(f_list,s1=0, s=10):
    for i in range(s1,s):
        plt.figure()
        plt.imshow(f_list[i], vmin=0, vmax=2)
 
"""main"""

f_list=game(Lx,Ly,1000,sick=True, random=False)

# Plot(f_list,0,15)


fig=plt.figure()

def updatefig(i):
    
    plt.imshow(f_list[i], vmin=0, vmax=2)
    plt.title("iteration: "+str(i))
    return 

anim = animation.FuncAnimation(fig, updatefig, frames=len(f_list), interval=20, blit=False) #interval=delay between frames in ms
plt.show()

# anim.save("animation_gameOfLife.mp4", writer=animation.FFMpegWriter(fps=10))

