# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 14:01:30 2021

@author: Sujeni

ESC202- Exercise 7- Traveling Salesman

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import copy
import math


""" Notes#############
tour =[0,1,2,... N-1]

E=0
for i in tour:
    E += sqrt(x[i]-x[i-1]**2..)

Swapping 2 cities:
    Choose a random pair i,j
    swap tour[i] with tour[j]
    
    tour=[...a, b, B, ... ,C,c, d ...]
                i      j
    new tour = [... a, c,B, ..,C,b, d...]
    
    d() is distance functino with sqrt
    
    E_old=...+distance(a,b)+..d(b,B)+...+d(C,c)+..d(c,d)+...
    E_new=d...+istance(a,c)+...d(c,B)+...+d(C,b) +d(b,d)+...
   deltaE = Enew - Eold = d(a,c) + d(c,B) + d(C,b) + d(b,d) - d(a,b) - d(b,B) - d(C,c) - d(c,d)
   
REverse a segment:
    tour=[...a,| b...,c,| d ...]
             i       j
    new_tour= tour=[...a,| c...,b| d ...]
    
    dE=d(a,c)+d(b,d)-d(a,b)-d(c,d)
            
"""

"""Parameters"""
N=10
T_factor=0.9

class tour:
    def __init__(self, number,r):
        self.n=number
        self.r=r #coordinates
        self.left=None
        self.right=None
        
        
"""TSP File"""


# Read instance header

def read_tour(filename):
    file = open(filename, 'r')
    Name = file.readline().strip().split()[1] # NAME
    FileType = file.readline().strip().split()[1] # TYPE
    Comment =file.readline().strip().split()[1] # COMMENT
    Dimension = file.readline().strip().split()[1] # DIMENSION
    EdgeWeightType = file.readline().strip().split()[1] # EDGE_WEIGHT_TYPE
    Display=file.readline().strip().split()[1] 
    Node=file.readline().strip().split()

    file.readline()

    N = int(Dimension)-1
    
    positions=np.zeros(shape=(N,2), dtype=float)
    
    tsp_array=np.zeros(N, dtype=tour)

    for i in range(0, N):
       
        z,x,y = file.readline().strip().split()[0:]     
        
        positions[i][0]=float(x)
        positions[i][1]=float(y)
        
        tsp_array[i]=tour(z, positions[i])
        tsp_array[i].n=int(z)
        
    file.close()
    return tsp_array
    
"""generate Tour"""
def generate_tour(N):
    positions=np.zeros(shape=(N,2), dtype='float')
    Tour_list=np.zeros(N, dtype=tour)
    
    for i in range(0,N):
        positions[i][0]=random.random()
        positions[i][1]=random.random()
        
        Tour_list[i]=tour(i, positions[i])
        Tour_list[i].left=(i-1)%N
        Tour_list[i].right=(i+1)%N
    
    return Tour_list

def Energy(array):
    sum=0.0
    for t in array:
        
        sum+=np.sqrt((t.r[0]-t.left.r[0])**2+(t.r[1]+t.left.r[0])**2 )#or use maybe np.roll function
    
    return sum
 

def distance(array,i,j):
    sum=(array[i].r[0]+array[i].r[1])**2+(array[i].r[1]+array[j].r[1])**2
    return sum
    
def dE(array, array_new, i, j):
    E_new=distance(array_new, (i-1)%N, i)+distance(array_new, i,(i+1)%N)+distance(array_new,j, (j-1)%N)+distance(array_new, j, (j+1)%N)
    E_old=distance(array, (i-1)%N, i)+distance(array, i,(i+1)%N)+distance(array,j, (j-1)%N)+distance(array, j, (j+1)%N)
     #maybe save E_old from previous calculations and reuse
    return E_new-E_old

def dE_segment(array, array_new, i, j):  
    E_new=distance(array_new, (i+1)%N, i)+distance(array_new, j, (j+1)%N)
    E_old=distance(array, (i+1)%N, i)+distance(array, j, (j+1)%N) #####check if indices are right
    
     #maybe save E_old from previous calculations and reuse
    
    return E_new-E_old

def swap_tour(array,i,j):
    array_new=copy.deepcopy(array)

    array_new[j]=array[i]
    array_new[i]=array[j]
    
    return array_new

def reverse_segment(array,i,j,T):
    
    
    end=j
    
    array_new=copy.deepcopy(array)
    for x in range(i,j+1):
       
        array_new[x]=array[end]
        
        array_new[end]=array[x]

    
        end-=1
    
    return array_new    
    
def main(array):
    
    Nend=200
    N=len(array)-1
    array_list=[array]
    T=Initial_T(array)
    
    count=0
    count_array=np.zeros(Nend)
    T_array=np.zeros(Nend, dtype=float)
    
    for n in range(0,Nend):
        #print("count")
        T_array[n]=T
        count_array[n]=count
        i=random.randint(0,N-1)
        j=random.randint(0,N-1)
        #array_new=copy.deepcopy(array)
        
        while i==j:
            j=random.randint(0,N-1)
        
        if i>j:
            j,i = i,j
      
        array_new=reverse_segment(array,i,j, T)

        deltaE=dE_segment(array, array_new, i, j)
    
        if deltaE<0:
            array=array_new
            
        
        elif deltaE>0:
             
             term=math.exp(-deltaE/T)
             if random.random()<term:
                 array=array_new
                 
      
            
        T*=T_factor
            
        
        count+=1
        array_list.append(array)   
        
    return array_list, count_array, T_array
             
        
    
def plot(array,ax):
       x=[t.r[0] for t in array]
       y=[t.r[1] for t in array]
      
       ax.plot(x,y, '-o')

def Initial_T(array):
    kB=1
    
    T0=0
    E0=0
    for i in range(0,100):
       T=T0
       i=random.randint(0,N-1)
       j=random.randint(0, N-1)
       
       array_new=reverse_segment(array, i,j,T)
       deltaE=abs(dE_segment(array, array_new, i, j))
       
       if (deltaE)>E0:
           E0=(deltaE)
           T0=(deltaE)*2/3*kB
           
    return T0
        
def order(array):
    for p in array:
        print(p.n)
       
"""test"""

#array=read_tour("berlin52.tsp")
array=read_tour("gr229.tsp")
#array=generate_tour(N)

array1=array

array_matrix, count, Temp=main(array1)

"""Images"""



fig3=plt.figure()
ax3=fig3.add_subplot(111)


x=[t.r[0] for t in array_matrix[0]]
y=[t.r[1] for t in array_matrix[0]]
z=[t.n for t in array_matrix[0]]


line, = ax3.plot(x,y, '-o', color='black', markerfacecolor='white',  markeredgecolor='gray')
text_step=ax3.text(60,-130, "step: "+ str(count[0])+ "/"+ str(count[len(count)-1]))
text_T=ax3.text(60,-150, "T:"+str(round(Temp[0],3)))

for z1, x1, y1 in zip(z, x,y):
    
    label = str(z1)
    ax3.annotate(label, # this is the text
                 (x1,y1), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', color='red', size='8') #

    
def animate(frame):
     x=[t.r[0] for t in array_matrix[frame]]
     y=[t.r[1] for t in array_matrix[frame]]
     
     line.set_ydata(y)
     line.set_xdata(x)
     
     text_step.set_text("step: "+ str(count[frame])+ "/"+ str(count[len(count)-1]))
     text_T.set_text( "T:"+str(round(Temp[frame],3)))
     
     return line,

anim = animation.FuncAnimation(fig3, animate,
                               frames=len(array_matrix)-1, interval=10, repeat=True, blit=False) 
    
#anim.save('Exercise7-Traveling_Salesman.mp4', fps=10,dpi=200, extra_args=['-vcodec', 'libx264'])