# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:04:14 2022

@author: Sujeni
"""
import matplotlib.pyplot as plt
import numpy as np
import time
import random

"""with N=6 and 10000 iterations-> sometimes I get 35"""
# Niter=15100

random.seed(100)

class cube():
    def __init__(self,N): #number of cubes
        self.N=N

        self.field=np.full((N,N),0)
        self.poss=[] #possiblities lists
        self.plot_poss=[] #to plot matrixes of same size
        self.variations=[]
        
        self.up=np.array([0, -1])
        self.left= np.array([-1,0])
        self.down=np.array([0,1])
        self.right=np.array([1,0])
        
        self.l=[self.up, self.down, self.left,self.right]

    
    def show_poss(self):
        for i in range(len(self.poss)):
            print(self.poss[i])
            print(" ")
           
    def attach_number(self, field, ix, iy, j,i, indexes, s):

        x=ix+self.l[j][0]
        y=iy+self.l[j][1]
     
        
        if ( ((x<0 or x>=self.N or y<0 or y>=self.N ) or field[x][y]==1 ) and (s<4) ):
            return self.attach_number(field, ix,iy,(j+1)%4, i,indexes, s+1) #try the other directions too
           
        elif ( ((x<0 or x>=self.N or y<0 or y>=self.N) or field[x][y]==1) and (s>=4 ) ): #if all 4 directions not possible, try new position
              
            i=(i+1)%len(indexes[0])
            
            ix=indexes[0][i]
            iy=indexes[1][i]
            
            return  self.attach_number(field, ix,iy,np.random.randint(4),i,indexes, 0)

        
        return x, y

    def cut(self,a): #to avoid that different locations in the matrix matter
        i=np.where(a==1)
        min_x=np.min(i[0])
        max_x=np.max(i[0])
        min_y=np.min(i[1])
        max_y=np.max(i[1])
        
        return a[min_x:max_x+1, min_y:max_y+1]
  
    def combinations(self, Niter):

         #set a starting point
        #number of cubes
         #create random 6x2 arrays avoiding repetition
        #check if they follow the rules
        #draw
        #save shape in dictionary or list with many possible transformations
        
        
        for n in range(Niter):

            field=np.zeros_like(self.field)
            
            x_start=np.random.randint(self.N)
            y_start=np.random.randint(self.N)
          
          
            field[x_start][y_start]=1 #what if in edge-> only 2 possibilities

            assert(len(np.where(field>0)[0])==1)
            C=1
            
            while C<self.N:   
                indexes=np.where(field>0)
                p=np.random.permutation(len(indexes[0])) #for permutation of indexes
                indexes1=indexes[0][p]
                indexes2=indexes[1][p]
                indexes=np.concatenate((indexes1, indexes2))
                indexes=indexes.reshape(2, len(indexes1))
                
                i=np.random.randint(len(indexes)-1)
                
                ix=indexes[0][i]
                iy=indexes[1][i]
                
                j=np.random.randint(4) #random direction
                x,y =self.attach_number(field, ix,iy, j,i,indexes, 0)

                assert field[x][y] != 1 # check if uncoccupied
                
                field[x][y]=1
                C+=1
            
            assert len(np.where(field==1)[0])==self.N #check if there are exactly desired number of cubes
            
            field1=self.cut(field)
            
            if any(np.array_equal(field1, i) for i in self.variations)==False:
                self.poss.append(field1)
                self.plot_poss.append(field)
                self.vary(field1)
        
    def vary(self, poss): #variations
        
        poss1=np.rot90(poss)
        poss2=np.rot90(poss1)
        poss3=np.rot90(poss2)
        
        flip=np.fliplr(poss)
        
        poss4=np.rot90(flip)
        poss5=np.rot90(poss4)
        poss6=np.rot90(poss5)
        
        self.variations.append(poss)
        self.variations.append(poss1)
        self.variations.append(poss2)
        self.variations.append(poss3)
        self.variations.append(poss4)
        self.variations.append(poss5)
        self.variations.append(poss6)
        self.variations.append(flip)
        
        return
    
    
def plot(cube): 
    assert len(cube.plot_poss)==len(cube.poss)
    figure = plt.figure(figsize=(10,10))
    
    """must be changed for other cube sizes"""

    L=len(cube.plot_poss)  
    n=np.ceil(np.sqrt(len(cube.plot_poss)))
    H=int(n)
    W=int(n)

    for i in range(1,L+1): 
        plt.subplot(H,W, i)
        plt.axis("off")
        plt.imshow(cube.plot_poss[i-1])
        
    
def main(N,Niter, p=True):
    CN=cube(N)
    
    start=time.time()
    CN.combinations(Niter)
    end=time.time()-start
    
    print("time: ", end, " s")
    print("combinations found for " , N, "cubes: ", len(CN.poss))
    
    if p==True:
        plot(CN)
    
    return CN
    
#C5=main(5, 1000)
#C6=main(6, 6000) #1.5 s
C7=main(7, 15500) #15000 iterations ->108, 6s
# C8=main(8, 80000) #57 s, 361



