# -*- coding: utf-8 -*-
"""
Created on Tue May  3 09:37:53 2022

@author: Sujeni

Computational Thinking
Naive Bayers Classifier
"""
import numpy as np
import matplotlib.pyplot as plt
#%%


#%%
"""generate Data"""
np.random.seed(100)
N=200


mean1 = [0, 3]
cov1 = [[70, 0], [0, 100]]  # diagonal covariance
x1, y1 = np.random.multivariate_normal(mean1, cov1, N).T
x1=x1.astype(int)
y1=y1.astype(int)

z1=np.ones_like(x1, dtype=int)


mean2=[20,-3]
cov2=[[20, 0], [0, 50]]
x2, y2 = np.random.multivariate_normal(mean2, cov2, N).T
x2=x2.astype(int)
y2=y2.astype(int)
z2=np.ones_like(x2, dtype=int)*2

# plt.figure(figsize=[10,10])
plt.figure()
# plt.subplot(1,2,1)
plt.scatter(x1, y1,label="1")
plt.scatter(x2,y2, label="2")
plt.title("training data")
plt.legend()
#plt.xlim(x1.min(), x1.max())
#plt.ylim(x2.min(), x2.max())
#plt.axis('equal')


#%%

D1=np.concatenate([x1.reshape(x1.shape[0],1), y1.reshape(x1.shape[0],1),z1.reshape(z1.shape[0],1)], axis=1)
D2=np.concatenate([x2.reshape(x2.shape[0],1), y2.reshape(x2.shape[0],1),z2.reshape(z2.shape[0],1)], axis=1)

data=np.concatenate([D1,D2])

"""Posteriors"""
def P(c):
    s=np.where(data[:,-1]==c)
    n=len(s[0])
    return n/len(data[:,0])

def P1(x):
    s=np.where(data[:,0]==x)
    n=len(s[0])
    return n/len(data[:,0])

def P2(x):
    s=np.where(data[:,1]==x)
    n=len(s[0])
    return n/len(data[:,1])
 
def P1cond(x,y):
    x_count=len( data[:,0]==x)
    
    c11, count11=np.unique(data[data[:,-1]==y][:,0], return_counts=True)
    
    x_index=np.where(data[:,0]==x)[0]
    p=0
    for i in x_index:
        if data[i,-1]==y:
            p+=1
            
    
    return p/data.shape[0]

def P2cond(x,y):
    x_count=len( data[:,1]==x)
    
    # c11, count11=np.unique(data[data[:,-1]==y][:,0], return_counts=True)
    
    x_index=np.where(data[:,1]==x)[0]
    p=0
    for i in x_index:
        if data[i,-1]==y:
            p+=1
            
    
    return p/data.shape[0]    
    
def Probability(x):
    x1,x2=x[0], x[1]
    k=np.unique(data[:,-1])
    
    P_list=[]
    for y in k:
       
         
        p=P1cond(x1,y)*P2cond(x2,y)*P(y)/(P1(x1)*P2(x2))
        P_list.append(p)
        
    return np.argmax(P_list)+1

    
xmax=data.max()
xmin=data.min()

#%%
"""classifying new data"""
n=800
d1=np.random.choice(data[:,0],size=(n))
d2=np.random.choice(data[:,1], size=(n))
d=np.concatenate([d1.reshape(d1.shape[0],1),d2.reshape(d2.shape[0],1)],axis=1) #np.random.randint(xmin, xmax, size=(200,2))


p=[]
for dx in d:
    p.append(Probability(dx))
# dx1,dx2=np.meshgrid(d,d)

# plt.subplot(2,2,2)
plt.figure()
plt.scatter(d[:,0], d[:,1],  c=p)
plt.title("new data")


