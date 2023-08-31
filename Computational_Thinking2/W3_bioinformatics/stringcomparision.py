# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 07:39:11 2022

@author: Sujeni
"""
import numpy as np
import matplotlib.pyplot as plt


s2="Du musst die Leute nehen wie sie sind, andere gibts nicht"
s1="Du muesch d'Lüüt neh wie sie sind, anderi gits net"
w=3

def dots(s1, s2, w): #string 1, string2, windowsize

    # if len(s1)>len(s2):
    #     diff=len(s1)-len(s2)
    #     s2=s2+diff*" "
    
    # if len(s1)<len(s2):
    #     diff= len(s2)-len(s1)
    #     s1=s1+diff*" "
     #Do I have to make the strings to same length?  
        
    d={} #dictionary
    
    #s1_array=np.zeros(len(s1))
    
    for i in range(len(s1)-w):
        d[s1[i:i+w]]=i
        
        #print(s1[i:i+5])
        
    k=list(d.keys())

    #print(k)
    v=np.array(list(d.values()))

    s2_array=np.zeros(len(k))
    
    for i in range(len(k)-w):
        #print(s2[i:i+w])
        if s2[i:i+w] in k:
            s2_array[i]= d[s2[i:i+w]] #return position in k list
       
    return v, s2_array

a1,a2=dots(s1,s2,w)


zero=(a2==0)
plt.figure()
plt.scatter(a1[~zero],a2[~zero]) #omit plotting values at 0
plt.title("Dot-Plot \n Windows Size="+str(w))
plt.xlabel( "frame shift \n"+ s1)
plt.ylabel(s2+ "\n frame shift")
# plt.savefig("MyDotPlot.png")