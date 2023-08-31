# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 07:39:11 2022

@author: Sujeni
"""
import numpy as np
import matplotlib.pyplot as plt

d={"B4F440": "Neanderthal", "O03169": "Coelacanth", "P00850":"Fruit fly", " Q33823":"Starfish", "Q95A26": " Bornean orangutan",
   "F8RBX8": "Pyhton","O21402": "Ostrich","P33507": "Malaria mosquito", "Q35920": "Atlantic salmon", "Q9T9W0" : "Chimpanze",
   "H6TQD1": "Great fruit-eating bat", "P00846": "Human","Q2I3G9": "Indian elephant", "Q38PR7": "Siberian woolly mammoth", "Q9TA24":"African elephant" }

w=3 #windowsize

def file_open(s1):
    f=open(s1, "r")
    s1=""
    lines=f.readlines()
    #print(lines)
    for i in range(1,len(lines)):
        s1=s1+lines[i]
        
    s1=s1.replace("\n", "") #Do I have to remove the linebreakers??
    f.close()
    
    return s1

g1="B4F440"
g2="P00846"

def dots(g1, g2, w, title="DotPlot"): #string 1, string2, windowsize

    f1="./Sequenz/"+g1+".fasta"
    f2="./Sequenz/"+g2+".fasta"
    
    s1=file_open(f1)
    s2=file_open(f2)

    # if len(s1)>len(s2):
    #     diff=len(s1)-len(s2)
    #     s2=s2+diff*" "
    
    # if len(s1)<len(s2):
    #     diff= len(s2)-len(s1)
    #     s1=s1+diff*" "
     #Do I have to make the strings to same length?  
        
    d={} #dictionary
    
    for i in range(len(s1)-w):
        d[s1[i:i+w]]=i
        
    k=list(d.keys())

    v=np.array(list(d.values()))

    s2_array=np.zeros(len(k))
    
    for i in range(len(k)-w):
        if s2[i:i+w] in k:
            s2_array[i]= d[s2[i:i+w]] #return position in k list
            
    plot(g1,g2, v, s2_array,title)
       
    return v, s2_array


def plot(g1,g2,a1,a2,title):
    zero=(a2==0)
    plt.figure()
    plt.scatter(a1[~zero],a2[~zero]) #omit plotting values at 0
    plt.title("Dot-Plot \n Windows Size="+str(w))
    plt.xlabel( "frame shft \n" + d[g1])
    plt.ylabel(d[g2]+ "\n frame shift")
    plt.savefig("./Bilder/"+title+".png")
    
dots(g1,g2,w, "humans")
dots("Q2I3G9", "Q9TA24",w, "elephant-m")
dots("Q2I3G9", "Q38PR7",w, "elephant")
dots("P00850", "P00846", w, "Fruchtfliege-Mensch")
# dots("B4F440", "Q95A26",w, "Neandertal-Orangutan")
# dots("P00850", "P33507", w, "Fruchtfliege-Moskito")
