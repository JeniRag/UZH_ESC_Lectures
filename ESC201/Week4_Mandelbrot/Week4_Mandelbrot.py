# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 12:25:36 2020

@author: Sujeni
"""
"""
ESC201-week4-exercise3-Mandelbrot Set
"""
import numpy as np
import matplotlib.pyplot as plt

def Mandelbrot(n,k):
    
    a=np.linspace(-2,1.5,n)
    b=np.linspace(-1.5,2,n)
    
    A,B=np.meshgrid(a,b)
    z0=A+B*1j
    z=np.copy(z0)
    c=np.copy(z0)
    arr=2*np.ones_like(z0, dtype=float) # array with int 2
    rmax=np.maximum(abs(z0), arr) #np.maximum() compares to arrays

    color_seq=np.zeros_like(z0, dtype=float)
    for i in range(k):
        cond=abs(z)<rmax
        color_seq[cond]+=0.1 #Sequenz of color
        z[cond]=z[cond]**2+c[cond]
        
    return z0, color_seq #returns tuple

M_set=Mandelbrot(400,100)

X = [x.real for x in M_set[0]]
Y = [y.imag for y in M_set[0]]

fig, ax=plt.subplots()
ax.set_xlim([-2,1.5])   
ax.set_ylim([-1.5,1.5])
ax.set_xlabel("real part")
ax.set_ylabel("imaginary part")
ax.set_title("Exercise3-Mandelbrot Set")
ax.scatter(X,Y, s=1,  c=M_set[1],cmap='magma')


#plt.show()
plt.savefig("./Mandelbrot.png", dpi=300)

"""
optional task: Julia sets
"""
def Julia(n,k,konst):
    a=np.linspace(-1.5,2,n)
    b=np.linspace(-1.5,1.5,n)
    A,B=np.meshgrid(a,b)
    z0=A+B*1j
    
    z=np.copy(z0)
    
    c=konst*np.ones_like(z0, dtype=complex)
   
    zmax=np.ones_like(z0, dtype=float)
    color_seq=np.zeros_like(z0, dtype=float)
    
    for i in range(k):
        
        cond=abs(z)<zmax
        
        z[cond]=z[cond]*z[cond]+c[cond]
        color_seq[cond]+=0.1
    return z0, color_seq

#Different konstants for Julia set:
konst1=0.355 + 0.355j
konst2=-0.4  -0.59j
konst3=-0.54 + 0.54j
konst4= 0.8 - 0.1j

Julia_set1=Julia(400,300, konst1)   
Julia_set2=Julia(400,300, konst2)
Julia_set3=Julia(400,300, konst3 )
Julia_set4=Julia(400, 300, konst4)


rows, cols = 2, 2 # SR: overall 4 elements
fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True) # They all have the same x and y axes tick-marks. Some more options: figsize=(5.0,5.0), dpi=200
	
fig.tight_layout() # Helps with subplot spacing. Alternatively, use: plt.subplots_adjust

plot1 = axs[0][0]
plot2 = axs[0][1]
plot3 = axs[1][0]
plot4 = axs[1][1]


plot1.set_xlim([-1.5,1.5])   
plot1.set_ylim([-1.5,1.5])
plot1.set_xlabel("real part")
plot1.set_ylabel("imaginary part")

plot1.set_title(konst1)
plot2.set_title(konst2)
plot3.set_title(konst3)
plot4.set_title(konst4)


plot1.scatter(Julia_set1[0].real,Julia_set1[0].imag,s=1, c=Julia_set1[1], cmap='viridis')
plot2.scatter(Julia_set2[0].real, Julia_set2[0].imag, s=1, c=Julia_set2[1], cmap='viridis')
plot3.scatter(Julia_set3[0].real, Julia_set3[0].imag, s=1, c=Julia_set3[1], cmap='viridis')
plot4.scatter(Julia_set4[0].real, Julia_set4[0].imag, s=1, c=Julia_set4[1], cmap='viridis')

#plt.show()	
plt.savefig("./Julia_Set.png", dpi=300)
