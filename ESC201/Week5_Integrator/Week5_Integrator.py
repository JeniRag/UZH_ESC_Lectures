# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 14:11:40 2020

@author: Sujeni
"""
"""
ESC 201- week5- ex.4
"""
import matplotlib.pyplot as plt
import numpy as np
    
"""Parameters"""
km=2
kmf=0.02
kfm=0.01
kf=1.06
h = 0.01
nSteps = 10000

"""initial values"""
t0 = 0
y0 =np.array([100, 15]) #initial condition at time=0. (mice, fox)


"""Lotka-Voltera-Model"""    
def LoVo (t,y): #Lotka-Voltera-Model
    
    m, fox = y[0], y[1]
    
    dm=(m*km-kmf*m*fox) #derivative
    df=(-kf*fox+kfm*fox*m) #derivative
    return np.array([dm, df])

"""Methods to solve the ODE"""
def MidPointRK (tn, yn, h, LoVo): # Midpoint Runge-Kutta
    ynew = yn + h* LoVo(tn + 1/2*h, yn + h/2*LoVo(tn, yn))
    return ynew

def eulerStep(tn, yn, h, LoVo):
    ynew = yn + h * LoVo(tn,yn)
    return ynew

def fourthorderRK(tn, yn, h, LoVo): # 4th order Runge-Kutta
    k1=h*LoVo(tn,yn)
    k2=h*LoVo(tn+h/2, yn+k1/2)
    k3=h*LoVo(tn+h/2, yn+k2/2)
    k4=h*LoVo(tn+h, yn+k3)
    ynew=yn+k1/6+k2/3+k3/3+k4/6 
    return ynew

def odeSolver(t0, y0, dfFunc, h, nSteps, solverStepFunc):
    """This is a general ODE solver that takes the
    derivative df/dt (dfFunc) and the algorithm for one time
    step (solverStepFunc) as function arguments.
    """
    yn = y0
    tn = t0
    tlist = [t0]
    ylist = [y0]
    for n in range(nSteps):
        ynew = solverStepFunc(tn, yn, h, dfFunc)
        tn += h
        tlist.append(tn)
        ylist.append(ynew)
        yn = ynew
    return (tlist, np.array(ylist))

t, y1 = odeSolver(t0, y0, LoVo, h, nSteps, eulerStep)
t, y2 = odeSolver(t0, y0, LoVo, h, nSteps, MidPointRK)
t, y3 = odeSolver(t0, y0, LoVo, h, nSteps, fourthorderRK)

#different Initial value with fourthorderRK
y01=np.array([10,34])
y02=np.array([200,150])
y03=np.array([125,100])
y04=np.array([30,200])
y05=np.array([106,100]) #Fixed Point

t, y41 =  odeSolver(t0, y01, LoVo, h, nSteps, fourthorderRK)
t, y42 =  odeSolver(t0, y02, LoVo, h, nSteps, fourthorderRK)
t, y43 =  odeSolver(t0, y03, LoVo, h, nSteps, fourthorderRK)
t, y44 =  odeSolver(t0, y04, LoVo, h, nSteps, fourthorderRK)
t, y45 =  odeSolver(t0, y05, LoVo, h, nSteps, fourthorderRK) 


"""Population over time plotting"""
rows, cols = 2, 2 # SR: overall 4 elements
fig, axs1 = plt.subplots(rows, cols, sharex=True, sharey=True)


plot1 = axs1[0][0]
plot2 = axs1[0][1]
plot3 = axs1[1][0]
plot4 = axs1[1][1]

for n in range(y1.shape[1]):
    plot1.plot(t, y1[:, n], label="y1"+str(n))
for n in range(y2.shape[1]):
    plot2.plot(t, y2[:, n], label="y2"+str(n))
for n in range(y3.shape[1]):
    plot3.plot(t, y3[:, n], label="y3"+str(n))
for n in range(y45.shape[1]):
    plot4.plot(t, y45[:, n], label="y45"+" "+str(n))
    
plot1.legend( loc='upper left')
plot1.set_title("Forward Euler Method")
plot1.set_xlabel('t')
plot1.set_ylabel('Population')

plot2.legend(loc='upper left')
plot2.set_title("Midpoint Runge-Kutta")
plot2.set_xlabel('t')
plot2.set_ylabel('Population')

plot3.legend(loc='upper left')
plot3.set_title("4th order Runge-Kutta")
plot3.set_xlabel('t')
plot3.set_ylabel('Population')

plot4.legend(loc='upper left')
plot4.set_title("4th order Runge-Kutta with Fixed Point as y0")
plot4.set_xlabel('t')
plot4.set_ylabel('Population')


plt.savefig("./Population.png", dpi=300)

"""Phase Diagramm Plotting"""
rows, cols = 2, 2 # SR: overall 4 elements
fig, axs2 = plt.subplots(rows, cols, sharex=True, sharey=True)

plot1p = axs2[0][0]
plot2p = axs2[0][1]
plot3p = axs2[1][0]
plot4p = axs2[1][1]

plot1p.set_title("Phasediagramm-Forward Euler Method")
plot1p.set_xlabel('mouse')
plot1p.set_ylabel('fox')

plot2p.set_title("Phasediagramm-Midpoint Runge-Kutta")
plot2p.set_xlabel('mouse')
plot2p.set_ylabel('fox')

plot3p.set_title("Phasediagramm-4th order Runge-Kutta")
plot3p.set_xlabel('mouse')
plot3p.set_ylabel('fox')

plot4p.set_title("Phasediagramm- different initial values")
plot4p.set_xlabel('mouse')
plot4p.set_ylabel('fox')


m1, f1 = zip(*y1)
m2, f2 = zip(*y2)
m3, f3 = zip(*y3)

m41, f41 = zip(*y41)
m42, f42 = zip(*y42)
m43, f43 = zip(*y43)
m44, f44 = zip(*y44)
m45, f45 = zip(*y45)

plot1p.scatter(m1,f1, s=1,c='Purple')
plot2p.scatter(m2,f2, s=1,c='Purple')
plot3p.scatter(m3,f3, s=1, c='Purple')
 
plot4p.scatter(m41,f41, s=1, c='Purple', label=y01)
plot4p.scatter(m42,f42, s=1, c='Blue', label=y02)
plot4p.scatter(m43,f43, s=1, c='Green', label=y03)
plot4p.scatter(m44,f44, s=1, c='Red', label=y04)
plot4p.scatter(m45,f45, s=1, c='Black', label=y05)
plot4p.legend()


plt.savefig("./Phasediagram.png", dpi=300)



