# -*- coding: utf-8 -*-
"""
Created on Sat Apr  7 20:49:14 2018
Solves the partial differential equation for heat in a current-carrying wire
using the Euler forward method. Includes functions to generate a mesh of position 
points, to generate an A matrix(populateAline) a and calculate the next time step
(calculatenexttime), as well as propagate through time to convergence (stepthrough). 
@author: InfiniteJest
"""
import numpy as np

kappa = 174   #thermal conductivity
rhozero = 5.6*10**-8
rho = 19300     #tungsten density
alpha = 0
c = 132     #tungsten heat capacity
L = 0.12    #length of wire
r = 0.0005  #radius of wire
Azero = np.pi*r**2  #Cross sectional area of wire
Izero = 1   #current

def generate_mesh(intx, step):
    """
    Set of mesh points for the position in the wire.
    """
    xrange = abs(intx[1] - intx[0])    #range of x values
    hx = xrange/float(step)     #change in x at each step
    xpoints = []
    for i in range(step-1):
        x1 = float(intx[0] + (i+1) * hx)   #set up x points
        xpoints.append(x1)
    return xpoints, hx

xpoints, hx = generate_mesh([-0.06,0.06], 10)

def AnalyticalT(x, Izero):
    """
    The analytic solution.
    """
    coef = np.sqrt(alpha*rhozero/kappa)*(Izero/Azero)
    return (np.cos((coef/2)*x)/np.cos(coef*0.12) - 1)/alpha

def AnalAlphaZero(x, Izero):
    """
    The analytic solution at alpha = 0
    """
    coef = ((Izero**2)*rhozero/(8*(Azero**2)*kappa))
    return coef*(0.12**2-4*x**2)


def populateAline(A, hx, index):
    """
    Sets up a line of the A matrix given an index and an hx, change in x.
    """
    Aline = np.zeros(np.shape(A)[-1])
    if index != 0:   #Set T = 0 at lower boundary
        Aline[index-1] = (1/(rho*c))*((kappa/hx**2))
    
    Aline[index] = (alpha*rhozero/((Azero**2)*rho*c))*Izero**2 - 2*(kappa)/(rho*c*hx**2)

    if index != np.shape(A)[-1]-1:  #Set T = 0 at higher boundary
        Aline[index+1] = (1/(rho*c))*((kappa/hx**2))
    return Aline

def calculatenexttime(A, Tarray, time, deltat):
    """
    Given matrix A, a temperature array of length x points, some time stamp, and a deltat,
    calculates the next time step using the Euler forward method.
    """
    Tnp1 = Tarray + deltat*(np.matmul(A, Tarray)) + rhozero*(Izero**2)/((Azero**2)*rho*c)
    time += deltat
    return Tnp1, time

def stepthrough(tsteps, Tinit, xpoints, hx, deltat = None):
    """
    Steps through time in tsteps using the Euler forward method over xpoints. 
    Supply an initial temperature array Tinit of length xpoints. If no deltat
    given, will calculate deltat based on the solution when alpha = 0. Sets up 
    an A matrix and returns the temperature array at the last time step.
    """
    time = 0
    A = np.zeros((len(xpoints), len(xpoints)))  #set up A matrix
    for i in range(len(A)):
        A[i] = populateAline(A, hx, i)  #populate matrix based on Lagrange interpolation

    if deltat is None:
        bestdeltaalpha0 = (hx**2)/(2*(kappa/(rho*c)))   #chooses the solution to delta t for alpha = 0 if none given
        deltat = bestdeltaalpha0
   
    T = Tinit
    times = [time]
    Ts = [T]
    for i in range(tsteps):
        T, time = calculatenexttime(A, T, time, deltat) #Step through
        Ts.append(T)
        times.append(time)
    return Ts[-1]


analytic = np.vectorize(AnalyticalT)
analalphazero = np.vectorize(AnalAlphaZero)

solution = analytic(xpoints, Izero)
solutionzero = analalphazero(xpoints, Izero)

#Test the solution matches at alpha = 0
tstep = 1000
Tinit = np.zeros(len(xpoints))
alphazeroT = stepthrough(tstep, Tinit, xpoints, hx)

import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(xpoints, alphazeroT)
plt.plot(xpoints, solutionzero)
fig.suptitle('Alpha = 0 Test')
plt.xlabel('Position (m)')
plt.ylabel('Temperature (K)')
fig.savefig('Alpha0Test.jpg')


#Calculate and plot the temperature as a function of current
alpha = 0.0045
TI0_1 = stepthrough(tstep, Tinit, xpoints, hx)
TI0_1 = max(TI0_1)

Izero = 2
TI0_2 = stepthrough(tstep, Tinit, xpoints, hx)
TI0_2 = max(TI0_2)

Izero = 5
TI0_5 = stepthrough(tstep, Tinit, xpoints, hx)
TI0_5 = max(TI0_5)

Izero = 10
TI0_10 = stepthrough(tstep, Tinit, xpoints, hx)
TI0_10 = max(TI0_10)

Izero = 15
TI0_15 = stepthrough(tstep, Tinit, xpoints, hx)
TI0_15 = max(TI0_15)

currents = [1, 2, 5, 10, 15]
temperatures = [TI0_1, TI0_2, TI0_5, TI0_10, TI0_15]

fig = plt.figure()
plt.plot(currents, temperatures)
fig.suptitle('Max Temperature vs. Current')
plt.xlabel('Current (Amps)')
plt.ylabel('Temperature (K)')
fig.savefig('tempvscurrent.jpg')
