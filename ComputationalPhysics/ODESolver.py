# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 15:11:46 2018

@author: InfiniteJest
"""

def v(t, x, xprime):    #the dx/dt function
    return xprime

def f(t, x, xprime):    #the dx2/dt2 function
    return -4.0*x -2.0*xprime

def RK4_onestep(fxn, tinterval, x, xprime, xorv='x'):
    """
    Calculates the next time step of a differential equation using the fourth order
    Runge Kutte method. fxn is the differential equation, tinterval a 2-ple list
    that includes the initial time and final time, xprime the initial value of the
    first derivative, and x the initial value of the solution. xorv switches between
    solving the equation in terms of the "x" or its derivative "v"
    """
    initt = tinterval[0]    #initial time step
    newt = tinterval[1]     #time step to find solution at
    deltat = abs(newt - initt)    #range of t values
    #the Runge Kutte Fourth order method
    kone = fxn(initt, x, xprime)
    ktwo = fxn((initt+deltat/2), (x+deltat*kone/2), (xprime+deltat*kone/2))
    kthree = fxn((initt+deltat/2), (x+deltat*ktwo/2), (xprime+deltat*ktwo/2))
    kfour = fxn((initt+deltat), (x+deltat*kthree), (xprime+deltat*kthree))
    if xorv == 'x':     #for integrating dx/dt
        x += (deltat/6)*(kone + 2*ktwo + 2*kthree + kfour)
        return x
    elif xorv == 'v':     #for integrating dv/dt
        xprime += (deltat/6)*(kone + 2*ktwo + 2*kthree + kfour)
        return xprime
    else:
        print("PROBLEM! Set 'xorv' to 'x' or to its derivative 'v' for integration.")

def SecondOrderRK4(f, v, trange, initx, initxprime, step):
    """
    Performs RK4_onesteps at a specific trange of times, given number of steps 'step',
    initial x 'initx', initial x prime 'initxprime' (the velocity), a second order
    differential equation f, and a first order v.
    """
    initt = trange[0]    #initial t over whole range
    finalt = trange[1]   #final t over whole range
    deltat = abs(finalt-initt)/float(step)  #change in time at each step
    x = initx
    t = initt
    xpoints = [x]
    xprime = initxprime
    xprimepoints = [xprime]
    tpoints = [t]
    timesteps = [i+1 for i in range(step)]   #list of time steps
    for timestep in timesteps:      
        #first do the dv/dt integration
        xprime = RK4_onestep(f, [t, (initt + deltat*timestep)], x, xprime, xorv='v')
        #the dx/dt integration (using the velocity function)
        x = RK4_onestep(v, [t, (initt + deltat*timestep)], x, xprimepoints[-1], xorv='x')             
        xpoints.append(x)
        xprimepoints.append(xprime)
        t = float(initt + timestep * deltat)
        tpoints.append(t)
    return xpoints, xprimepoints, tpoints

xpoints1, xprimepoints1, tpoints1 = SecondOrderRK4(f, v, [0,20], 2, 30, 20)
xpoints2, xprimepoints2, tpoints2 = SecondOrderRK4(f, v, [0,20], 2, 30, 100)
xpoints3, xprimepoints3, tpoints3 = SecondOrderRK4(f, v, [0,20], 2, 30, 1000)
xpoints4, xprimepoints4, tpoints4 = SecondOrderRK4(f, v, [0,20], 2, 30, 2000)
xpoints5, xprimepoints5, tpoints5 = SecondOrderRK4(f, v, [0,20], 2, 30, 5000)
xpoints6, xprimepoints6, tpoints6 = SecondOrderRK4(f, v, [0,20], 2, 30, 10000)

import numpy as np
xpoints = np.expand_dims(np.array([xpoints1[-1]]+[xpoints2[-1]]+[xpoints3[-1]] \
                        +[xpoints4[-1]]+[xpoints5[-1]]+[xpoints6[-1]]), axis=1)

deltat = np.expand_dims(np.array([tpoints1[1]-tpoints1[0]]+[tpoints2[1]-tpoints2[0]]+ \
                        [tpoints3[1]-tpoints3[0]]+[tpoints4[1]-tpoints4[0]]+ \
                        [tpoints5[1]-tpoints5[0]]+[tpoints6[1]-tpoints6[0]]), axis=1)


#Write the raw data to a file
import pandas as pd
data = np.concatenate((deltat, xpoints), axis=1)
df = pd.DataFrame(data)
df.columns = ['Delta t', 'Final x']

file = open('Module4HW2', 'a')
file.write(str(df))
file.write('\n')
file.write('Distance over Time\n')
data = np.concatenate((np.expand_dims(np.array(tpoints6), axis=1), np.expand_dims(np.array(xpoints6), axis=1)), axis=1)
df = pd.DataFrame(data)
df.round(decimals=6)
df.columns = ['Time', 'Position']
with pd.option_context('display.max_rows', None, 'display.max_columns', 3):
    file.write(str(df))
file.close()

import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(tpoints6, xpoints6)
fig.suptitle('Damped Spring Position over Time')
plt.xlabel('Time (s)')
plt.ylabel('Distance (m)')
fig.savefig('DampedSpring.jpg')