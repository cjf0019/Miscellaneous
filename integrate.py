# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 12:33:17 2018
Calculates the Moment of Inertia of an ellipsoid via the integral of r^2 * dm,
where r is the distance from the axis of rotation and dm the mass differential.
Here, we will assume the azimuthal to be integrated out (to 2pi). So, we are
really integrating the area element of 1 = (x^2/a^2 + z^2/c^2) times 2pi for the
volume. The distance is x^2 from the axis in this case. We can also see that 
z will be values from +sqrt(c^2 * (1- (x^2/a^2))) to negative of that value, or
twice the positive of the value with zero as the lower bound. So, we end up 
getting (M / V) *  4pi * integral[from 0 to a](x^2 * sqrt( c^2 (1 - x^2/a^2) dx)
@author: Connor Favreau
"""
import numpy as np
import os
os.chdir('C:\\Users\\InfiniteJest\\Documents\\Python_Scripts')

def ellipsoidvolume(x, y, a, b, c):
    f = 8 * float(c) * np.sqrt(1-float(x)**2/float(a)**2-float(y)**2/float(b)**2)    #the 8 comes from the axis symmetry for all variables (2*2*2)
    return f

def ellipsoidx2(x, y, a, b, c):
    f = (x**2) * ellipsoidvolume(x, y, a, b, c)
    return f


#for a sphere:
def spherevol(x, y):
    """
    a = c = 1
    """
    return ellipsoidvolume(x, y, 1, 1, 1)

def spherex2(x, y):
    return ellipsoidx2(x, y, 1, 1, 1)

def ellipsvolinstance(x, y):
    f = ellipsoidvolume(x, y, 1, 1, 2)
    return f

def ellipsx2instnace(x, y):
    return ellipsoidx2(x, y, 1, 1, 2)



def dbl_int(f, intx, inty, step):
    """
    Calculate Double Integrals. Intx and inty are the x and y ranges to integrate
    over. f is the function and step the number of slices in each dimension.
    Implements the midpoint rule.
    """
    xrange = abs(intx[1] - intx[0])
    yrange = abs(inty[1] - inty[0])
    hx = xrange/float(step)    #change per step in each direction
    hy = yrange/float(step)
    integral = 0
    for i in range(step):
        x1 = float(intx[0] + (i  + 0.5) * hx)  #midpoint of x
        for j in range(step):
            y1 = float(inty[0] + (j + 0.5) * hy)   #midpoint of y
            func = f(x1, y1)
            if np.isnan(func):    #get rid of invalid points
                pass
            else:
                integral += func * hx * hy   #multiple value at midpoint by dx and dy
    return integral



with open('Ellipsoid.txt', 'a') as file:
    file.write("Moment of Inertia of a Sphere and Ellipsoid of Mass 0.01 kg, uniform density"+"\n")
    #Calculate sphere volume and moi
    spherevolume = dbl_int(spherevol, [0, 1], [0, 1], 1500)
    spheremoi = (0.01/spherevolume) * dbl_int(spherex2, [0, 1], [0, 1], 1500) #first term is the density
    file.write("The sphere volume is: "+ repr(spherevolume) + ". Done in 1500 Steps."+"\n")
    file.write("The sphere moment of inertia is: " + repr(spheremoi) + ". Done in 1500 Steps."+"\n")
    
    #Convergence Check
    concheck = (0.01/spherevolume) * dbl_int(spherex2, [0, 1], [0, 1], 1300)
    file.write("Convergence check: Step at 1300 within"+ repr(100*abs(spheremoi-concheck)/spheremoi)+"percent"+"\n")
    
    #Calculate ellipsoid volume and moi
    ellvol = dbl_int(ellipsvolinstance, [0, 1], [0, 1], 1500)
    ellmoi = (0.01/ellvol) * dbl_int(ellipsx2instnace, [0, 1], [0, 1], 1500)
    file.write("The Ellipsoid volume is: "+ repr(ellvol)+ ". Done in 1500 Steps."+"\n")
    file.write("The moment of inertia is: "+repr(ellmoi) +". Done in 1500 Steps."+"\n")







#############################################################################
### The following code wasn't used for the final moment of inertia, but could
### could have been if I had done the triple integral or simpsons rule.

def simpsons_rule(f, interval, step):
    """
    Simpson's Rule Integration for 1D. f is the function, interval the integration range,
    and step the number of slices to add.
    """
    integral = 0
    begin = float(interval[0])   #integration start point
    end = float(interval[1])    #integration end point
    h = abs(end - begin)/float(step)    # change in x at each step
    for i in range(step):   
        x1 = float(begin + i * h)
        x2 = float(begin + (i + 1) * h)
        integral += (h/6)*(f(x1) + 4*f(0.5*(x1+x2)) + f(x2))   #SImpson's Rule
    return integral

def generate_mesh(intx, inty, step):
    """
    Not implemented, but can be used to iterate through a set of mesh points
    in a square grid, i.e., (x1, y1), (x2, y1), (x1, y2), (x2, y2), etc...
    """
    xrange = abs(intx[1] - intx[0])    #range of x values
    yrange = abs(inty[1] - inty[0])    #range of y values
    hx = xrange/float(step)     #change in x at each step
    hy = yrange/float(step)     #change in y
    xpoints = []
    ypoints = []
    for i in range(step):
        x1 = float(intx[0] + (i  + 0.5) * hx)   #set up x points, at midpoint
        y1 = float(inty[0] + (i + 0.5) * hy)    #set up y points, at midpoint
        xpoints.append(x1)
        ypoints.append(y1)
    for i in xpoints:
        for j in ypoints:
            yield (i, j)   #acts as a generator for memory efficiency, yields each point

def dbl_int(f, intx, inty, step):
    """
    Calculate Double Integrals. Intx and inty are the x and y ranges to integrate
    over. f is the function and step the number of slices in each dimension.
    Implements the midpoint rule.
    """
    xrange = abs(intx[1] - intx[0])
    yrange = abs(inty[1] - inty[0])
    hx = xrange/float(step)
    hy = yrange/float(step)
    integral = 0
    for i in range(step):
        x1 = float(intx[0] + (i  + 0.5) * hx)  #midpoint of x
        for j in range(step):
            y1 = float(inty[0] + (j + 0.5) * hy)   #midpoint of y
            func = f(x1, y1)
            if np.isnan(func):    #get rid of invalid points
                pass
            else:
                integral += func * hx * hy   #multiple value at midpoint by dx and dy
    return integral


def trpl_int(f, intx, inty, intz, step):
    """
    Calculates triple integrals.
    """
    xrange = abs(intx[1] - intx[0])
    yrange = abs(inty[1] - inty[0])
    zrange = abs(intz[1] - intz[0])
    hx = xrange/float(step)
    hy = yrange/float(step)
    hz = zrange/float(step)
    integral = 0
    for i in range(step):
        x1 = float(intx[0] + (i  + 0.5) * hx)  #midpoint of x
        for j in range(step):
            y1 = float(inty[0] + (j + 0.5) * hy)   #midpoint of y
            for k in range(step):
                z1 = float(intz[0] + (k + 0.5) * hz)
                integral += f(x1, y1, z1) * hx * hy * hz  #multiple value at midpoint by dx and dy
    return integral




