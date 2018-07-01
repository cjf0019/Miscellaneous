# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 18:31:26 2018

@author: InfiniteJest
"""

import numpy as np

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



def density(x, y, z):
    return 0.01*(2 + x**2 + y * (z + 1))

def I11int(x, y, z):
    return (y**2 + z**2) * density(x, y, z)

def I22int(x, y, z):
    return (x**2 + z**2) * density(x, y, z)

def I33int(x, y, z):
    return (x**2 + y**2) * density(x, y, z)

def I12int(x, y, z):
    return -x*y * density(x, y, z)

def I13int(x, y, z):
    return -x*z * density(x, y, z)

def I23int(x, y, z):
    return -y*z * density(x, y, z)

        
xint = [-1, 1]
yint = [-1.5, 1.5]
zint = [-2, 2]
I11 = trpl_int(I11int, xint, yint, zint, 250)
print("I11 = ", I11)
I12 = trpl_int(I12int, xint, yint, zint, 250)
print("I12 = ", I12)
I13 = trpl_int(I13int, xint, yint, zint, 250)
print("I13 = ", I13)
I22 = trpl_int(I22int, xint, yint, zint, 250)
print("I22 = ", I22)
I23 = trpl_int(I23int, xint, yint, zint, 250)
print("I23 = ", I23)
I33 = trpl_int(I33int, xint, yint, zint, 250)
print("I33 = ", I33)


###FIND THE ROOTS
#First define the eigenvalue root problem from the determinant
def eigenvalueeq(x):
    """
    Calculates the determinant of the eigenvalue matrix. Assumes Iij = Iji.
    """
    return (I11 - x) * ((I22 - x)*(I33 - x) - I23**2) - I12 * (I12 * (I33 - x) \
           - I13 * I23) + I13 * (I12 * I23 - I13 * (I22 - x))

def generatepointset(interval, steps):
    """
    Creates a list of uniformly distributed points along an axis. Interval is
    a list of two elements, the first the beginning and the second the end of 
    the interval. The stepsize is how much distance between each point.
    """
    length = float(interval[1] - interval[0])
    stepsize = float(length/steps)   #evenly distribute a set of x points
    pointset = []
    start = float(interval[0])     #initial point
    pointset.append(float(start))
    for i in range(steps):
        start = start + stepsize     #add point
        pointset.append(start)
    return np.array(pointset)
    

def scanf(f, points):
    """
    Scans the function for areas where it switches signs. Returns list of iterables,
    each representing the two points of a sign switch. 
    """
    yvalues = np.array([f(point) for point in points])    #calculates function over given points
    roothits = []
    for i in (range(len(yvalues)-1)):
        if (yvalues[i] < 0 and yvalues[i+1] > 0) or (yvalues[i] > 0 and yvalues[i+1] < 0):
            roothits.append((points[i], points[i+1]))
    return roothits


def differentsigns(x1, x2, f):
    return (f(x1) < 0 and f(x2) > 0) or (f(x1) > 0 and f(x2) < 0) #returns true if different signs

def secant(x1, x2, f, steps):
    """
    Finds the root of f between two points x1 and x2, in steps with the secant
    method.
    """
    for step in range(steps):
        x3 = (x1 + x2)/2     #calculate midpoint
        if differentsigns(x3, x1, f):
            x2 = x3   #move closer to x1!
        else:
            x1 = x3   #move close to x2!
    return x3
        
def falsepos(x1, x2, f, steps):
    """
    Finds the root with the false positive method.
    """
    for step in range(steps):
        invslope = (f(x2)-f(x1))/(x2-x1)
        x3 = x2 - f(x2)*invslope
        if differentsigns(x1, x3, f):
            x2 = x3
        else:
            x1 = x3
    return x3

def findroots(f, interval, intstep, steps, method='secant'):
    """
    Determies the roots over a defined interval (a list of 2 values for lower
    and upper bound). First calculates a set of intstep points, then determines
    areas where a sign change has occurred. Then, it applies either the secant
    or false positive ('falsepos') method at each region to find the root.
    """
    points = generatepointset(interval, steps)   #generate points over given interval
    hits = scanf(f, points)    #scan for sign changes
    roots = []
    if method == 'secant':
        for hit in hits:    #iterate through sign changes for each root
            roots.append(secant(hit[0], hit[1], f, steps))
    elif method == 'falsepos':
        for hit in hits:
            roots.append(falsepos(hit[0], hit[1], f, steps))
    return roots


roots = findroots(eigenvalueeq, [0, 10], 5000, 50000)
print("The roots are: ", roots)

with open('MomentofInertiaRoots.txt', 'a') as file:
    file.write("Calculates the Inertia Tensor Eigenvalues of a Box of Density Function"+"\n" \
               + "rho = 0.01*(2 + x**2 + y * (z + 1))"+"\n")
    #Give back the tensor element integrals
    file.write("The tensor elements are:" + "\n" + "I11: " + str(I11) + "\n" + \
               "I12: " + str(I12) + "\n" + "I13: " + str(I13) + "\n" + "I22: " + \
               str(I22) + "\n" + "I23: " + str(I23) + "\n" + "I33: " + str(I33) + "\n")
    
    #Give back roots
    file.write("Using a range of " + str(0) + "to" + str(10) + ", scanning 5000 points," + "\n" \
               + " performing 50000 steps of the secand method, the roots are: "+ "\n"
                + str(roots))