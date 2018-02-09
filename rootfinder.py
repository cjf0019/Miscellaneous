# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 18:31:26 2018

@author: InfiniteJest
"""

import numpy as np

def f(x):
    return np.log(x) + np.exp(x)

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
    yvalues = np.array([f(point) for point in points])
    roothits = []
    for i in (range(len(yvalues)-1)):
        if (yvalues[i] < 0 and yvalues[i+1] > 0) or (yvalues[i] > 0 and yvalues[i+1] < 0):
            roothits.append((points[i], points[i+1]))
    return roothits

def differentsigns(x1, x2, f):
    return (f(x1) < 0 and f(x2) > 0) or (f(x1) > 0 and f(x2) < 0)

def secant(x1, x2, f, steps):
    for step in range(steps):
        x3 = (x1 + x2)/2
        if differentsigns(x3, x1, f):
            x2 = x3   #move closer to x1!
        else:
            x1 = x3 
    return x3
        
def falsepos(x1, x2, f, steps):
    for step in range(steps):
        invslope = (f(x2)-f(x1))/(x2-x1)
        x3 = x2 - f(x2)*invslope
        if differentsigns(x1, x3, f):
            x2 = x3
        else:
            x1 = x3
    return x3

def findroots(f, interval, intstep, steps, method='secant'):
    points = generatepointset(interval, steps)
    hits = scanf(f, points)
    roots = []
    if method == 'secant':
        for hit in hits:
            roots.append(secant(hit[0], hit[1], f, steps))
    elif method == 'falsepos':
        for hit in hits:
            roots.append(falsepos(hit[0], hit[1], f, steps))
    return roots

