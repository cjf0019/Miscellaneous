# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 11:12:14 2018

@author: InfiniteJest
"""

import os
import re
import collections
os.chdir("C:\\Users\\InfiniteJest\\Documents\\Physics\\Neon")
import numpy as np

pectext = open('pec.passfeb22')

numtrans = int(pectext.readline().split()[0])
pectext.readline()
pectext.readline()
pectext.readline()
pectext.readline()
pectext.readline()
pectext.readline()

def read_one_pec(pectext, numdens, numtemp):
    densities = [float(i) for i in pectext.readline().strip('\n').split()]
    while len(densities) < int(numdens):
        densities = densities + [float(i) for i in pectext.readline().strip('\n').split()]
    densities = np.array(densities)
    
    temps = [float(i) for i in pectext.readline().strip('\n').split()]
    while len(temps) < int(numtemp):
        temps = temps + [float(i) for i in pectext.readline().strip('\n').split()]
    temps = np.array(temps)
    matrix = []
    for i in range(int(numdens)):
        pecs = [float(j) for j in pectext.readline().strip('\n').split()]
        while len(pecs) < int(numtemp):
            pecs = pecs + [float(j) for j in pectext.readline().strip('\n').split()]
        matrix.append(pecs)
    return temps, densities, matrix

def collect_pecs(pectext, numwavelengths, remove_empties=False):
    wavelengths = []
    pecs = []
    for i in range(numwavelengths-1):
        info = pectext.readline()
        wavelength, numdens, numtemp = info.split()[0:3]
        temps, densities, pec = read_one_pec(pectext, numdens, numtemp)
        if remove_empties == True:
            if np.any(np.array(pec[:][:]) == 1e-74):
                pass
            else:
                pecs.append(pec)
                wavelengths.append(wavelength)
        else:
            pecs.append(pec)
            wavelengths.append(wavelength)
    return temps, densities, wavelengths, pecs
            
temps, densities, wavelengths, pecs = collect_pecs(pectext, 100, remove_empties=True)
wavetopec = dict(list(zip(wavelengths, pecs)))

def calculate_line_ratio(wavelength1, wavelength2, pecdict, temps, densities):
    """
    Calculates all possible Te and Ne Line ratios between two wavelengths and 
    selects the highest sloping and lowest sloping of each. Returns tuples of
    (temax, temin) and (nemax, nemin).
    """
    ratiomatrix = np.divide(pecdict[wavelength1], pecdict[wavelength2])
    telineratios = np.apply_along_axis(lambda row: np.polyfit(temps, row, 1)[0], 1, ratiomatrix)
    telineratios = (np.max(telineratios), np.min(telineratios))
    nelineratios = np.apply_along_axis(lambda column: np.polyfit(densities, column, 1)[0], 0, ratiomatrix)
    nelineratios = (np.max(nelineratios), np.min(nelineratios))
    return telineratios, nelineratios


def find_best_and_worst_ratio(ratiomat):
    