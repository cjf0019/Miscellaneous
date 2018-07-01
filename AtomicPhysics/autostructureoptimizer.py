# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 15:14:12 2018

@author: InfiniteJest
"""
import os
import re
os.chdir("C:\\Users\\InfiniteJest\\Documents\\Python_Scripts\\Physics")
import subprocess
import shutil
from scipy.optimize import minimize
import numpy as np

class DAS:
    def __init__(self, file):
        self.file = file
        file = open(file)
        self.firstline = file.readline()
        self.inputs = file.readline()
        self.nconfs = int(re.search('[0-9]+', re.search('MXCONF(\=| \=)[0-9]*', self.inputs).group(0)).group(0))
        self.norbs = int(re.search('[0-9]+', re.search('MXVORB(\=| \=)[0-9]*', self.inputs).group(0)).group(0))
        self.orbs = file.readline()
        self.configs = file.readline()
        line = ''
        while '&SMINIM' not in line:
            self.configs = self.configs + line
            line = file.readline()
        self.sminim = line
        self.nlam = int(re.search('[0-9]+', re.search('NLAM(\=| \=)[0-9]*', self.sminim).group(0)).group(0))
        self.lambdas = file.readline()
        line = ''
        while '&END' not in line:
            self.lambdas = self.lambdas + line
            line = file.readline()
        self.sradwin = line
        return

    def changelambda(self, lambdanum, newvalue):
        """
        Lambdanum is the actual number, not the index.
        """
        lambdas = self.lambdas.split()
        lambdas[lambdanum-1] = str(newvalue)
        for j in [i + 1 for i in range(len(lambdas) // 10)]:
            lambdas.insert(10*j, '\n')  #line break every 10 values
        self.lambdas = '  '.join(lambdas)
        print('Lambda Number ', lambdanum, "replaced with value", newvalue)
        return 
    
    def writedas(self, newdas):
        file = open(newdas, 'a')
        file.write(self.firstline)
        file.write(self.inputs)
        file.write(self.orbs)
        file.write(self.configs)
        file.write(self.sminim)
        file.write(self.lambdas)
        file.write(self.sradwin)
        file.close
        return
    
das = DAS('dasd2l5')

def autofxn(das, orbs, newvalues, nistenergies, trmorlvl = 'term'):
    """
    Inserted into the Nelder Mead algorithm to produce an overall scalar cost,
    based on the differences between the calculated terms and a list of NIST energies.
    das is a das file to be used
    orbs is a list of orbital numbers to be modified, with newvalues a list of updated lambdas
    nistenergies is a list of the NIST energies corresponding to the first n terms/levels.
    """
    newdas = DAS(das)
    for i, j in zip(orbs, list(newvalues)):
        newdas.changelambda(i, j)
    newdas.writedas('das')
    subprocess.run(['./asdeck.x < das'])

    if trmorlvl == 'term':
        terms = open('TERMS')
    else:
        terms = open('LEVELS')
    terms.readline()
    energydiff = []
    for i in range(len(nistenergies)):
        energydiff.append(abs(float(terms.readline().split()[-1])-nistenergies[i]))
    return sum([i**2 for i in energydiff])

    
def optimizelambdas(das, orbs, nistenergies, trmorlvl = 'term'):
    initialdas = DAS(das)
    initlambdas = np.array([float(initialdas.lambdas.split()[i-1]) for i in orbs])
    def specificautofxn(newvalues):
        return autofxn(initialdas, orbs, newvalues, nistenergies, trmorlvl)
    res = minimize(specificautofxn, initlambdas, method='nelder-mead', options={'xtol': 1e-3, 'maxiter':30})    
    return res


    
    
    
    
    
shutil.copy2(os.getcwd()+'das', 'newdas')

        
        
        