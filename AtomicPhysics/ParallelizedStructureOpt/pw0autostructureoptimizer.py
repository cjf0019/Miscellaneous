# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 15:14:12 2018

@author: InfiniteJest
"""
import os
import re
import subprocess
import shutil
import comparenistenergies, compareavalues
from scipy.optimize import minimize, basinhopping, brute
import numpy as np
import multiprocessing



das = 'das'
orbs = [17, 18]
bounds = [(0.95, 1.1) for i in orbs]
workers=10
nistfile = 'w1nistenergies.csv'
nistavalues = {'4 7 1 7 6 7 0 3': 12400000, '4 7 1 81 6 7 0 3': 108000000, \
        '2 5 1 7 0 5 2 1': 3200000, '4 5 3 7 2 5 2 1': 2140000}



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
        lambdas = ''
        line = ''
        while '&END' not in line:
            line = file.readline()
            lambdas = lambdas + line
        self.lambdas = re.sub('&END|&SRADWIN', '', lambdas)
        self.sradwin = ' &SRADWIN &END'
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
        file = open(newdas, 'w')
        file.write(self.firstline)
        file.write(self.inputs)
        file.write(self.orbs)
        file.write(self.configs)
        file.write(self.sminim)
        file.write(self.lambdas)
        file.write(self.sradwin)
        file.close
        return



def confautofxn(newvalues):
    pid = str(os.getpid())
    path = "./"+pid+"/"
    print("FIRSTPID",str(os.getpid()))
    print("PID",pid)
    if (not os.path.isdir("./"+pid)) and (os.getcwd().split('/')[-1] != pid):
        os.mkdir(pid)
        shutil.copy('./asdeck.x',path[:-1])
        shutil.copy('./das',path[:-1])
        os.chdir(pid)

    newdas = DAS(das)
#    for i, j in zip(orbs, list(newvalues)):
    if len(np.shape(newvalues)) == 0:
        newvalues = [float(newvalues)]
    for i, j in zip(orbs, newvalues):
        newdas.changelambda(i, j)
    newdas.writedas("das")
    print("CWD",os.getcwd())

    os.system("./asdeck.x <"+" "+das)    #run Autostructure
#    subprocess.run(['./asdeck.x < das'])
    print('AUTO RUN COMPLETED')
#    if trmorlvl == 'term':
#        energies = open('TERMS')
#    else:
    energies = comparenistenergies.LEVELS('LEVELS').levels
    print("ENERGIES",energies)
    nistenergies = comparenistenergies.LEVELS("../"+nistfile,sep=',').levels

    return comparenistenergies.get_sum_diff(energies,nistenergies,purity_thresh=65)



def avalautofxn(newvalues):
    pid = str(os.getpid())
    path = "./"+pid+"/"
    print("FIRSTPID",str(os.getpid()))
    print("PID",pid)
    if (not os.path.isdir("./"+pid)) and (os.getcwd().split('/')[-1] != pid):
        os.mkdir(pid)
        shutil.copy('./asdeck.x',path[:-1])
        shutil.copy('./das',path[:-1])
        os.chdir(pid)

    newdas = DAS(das)
#    for i, j in zip(orbs, list(newvalues)):
    if len(np.shape(newvalues)) == 0:
        newvalues = [float(newvalues)]
    for i, j in zip(orbs, newvalues):
        newdas.changelambda(i, j)
    newdas.writedas("das")
    print("CWD",os.getcwd())

    os.system("./asdeck.x <"+" "+das)    #run Autostructure
#    subprocess.run(['./asdeck.x < das'])
    print('AUTO RUN COMPLETED')

    diff = compareavalues.compareavals(nistavalues)
    return diff


def autofxn(das, orbs, newvalues, nistenergies, trmorlvl = 'term'):
    """
    Inserted into the optimization algorithm to produce an overall scalar cost,
    based on the differences between the calculated terms and a list of NIST energies.
    das is a das file to be used
    orbs is a list of orbital numbers to be modified, with newvalues a list of updated lambdas
    nistenergies is a list of the NIST energies corresponding to the first n terms/levels.
    """
    pid = str(os.getpid())
    path = "./"+pid+"/"
    print("PID",pid)
    if not os.path.isdir("./"+pid):
        os.mkdir(pid)
        shutil.copy('./asdeck.x',path[:-1])
        shutil.copy('./das',path[:-1])
    newdas = DAS(path+das)
#    for i, j in zip(orbs, list(newvalues)):
    if len(np.shape(newvalues)) == 0:
        newvalues = [float(newvalues)]
    for i, j in zip(orbs, newvalues):
        newdas.changelambda(i, j)
    newdas.writedas(path+"das")
    os.system(path+"asdeck.x <"+" "+path+das)    #run Autostructure
#    subprocess.run(['./asdeck.x < das'])
    print('AUTO RUN COMPLETED')
    if trmorlvl == 'term':
        terms = open(path+'TERMS')
    else:
        terms = open(path+'LEVELS')
    terms.readline()  #Read the header line
    energydiff = []
    nistcheck = 0
    while nistcheck != len(nistenergies):    #pick out NIST energies from TERMS file
#    for i in range(len(nistenergies)):
        termline = terms.readline().split()
        print(termline)
        termconf = ' '.join([termline[i] for i in (0,1,3)])
        print(termconf, 'TERMCONF')
        print(nistenergies.keys(), 'NIST')
        if termconf in nistenergies.keys():
            energydiff.append(abs(float(termline[-1]) - \
            nistenergies[termconf]))
            nistcheck += 1
    print('COMPARISON MADE')
    return sum([i**2 for i in energydiff])

   



def optimizelambdas(das, orbs, bounds = None, trmorlvl = 'term', \
			workers=1):
    initialdas = DAS(das)
    print(initialdas.lambdas)
    initlambdas = np.array([float(initialdas.lambdas.split()[i-1]) for i in orbs])
#    res = minimize(specificautofxn, initlambdas, method='nelder-mead', bounds=bounds, options={'xtol': 1e-6, 'fatol': 1e-6,  'maxiter':30})    
#    bounded_step = RandomDisplacementBounds(np.array([b[0] for b in bounds]), \
#					    np.array([b[1] for b in bounds]))
#    res = basinhopping(specificautofxn, initlambdas, minimizer_kwargs={"method":'nelder-mead'}, \
#                       niter=20, take_step=bounded_step)
#    res = basinhopping(specificautofxn, initlambdas, minimizer_kwargs={"method":'L-BFGS-B', \
#			'bounds':bounds}, niter=40, disp=True, stepsize=0.025, T=0.002) 
#    res = shgo(specificautofxn, bounds, minimizer_kwargs={"method":'L-BFGS-B'}) 
    res = brute(avalautofxn, Ns=10, ranges=bounds, disp=True, \
		workers=workers, full_output=True, finish=None) 
    return res



if __name__ == "__main__":
    print(bounds)
#nistenergies = {'1 0 1': 0.000, '3 1 3': 0.2002979, '1 1 3': 0.3878840, '3 0 4': 0.4745964, \
#        '1 2 2': 0.5183508}

#optimized = optimizelambdas(dasfile, orbnums, 'w2nistenergies.csv', bounds=bounds)
#optimized = optimizelambdas('das', [2,3], {'1 0 1': 0.000, '3 1 2':0.2002979})

    def returnopt():
        return optimizelambdas(das,orbs,bounds=bounds,trmorlvl='level', \
                               workers=workers)



    optimized = returnopt()
    print("OPT",optimized[0])

    np.save('bruteopt.npy', optimized[0], allow_pickle=False)
    np.save('bruteoptcost.npy',optimized[1], allow_pickle=False)
    np.save('brutegrid.npy',optimized[2], allow_pickle=False)
    np.save('brutecost.npy',optimized[3], allow_pickle=False)

#print("GRID",optimized.grid)
#print("JOUT",optimized.Jout)

#print("FINAL VALUE",optimized)

