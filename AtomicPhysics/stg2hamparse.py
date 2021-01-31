# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:29:55 2019

@author: InfiniteJest
"""

import re
import os
import numpy as np


os.chdir('C:\\Users\\InfiniteJest\\Documents\\Python_Scripts\\Physics')

rout2r = open('rout2r002')

def parse1cc(rout2r, nrang2):
    ccblock = np.empty((nrang2,nrang2))
    mnp1chunks = nrang2 // 8
    line = ''
    while "CONTINUUM-CONTINUUM CONTRIBUTION" not in line:
        line = rout2r.readline()
    print(line)
    channels = re.search('[0-9]+',re.search('CHANNEL.*[0-9]+',line).group(0))
    print(channels.group(0))
    ch1 = int(channels.group(0))
#    ch2 = int(channels.group(1))   GETTING ISSUE WITH RE AND THIS
    rout2r.readline()
    indcol = 0
    indrow = 0
    while mnp1chunks != -1:
        mnp1chunks -= 1
        for i in range(nrang2):
            ccblock[indrow][indcol:indcol+8] = [np.float(i) for i in rout2r.readline().split()]
            indrow += 1
        indrow = 0
        indcol += 8
        rout2r.readline()
        rout2r.readline()
        rout2r.readline()
    return ccblock


def gettosetmx(rout2r):
    line = ''
    while "NRANG2 =" not in line:
        line = rout2r.readline()
    nrang2 = int(re.search('[0-9]+',re.search('NRANG2 =.*[0-9]+',line).group(0)).group(0))
    while "L =" not in line:
        line = rout2r.readline()
    syminfo = line
    while "ENTER CONTINUUM-CONTINUUM LOOP" not in line:
        line = rout2r.readline()


from sklearn.decomposition import PCA, TruncatedSVD
pca = PCA(n_components=50)
reduced = pca.fit_transform(test)
lsa = TruncatedSVD(n_components=50)
lsareduced = lsa.fit_transform(test)

import matplotlib.pyplot as plt

plt.imshow(reduced, cmap='hot', interpolation='nearest')
