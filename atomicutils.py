# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 14:55:44 2018

@author: InfiniteJest
"""

def transnumber(numterms, termnum):
    transnum = 0
    for i in range(termnum):
        numterms -= 1
        transnum += numterms
    return transnum

