# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:39:39 2018

@author: InfiniteJest
"""

import pandas as pd 
import os
import re
os.chdir("C:\\Users\\InfiniteJest\\Documents\\Physics\\Neon")


adf04text = open('adf04ne0')

firstline = adf04text.readline().split()
charge = firstline[1]
atomnum = firstline[2]
ionpot = firstline[4]


#Read in the energies section of the adf04 file as a block of text
line = ''
energies = ''
while line.strip() != '-1':
    line = adf04text.readline()
    energies = energies + line
    
energies = '\n' + energies    #add a newline for first energy
energies = energies.replace('\n   -1\n', '')   #get rid of that last -1

def extract_level_indices(energies):
    splitit = re.split('(?<=\.0\))[\s]', energies)   #split into list
    splitit2 = [re.sub('.*\n[\s]+(?=[0-9])', '', i) for i in splitit]  #leave only lvl number and lvl
    splitit2.pop(-1)
    numtolvl = {}  #create dictionary for number to lvl
    for level in splitit2:
        numandlvl = level.split(' ', 1)   #split the index from the level
        numtolvl[numandlvl[0]] = numandlvl[1]   #index as key, level as value
    return numtolvl

numtolvl = extract_level_indices(energies)
lvltonum = dict(zip(numtolvl.values(),numtolvl.keys()))

#Extract the temperatures
tempsline = adf04text.readline()
tempsline = tempsline.strip('\n')
temps = tempsline.split()[2:]

line = ''
rates = ''
while line.strip() != '-1':
    line = adf04text.readline()
    rates = rates + line
    
rates = '\n' + rates    #add a newline for first energy
rates = rates.replace('\n  -1\n', '')   #get rid of that last -1

def convert_to_scientific(text):
    #adds the extra 'e' to numbers that is missing in adf04 files, for formatting and so they can be used by numpy
    converted = text.replace('+', 'e+')
    converted = converted.replace('-', 'e-')
    return converted

def extract_trans_indices(rates):
    ratetext = re.sub('(?<=[\n])[\s]+(?=[0-9])', '', rates) #remove extra spacing before indices
    transitions = ratetext.split('\n')[1:]
    numtorate = {}
    for trans in transitions:
        numandrate = re.split('(?<=[0-9]) (?=[0-9])', trans, maxsplit=1)
        numtorate[numandrate[0]] = numandrate[1]
#        numandrate = trans.split(' ', 2)
#        numtorate[' '.join(numandrate[:2])] = numandrate[2]
    return numtorate

numtorate = extract_trans_indices(rates)
ratetonum = dict(zip(numtorate.values(),numtorate.keys()))

###COULD ADD CODE HERE TO SWITCH BETWEEN INDICES OF DIFFERENT ADF04 FILES


def replace_A_coeff(adf041dict, adf042dict, onlyempty=False):
    #will write out the column of adf042 into adf041
    #use the "numtorate" for both adf04 files as the dicts
    A_coeff = []
    if onlyempty:
        for key in adf041dict.keys():
            if key in adf042dict.keys():    #might be missing a transition in second adf04 file
                if adf041dict[key].split(' ')[0] == '1.00-30':
                    A_coeff.append(adf042dict[key].split(' ')[0])
                else:
                    A_coeff.append(adf041dict[key].split(' ')[0]) #keep original values when not zero
            else:
                A_coeff.append(adf041dict[key].split(' ')[0])
    else:
        for key in adf041dict.keys():
            if key in adf042dict.keys():
                A_coeff.append(adf042dict[key].split(' ')[0])
            else: 
                A_coeff.append(adf041dict[key].split(' ')[0])
    return A_coeff

def write_column_to_file(collist, file):
    with open(file, 'w') as f:
        f.write('\n'.join(collist))
    return

