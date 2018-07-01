# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:39:39 2018

@author: InfiniteJest
"""

import pandas as pd 
import os
import re
import collections
os.chdir("C:\\Users\\InfiniteJest\\Documents\\Physics\\Neon")

###The following indices represent a change in the infinite omega file to accomodate NIST shifting
#newindices = [i + 1 for i in range(27)] + [29, 30, 31, 32, 28, 33, 34, 35, 36, \
##             37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 53, 58, 59, \
##             51, 52, 54, 55, 56, 57, 60, 61, 62, 63, 64, 65, 66, 68, 70, 69, 71, 72, \
 #            73, 74, 75, 67, 76, 77, 78, 79] + [i for i in range(80, 634)]
#oldindices = [i + 1 for i in range(633)]
#oldtonew = dict(zip(oldindices, newindices))
newindices = [i + 1 for i in range(57)] + [59, 58] + [i for i in range(60, 73)] + \
            [74, 73] + [i for i in range(75, 80)]
oldindices = [i + 1 for i in range(79)]
oldtonew = dict(zip(oldindices, newindices))


adf04text = open('adf04apr17')

firstline = adf04text.readline()
charge = firstline.split()[1]
atomnum = firstline.split()[2]
ionpot = firstline.split()[4]


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
    numtolvl = collections.OrderedDict()  #create dictionary for number to lvl
    for level in splitit2:
        numandlvl = level.split(' ', 1)   #split the index from the level
        numtolvl[numandlvl[0]] = numandlvl[1]   #index as key, level as value
    return numtolvl

def reorderenergies(oldtonew, numtolvl):
    new = collections.OrderedDict()
    for i in numtolvl:
        if len(i) == 1:
            new['0'+i] = numtolvl[str(oldtonew[int(i)])]
        else:
            new[i] = numtolvl[str(oldtonew[int(i)])]
            print(new[i])
    return new

numtolvl = extract_level_indices(energies)
lvltonum = collections.OrderedDict(zip(numtolvl.values(),numtolvl.keys()))

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
    numtorate = collections.OrderedDict()
    for trans in transitions:
        splittrans = trans.split()
        if (len(splittrans[0]) + len(splittrans[1]) == 2) or \
            (len(splittrans[0]) + len(splittrans[1]) == 3):
            numtorate['   '.join(splittrans[:2])] = ' '.join(splittrans[2:])
        elif (len(splittrans[0]) + len(splittrans[1]) == 4) or \
         (len(splittrans[0]) + len(splittrans[1]) == 5):
            numtorate['  '.join(splittrans[:2])] = ' '.join(splittrans[2:])
        else:
            numtorate[' '.join(splittrans[:2])] = ' '.join(splittrans[2:])
#        numandrate = re.split('(?<=[0-9]) (?=[0-9])', trans, maxsplit=1)
#        numtorate[numandrate[0]] = numandrate[1]
    return numtorate

numtorate = extract_trans_indices(rates)


def convertindices(oldtonew, numtorate):
    """
    reorders the adf04 rates to correspond to a new set of indices. 'oldtonew' is
    a dictionary in which oldtonew[old] = new . 
    """
    newnumtorate = collections.OrderedDict()
    for i in numtorate:
        rate = numtorate[i]
        indices = i.split()
        indices = [str(oldtonew[int(i)]) for i in indices]
        if (len(indices[1]) + len(indices[0]) == 2) or \
            (len(indices[1]) + len(indices[0]) == 3):
            newnumtorate['   '.join(indices)] = rate
        elif (len(indices[1]) + len(indices[0]) == 4) or \
         (len(indices[1]) + len(indices[0]) == 5):
            newnumtorate['  '.join(indices)] = rate
        else:
            newnumtorate[' '.join(indices)] = rate
    return newnumtorate
            
newnumtorate = convertindices(oldtonew, numtorate)


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


def replace_inf_point(adf041dict, adf042dict, onlyempty=False):
    #will write out the column of adf042 into adf041
    #use the "numtorate" for both adf04 files as the dicts
    inf_points = []
    if onlyempty:
        for key in adf041dict.keys():
            if key in adf042dict.keys():    #might be missing a transition in second adf04 file
                if adf041dict[key].split(' ')[-1][-7:] == '0.00+00':
                    print('yes')
                    inf_points.append(adf042dict[key].split(' ')[-1][-8:])
                else:
                    inf_points.append(adf041dict[key].split(' ')[-1][-8:]) #keep original values when not zero
            else:
                inf_points.append(adf041dict[key].split(' ')[-1][-8:])
    else:
        for key in adf041dict.keys():
            if key in adf042dict.keys():
                inf_points.append(adf042dict[key].split(' ')[-1][-8:])
            else: 
                inf_points.append(adf041dict[key].split(' ')[-1][-8:])
    return inf_points


def write_column_to_file(collist, file):
    with open(file, 'w') as f:
        f.write('\n'.join(collist))
    return



#############################################################################
#Will also need to reorder the adasexj.in file generated during inf omega calculation.
#Then use in adasexj.x.

adasexjin = open('adasexj.in')

def adasexjlvlreorder(oldtonew, adasexjin):
    adasexjin.readline()
    line = ''
    energies = ''
    while not re.search("NAME", line.strip()):   #add all the energies
        line = adasexjin.readline()
        energies = energies + line
    nmtolvl = extract_level_indices(energies)   #extract index to level dictionary
    new = reorderenergies(oldtonew, nmtolvl)    #change energy indices (doesn't change list order yet)
    ind = [str(i+1) for i in range(79)]
    ind2 = []
    for i in ind:
        if len(i) == 1:
            ind2.append('0' + i)   #adds zero to single digits, as done in reorder energies
        else:
            ind2.append(i)
    new2 = []
    for i in ind2:
        new2.append(new[i])
    return new2
        

def addNISTAvalues(adf04dict, nistfile, energies):
    #use the "numtorate" for both adf04 files as the dicts
    A_coeff = []
    energies = energies.split('\n')[1:]
    adf04energies = []
    for i in energies:
        adf04energies.append(i.split()[-1])
        
    nist = open(nistfile)    #Assumes a csv file from the NIST webpage, search of all A values
    nistavalues = []
    nistlowerenergy = []
    nisthigherenergy = []
    nist.readline()
    nist.readline()
    nist.readline()
    for line in nist:
        if re.search('[0-9]', line):
            split = line.strip('\n').split(',')
            nistavalues.append(split[3].replace('E',''))
            nistenergies = re.split('[\s]+\-[\s]+', split[4])
            nistlowerenergy.append(nistenergies[1])
            nisthigherenergy.append(nistenergies[0])
    transitions = []
    for i in range(len(nistlowerenergy)):
        if nistlowerenergy[i] in adf04energies:
            if nisthigherenergy[i] in adf04energies:
                transitions.append([str(adf04energies.index(nistlowerenergy[i])+1), \
                                  str(adf04energies.index(nisthigherenergy[i])+1)])
            else:
                "HIGHER ENERGY MISSING!!!"
        else:
            "LOWER ENERGY MISSING!!!"
    print(transitions)
    
    for key in adf04dict.keys():
        if key.split() in transitions:
            print(nistavalues[transitions.index(key.split())])
            A_coeff.append(nistavalues[transitions.index(key.split())])
        else:
            A_coeff.append(adf04dict[key].split(' ')[0])
    return A_coeff


def getlvl2termnums(energies):
    """
    generates a list whose elements correspond to the term number of a level in the
    adf04 file. So, the 0th index would correspond to the term number (starting at 1)
    of the first level.
    """
    energies = energies.split('\n')[1:]
    terms = {}
    termnums = [0]
    for energy in energies:
        term = energy.split()[-3].strip('(')
        term = energy.split()[-5] + energy.split()[-4] + term
        if term not in terms.keys():
            terms[term] = max(termnums)+1
            termnums.append(max(termnums)+1)
        else:
            termnums.append(terms[term])
    termnums.pop(0)
    return termnums
        
def write_adf04(file, numtorate):
    """
    NOTE: Currently takes the global firstline, energies and tempsline used...
    In the future should put these into a class, and also add in updated energies.
    Right now, will still need to manually switch order in energy list if changing
    them around.
    """
    newadf04 = open(file, 'a')
    newadf04.write(firstline)
    newadf04.write(energies[2:])
    newadf04.write('\n   -1\n')
    newadf04.write(tempsline+'\n')
    for i in numtorate:
        if len(i) == 5:
            newadf04.write('   '+i+' '+numtorate[i]+'\n')
        else:
            newadf04.write('  '+i+' '+numtorate[i]+'\n')
    newadf04.write('  -1\n')
    newadf04.write('  -1  -1\n')
    newadf04.close()
    return

    
    


#### TO GENERATE NEW REORDERED ADF04 FILE:
newindices = [i + 1 for i in range(57)] + [59, 58] + [i for i in range(60, 73)] + \
            [74, 73] + [i for i in range(75, 80)]
oldindices = [i + 1 for i in range(79)]
oldtonew = dict(zip(oldindices, newindices))

