"""
Compares the quantum numbers of terms in two OMEGA files. User inputs the two file names.
"""


import numpy as np
from sys import argv

def omgterms(omega):
    file = open(omega)

    file.readline()
    file.readline()
    line = file.readline().split()
    terms = []

    for i in range(15960):
        one = int(line.pop(0))
        two = int(line.pop(0))
        terms += [[abs(one), two]]
    return terms

def compare(omgterms1,omgterms2):
    if len(omgterms1) != len(omgterms2):
        print("TERMS NOT EQUAL LENGTHS")
    else:
        for i in range(len(omgterms1)):
            if (omgterms1[i][0]!=omgterms2[i][0]) or (omgterms1[i][1]!=omgterms2[i][1]):
                print(i, '\t', omgterms1[i], '\t', omgterms2[i])
    print(omgterms1[0][0])


compare(omgterms(argv[1]), omgterms(argv[2]))

