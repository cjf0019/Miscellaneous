# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 12:30:29 2018
Reads in a TERMS file and separates pseudostate configs from nonpseudostate
configs (User must manually indicate which are not pseudostates). 
@author: InfiniteJest
"""

import os
os.chdir("C:\\Users\\InfiniteJest\\Documents\\Physics\\Tungsten")
import pandas as pd

terms = pd.read_csv('TERMSn11L7', sep='\s+')
terms = terms.drop(terms.index[[-1]])
terms.index = [i + 1 for i in terms.index]

class TERMS(object):
    def __init__(self, file, ip=None, pseudolist=None, opticallist=None):
        self.terms = pd.read_csv(file, sep='\s+')
        self.terms = self.terms.drop(self.terms.index[[-1]])
        self.terms.index = [i + 1 for i in self.terms.index]
        self.numterms = max(self.terms.index)
        if ip is not None:
            self.aboveip = self.terms[self.terms['ENERGY(RYD)'] > ip]

        if pseudolist is not None and opticallist is not None:
            self.pseudolist = pseudolist            
            self.opticallist = opticallist
            self.pseudo = self.terms.query('CF in '+str(self.pseudolist))
            self.optical = self.terms.query('CF not in '+str(self.opticallist))
            if len(pseudolist+opticallist) != self.numterms:
                print('WARNING: PSEUDOSTATE AND OPTICAL LIST DOES NOT MATCH TOTAL # TERMS.')
            missing = [i+1 for i in list(range(self.numterms)) if i+1 not in pseudolist+opticallist]
            if len(missing) > 0:
                print('THE FOLLOWING TERM #s NOT ACCOUNTED FOR: ',' '.join(missing))
            shared = [i for i in pseudolist if i in opticallist]
            if len(shared) > 0:
                raise ValueError("TERMS ", " ".join(shared), "ARE PRESENT IN BOTH OPTICAL AND PSEUDO LISTS.")
                      
        elif pseudolist is not None:
            self.pseudolist = pseudolist
            self.pseudo = self.terms.query('CF in '+str(self.pseudolist))
            self.opticallist = [i+1 for i in list(range(self.numterms)) if i+1 not in pseudolist]
            self.optical = self.terms.query('CF in '+str(self.opticallist))
        
        elif opticallist is not None:
            self.opticallist = opticallist
            self.optical = self.terms.query('CF in '+str(self.opticallist))
            self.pseudolist = [i+1 for i in list(range(self.numterms)) if i+1 not in pseudolist]
            self.pseudo = self.terms.query('CF in '+str(self.pseudolist))
        return

    def makedirionsnlist(self, outfile):
        return
    
    def makeenergyvsconfnumlist(self, termdf, outfile):
        """
        Makes a csv of the energy vs. configuration number... useful for assesssing the
        quality of pseudostates.
        """
        df = pd.concat((termdf['ENERGY(RYD)'], termdf['CF']), axis=1)
        df.columns = ['#ENERGY(RYD)', 'CF']   #make easier to plot in gnuplot
        df.to_csv(outfile, sep='\t', index=False)
        return

    def gettermenergies(self, L=None, S=None, P=None, CF=None):
        if L==None and S==None and P==None and CF==None:
            raise ValueError('MUST SPECIFY AN L, S, P, OR CONF NUMBER!!!')
        else:
            qnums = locals()
            toquery = ' and '.join(['%s in %s' % (qnum, qnums[qnum]) for qnum \
                                                  in qnums if qnums[qnum] is not None and qnum != 'self'])
            result = self.terms.query(toquery)
            result = result['ENERGY(RYD)']
        return result




ionpotterm = 433    #at what term number does the ionization start?
nonpseudoconfigs = [1, 2, 3, 4, 7, 8]

nonpseudo = terms.query('CF in '+str(nonpseudoconfigs))
aboveconfnums = [i for i in list(nonpseudo.index) if i > ionpotterm]   #term nums above ip

pseudo = terms.query('CF not in '+str(nonpseudoconfigs))
pseudoionconfnums = [i for i in list(pseudo.index) if i > ionpotterm]  #pseudo term nums above ip


#iterate over pseudoterms above the ip. 
pseudoterms = []
begin = pseudoionconfnums[0]
for i in range(len(pseudoionconfnums)-1)[1:]:
    if pseudoionconfnums[i+1] - pseudoionconfnums[i] != 1:
        print('here')
        end = pseudoionconfnums[i]
        pseudoterms = pseudoterms + list(range(begin, end+1))
        begin = pseudoionconfnums[i+1]
if pseudoionconfnums[-1] - pseudoionconfnums[-2] == 1:
    pseudoterms = pseudoterms + list(range(begin, pseudoionconfnums[-1]+1))
else:
    pseudoterms = pseudoterms + list(range(pseudoionconfnums[-1], pseudoionconfnums[-1]+1))



pseudotermchunks = []
begin = pseudoionconfnums[0]
for i in range(len(pseudoionconfnums)-1)[1:]:
    if pseudoionconfnums[i+1] - pseudoionconfnums[i] != 1:
        end = pseudoionconfnums[i]
        pseudotermchunks.append((begin, end))
        begin = pseudoionconfnums[i+1]
if pseudoionconfnums[-1] - pseudoionconfnums[-2] == 1:
    pseudotermchunks.append((begin, pseudoionconfnums[-1]))
else:
    pseudotermchunks.append((pseudoionconfnums[-1], pseudoionconfnums[-1]))

stgsiglist = open('stgsign11L7directionlist.dat', 'a')

forstgsig = [0]
forstgsig = forstgsig*(len(pseudoionconfnums) + len(aboveconfnums))
for i in pseudoterms:
    forstgsig[i-ionpotterm] = 1

for term in forstgsig:
  stgsiglist.write("%s\n" % term)


forstgsig = []
for i in range(len(pseudoterms) + len(aboveionconfnums)):
    forstgsig.append(0)

