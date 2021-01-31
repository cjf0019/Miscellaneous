# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 12:30:29 2018
Reads in a TERMS file and separates pseudostate configs from nonpseudostate
configs (User must manually indicate which are not pseudostates). 
@author: InfiniteJest
"""

import os
import pandas as pd
import sys

#levels = pd.read_csv(sys.argv[1], sep='\s+')
#levels = pd.read_csv('LEVELS', sep='\s+')
#levels = levels.drop(levels.index[[-1]])
#levels.index = [i + 1 for i in levels.index]

class LEVELS(object):
    def __init__(self, file, ip=None, pseudolist=None, opticallist=None, sep='\s+'):
        self.levels = pd.read_csv(file, sep=sep)
        self.levels = self.levels.drop(self.levels.index[[-1]])
        self.levels.index = [i + 1 for i in self.levels.index]
        self.numlevels = max(self.levels.index)
        if ip is not None:
            self.aboveip = self.levels[self.levels['ENERGY(RYD)'] > ip]

        if pseudolist is not None and opticallist is not None:
            self.pseudolist = pseudolist            
            self.opticallist = opticallist
            self.pseudo = self.levels.query('CF in '+str(self.pseudolist))
            self.optical = self.levels.query('CF not in '+str(self.opticallist))
            if len(pseudolist+opticallist) != self.numlevels:
                print('WARNING: PSEUDOSTATE AND OPTICAL LIST DOES NOT MATCH TOTAL # LEVELS.')
            missing = [i+1 for i in list(range(self.numlevels)) if i+1 not in pseudolist+opticallist]
            if len(missing) > 0:
                print('THE FOLLOWING LEVEL #s NOT ACCOUNTED FOR: ',' '.join(missing))
            shared = [i for i in pseudolist if i in opticallist]
            if len(shared) > 0:
                raise ValueError("LEVELS ", " ".join(shared), "ARE PRESENT IN BOTH OPTICAL AND PSEUDO LISTS.")
                      
        elif pseudolist is not None:
            self.pseudolist = pseudolist
            self.pseudo = self.levels.query('CF in '+str(self.pseudolist))
            self.opticallist = [i+1 for i in list(range(self.numlevels)) if i+1 not in pseudolist]
            self.optical = self.levels.query('CF in '+str(self.opticallist))
        
        elif opticallist is not None:
            self.opticallist = opticallist
            self.optical = self.levels.query('CF in '+str(self.opticallist))
            self.pseudolist = [i+1 for i in list(range(self.numlevels)) if i+1 not in pseudolist]
            self.pseudo = self.levels.query('CF in '+str(self.pseudolist))
        return

    
    def makeenergyvsconfnumlist(self, leveldf, outfile):
        """
        Makes a csv of the energy vs. configuration number... useful for assesssing the
        quality of pseudostates.
        """
        df = pd.concat((leveldf['ENERGY(RYD)'], leveldf['CF']), axis=1)
        df.columns = ['#ENERGY(RYD)', 'CF']   #make easier to plot in gnuplot
        df.to_csv(outfile, sep='\t', index=False)
        return

    def getlevelenergies(self, TWOJ=None, L=None, S=None, P=None, CF=None):
        if L==None and S==None and P==None and CF==None:
            raise ValueError('MUST SPECIFY A 2J, L, S, P, OR CONF NUMBER!!!')
        else:
            qnums = locals()
            toquery = ' and '.join(['%s in %s' % (qnum, qnums[qnum]) for qnum \
                                                  in qnums if qnums[qnum] is not None and qnum != 'self'])
            result = self.levels.query(toquery)
            result = result['ENERGY(RYD)']
        return result

    def levelnumsfromconfs(self, conflist):
        return self.levels.query('CF in '+str(conflist)).index.tolist()

    def calculatestgsigminmax(self, initlevel, finallevelrange, elastic=False):
        level = 1
        count = 0
        if elastic == True:
            inc = len(self.levels)
            ntrsubt = initlevel + 1
        else:
            inc = len(self.levels) - 1
            ntrsubt = initlevel
        while level != initlevel+1:
            count += inc
            inc -= 1
            print(count)
            level += 1
        ntrmn = count - (len(self.levels) - finallevelrange[0])
        ntr = (len(self.levels) - (finallevelrange[1])) + count
        return ntrmn, ntr

    def makestgsiglist(self, conflist, initlevel, finallevelrange, outfile, elastic=False):
        levelnums = self.levelnumsfromconfs(conflist)
        ntrmn, ntr = self.calculatestgsigminmax(initlevel, finallevelrange, elastic=elastic)
        siglevelrange = list(range(ntrmn, ntr+1))
        forstgsig = [0]*len(siglevelrange)
        levelrange = list(range(finallevelrange[0], finallevelrange[1]+1))
        for levelnum in levelnums:
            if levelnum in levelrange:
                forstgsig[levelnum-finallevelrange[0]] = 1
        stgsiglist = open(outfile, 'w')
        for level in forstgsig:
            stgsiglist.write("%s\n" % level)
        print('dstgsig list printed.')
        return

    def makedirionsnlist(self, outfile):
        return


#configs = [1, 2, 3, 4, 7, 8]
#configs = [8]

#confs = levels.query('CF in '+str(configs))

#confavg = levels.groupby('CF')['ENERGY(RYD)'].mean()

#confs = [confavg.iloc[4], confavg.iloc[0], confavg.iloc[7], confavg.iloc[50]]

#print("CONFS",confs)

def get_term_egy_diff(energies, nistenergies):
    energies['S'] = energies['S'].apply(lambda x: abs(x))
    compare = energies.merge(nistenergies, on=['2J','S','L','CF'])
    termenergies = compare.loc[:,['CF','ENERGY(RYD)','NIST(RYD)']].groupby(by='CF').mean()
#    termenergies['DIFFSQ'] = termenergies.apply(lambda row: (row['NIST(RYD)'] - row['ENERGY(RYD)'])**2)
    termenergies['DIFFSQ'] = (termenergies['NIST(RYD)'] - termenergies['ENERGY(RYD)'])**2
    return termenergies

def get_sum_diff(energies, nistenergies):
    termenergies = get_term_egy_diff(energies, nistenergies)
    return termenergies['DIFFSQ'].sum()


#levels, nistlevels = sys.argv[1:3]


#levels = LEVELS(levels)
#nistlevels = LEVELS(nistlevels, sep=',')

#result = get_sum_diff(levels.levels, nistlevels.levels)

