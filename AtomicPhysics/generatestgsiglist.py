# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 12:30:29 2018
Reads in a TERMS file and separates pseudostate configs from nonpseudostate
configs (User must manually indicate which are not pseudostates).
Must specify the TERMS file and list of configuration numbers (from the AUTOSTRUCTURE run) at the bottom of this code.
@author: Connor Favreau
"""

import pandas as pd

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

    def termnumsfromconfs(self, conflist):
        return self.terms.query('CF in '+str(conflist)).index.tolist()

    def calculatestgsigminmax(self, initterm, finaltermrange, elastic=False):
        term = 1
        count = 0
        if elastic == True:
            inc = len(self.terms) 
            ntrsubt = initterm + 1
        else:
            inc = len(self.terms) - 1
            ntrsubt = initterm
        while term != initterm+1:
            count += inc
            inc -= 1
            print(count)
            term += 1
        ntrmn = count - (len(self.terms) - finaltermrange[0])
        ntr = (len(self.terms) - (finaltermrange[1])) + count
        return ntrmn, ntr

    def makestgsiglist(self, conflist, initterm, finaltermrange, outfile, elastic=False):
        termnums = self.termnumsfromconfs(conflist)
        ntrmn, ntr = self.calculatestgsigminmax(initterm, finaltermrange, elastic=elastic)
        print(ntrmn, ntr)
        sigtermrange = list(range(ntrmn, ntr+1))
        print(sigtermrange)
        forstgsig = [0]*len(sigtermrange)
        termrange = list(range(finaltermrange[0], finaltermrange[1]+1))
        for termnum in termnums:
            if termnum in termrange:
                forstgsig[termnum-finaltermrange[0]] = 1
        stgsiglist = open(outfile, 'w')
        for term in forstgsig:
            stgsiglist.write("%s\n" % term)
        print('dstgsig list printed.')
        return

    def makedirionsnlist(self, outfile):
        return




#conflist = list(range(46,87))
conflist = list(range(5,42))
#19325
terms = TERMS('TERMSw0')
#terms.makestgsiglist(conflist, 2, [1486,15756], 'stgsiglist', elastic=False)
terms.makestgsiglist(conflist, 2, [16326,31917], 'stgsiglist', elastic=False)


