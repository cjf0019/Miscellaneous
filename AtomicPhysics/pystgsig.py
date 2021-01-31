# -*- coding: utf-8 -*-
"""
Extract cross sections from an OMEGA file.

Created on Thu Mar 29 15:09:24 2018
Args:
1) The OMEGA file name
2) The lower transition number
3) The upper transition number
4) The lower term number
5) 1 for xsec, 0 for collision strength

@author: Connor Favreau
"""
import sys
import numpy as np

omegafile = sys.argv[1]
lowertrans = int(sys.argv[2])
uppertrans = int(sys.argv[3])
lowerterm = int(sys.argv[4])
xsecorcol = int(sys.argv[5])   #1 for xsec, 0 for col
if xsecorcol == 1:
    xsecorcol = True
else:
    xsecorcol = False

def xtrctoneenergycolstrength(file, transrange, ntr, units=13.6058, useXsec=True, lowerenergy=0.0, ion=1, w=1):
    """
    Inputs: file: the omega file to read
    transrange: a two-element list, giving the first and last transition to be extracted
    units: default is in eV
    lowerenergy: set to the lower term energy
    (default is 0.0 for incident energy relative to ground)
    ion: the charge of the atom + 1... used in the energy scaling (default is 1 for neutral)
    w: the statistical weight of the lower transition... only used if converting to xsec
    """

    colstrength = 0
    linecount = 0
    begin = transrange[0]
    end = transrange[1]
    line = file.readline()
    print(line)
    linecount += 1
    firsttrans = line.split()[1:]
    energy = float(line.split()[0])
    line = ' '.join(firsttrans)
    transcount = len(firsttrans)  #from first line
    print(transcount)
    numlinestobegin = int(np.floor(float(begin)/float(transcount))) #how many lines to skip before the beginning trans
    if begin % transcount == 0:
        numlinestobegin -= 1
    beginshift = begin - (begin // transcount)*transcount - 1  #indicates where inside a line indexwise to start counting the transition
    if beginshift == -1:
        beginshift = 5
    numlinestoend = int(np.floor(end/float(transcount)))
    endshift = end - (end // transcount)*transcount
    if end % transcount == 0:
        numlinestoend -= 1
    if endshift == 0:
        endshift = 6
    numtrlines = int(np.ceil(ntr/float(transcount)))

    #Start iteration over each line
    for i in range(numlinestobegin):  #skip to the beginning transition
        line = file.readline()
        linecount += 1
        transcount += len(line.split())
 
    if numlinestobegin - numlinestoend == 0:   #The beginning and end are on same line
        colstrength += sum([float(i) for i in line.split()[beginshift:endshift]])
        addone = 1

    else:
        colstrength += sum([float(i) for i in line.split()[beginshift:]])
        for i in range(numlinestoend - numlinestobegin - 1):
            line = file.readline()
            transcount += len(line.split())
            colstrength += sum([float(i) for i in line.split()])
            linecount += 1
        line = file.readline()
        linecount += 1
        colstrength += sum([float(i) for i in line.split()[:endshift]])
        addone = 1
    for i in range(numtrlines - numlinestoend - addone):  #finish reading through entire energy point
        line = file.readline()
        linecount += 1
        transcount += len(line.split())
    if energy == 0:
        pass
    else:
        energyconv = (units*ion**2)*(float(energy)-float(lowerenergy))   #returns incident energy relative to term threshold
        energy = (float(energy)-float(lowerenergy))*units
        if useXsec == True:
            #WILL OUTPUT IN Mb!
            colstrength = colstrength*(87.9735/(w*energyconv/units))    #convert back to a cross section from collision strength
    return energy, colstrength  


def stgsigspecify(file, stgsigfile, ntr, units=13.6058, lowerenergy=0.0, ion=1, w=1):
    """
    Similar to the above, extracts a cross section from an omega file, but is geared towards
    including only specific transitions.
    
    Inputs: file: the omega file to read
    transrange: a two-element list, giving the first and last transition to be extracted
    units: default is in eV
    lowerenergy: set to the lower term energy
    (default is 0.0 for incident energy relative to ground)
    ion: the charge of the atom + 1... used in the energy scaling (default is 1 for neutral)
    w: the statistical weight of the lower transition... only used if converting to xsec
    """
    dstgsig = open(stgsigfile)
    settings = dstgsig.readline()
    begin = int(re.search('\-[0-9]+', re.search('(ntrmn|NTRMN)(\=| \=)\-[0-9]*', settings).group(0)).group(0))
    end = int(re.search('[0-9]+', re.search('(ntran|NTRAN)(\=| \=)\-[0-9]*', settings).group(0)).group(0))
    numtrans = abs(end) - abs(begin) 
    useXsec = False
    if begin < 0:
        useXsec = True
        begin = abs(begin)
    
    linecount = 0
    line = file.readline()
    linecount += 1
    firsttrans = line.split()[1:]
    energy = float(line.split()[0])
    line = ' '.join(firsttrans)
    transcount = len(firsttrans)  #from first line
    numlinestobegin = int(np.floor(begin/transcount)) #how many lines to skip before the beginning trans
    if begin % transcount == 0:
        numlinestobegin -= 1
    beginshift = begin - (begin // transcount)*transcount - 1  #indicates where inside a line indexwise to start counting the transition
    if beginshift == -1:
        beginshift = 5
    numlinestoend = int(np.floor(end/transcount))
    endshift = end - (end // transcount)*transcount
    if end % transcount == 0:
        numlinestoend -= 1
    if endshift == 0:
        endshift = 6
    numtrlines = int(np.ceil(ntr/transcount))

    #Start iteration over each line
    for i in range(numlinestobegin):  #skip to the beginning transition
        line = file.readline()
        linecount += 1
        transcount += len(line.split())

    transitions = []    
    if numlinestobegin - numlinestoend == 0:   #The beginning and end are on same line
        transitions += [i for i in line.split()[beginshift:endshift]]
        addone = 1

    else:
        transitions += [i for i in line.split()[beginshift:]]
        for i in range(numlinestoend - numlinestobegin - 1):
            line = file.readline()
            transcount += len(line.split())
            transitions += [i for i in line.split()]
            linecount += 1
        line = file.readline()
        linecount += 1
        transitions += [i for i in line.split()[:endshift]]
        addone = 1
    for i in range(numtrlines - numlinestoend - addone):  #finish reading through entire energy point
        line = file.readline()
        linecount += 1
        transcount += len(line.split())

    #sum only those transitions marked with a 1 in the dstgsig file
    colstrength = 0
    for i in range(numtrans):
        keep = int(dstgsig.readline().strip())
        if keep == 1:
            colstrength += float(transitions[i])
        else:
            pass

    if energy == 0:
        pass
    else:
        energyconv = (units*ion**2)*(float(energy)-float(lowerenergy))   #returns incident energy relative to term threshold
        energy = (float(energy)-float(lowerenergy))*units
        if useXsec == True:
            #WILL OUTPUT IN Mb!
            colstrength = colstrength*(87.9735/(w*energyconv/units))    #convert back to a cross section from collision strength
    return energy, colstrength  



class OMEGA:
    def __init__(self, file):
        self.file = file
        omega = open(file)        
        self.nzion, self.nelc = [int(i.strip()) for i in omega.readline().split()]
        self.ion = self.nzion - self.nelc
        if self.ion == 0:
            self.ion = 1
        self.terms, self.nmpts, self.ntr = [int(i.strip()) for i in omega.readline().split()]
        self.qnums = []
        qnums = omega.readline().split()
        for term in np.array([2*(i + 1) for i in range(self.terms)]):
            self.qnums.append([int(i) for i in qnums[term-2:term]])     #store each term's 2 quantum numbers as tuples
        self.numtermlines = int(np.ceil(self.terms/float(5)))
        self.energies = []
        for line in range(self.numtermlines):
            self.energies = self.energies + [float(i) for i in omega.readline().split()]
        omega.close()
        return
        
    def xtrctxsec(self, transrange, lowerterm, useXsec=True):
        omega = open(self.file)
        for line in range(self.numtermlines + 3):   #skip to the cross sections
            omega.readline()
            
        incidentenergies = []
        xsec = []
        for i in range(self.nmpts):
            energy, onexsec = xtrctoneenergycolstrength(omega, transrange, self.ntr, useXsec=useXsec, \
                                  lowerenergy=self.energies[lowerterm-1], ion=self.ion, \
                                  w = (2*abs(self.qnums[lowerterm-1][1])+1)*abs(self.qnums[lowerterm-1][0]))
            incidentenergies.append(energy)
            xsec.append(onexsec)
        return incidentenergies, xsec
        

def rost(x,p, ip): # the rost pittard function
    a,b,c,d=p
    return b*((float(1)/(((x-ip)/a)+0.88731))*((x-ip)/a/(((x-ip)/a)+ \
           0.88731))**1.127)/0.25899 +((a/c)*d*((float(1)/(((a/c)*(x-ip)/c)+ \
            0.29577))*((a/c)*(x-ip)/c/(((a/c)*(x-ip)/c)+0.29577))**3.381))/ \
            ((float(1)/((a/c)+0.29577))*((a/c)/((a/c)+0.29577))**3.381)

def cheb(x, umin, umax, ncheb, a, ip):
    """
    x is the cross section, umin the threshold x/ip, umax the max x/ip, ncheb the
    number of parameters, a the parameter list, and ip the ionization potential.
    """
    u = x/ip
    z = ((np.log(u)-np.log(umin))-(np.log(umax)-np.log(u)))/(np.log(umax)- \
         np.log(umin))
    t = [z, 2*z**2 - 1]
    nply = ncheb
    for i in range(3, ncheb):
        t.append(2*z*t[-1]-t[-2])
    usg = a[0]/2
    for i in range(2, ncheb+1):
        usg += a[-1]*t[-2]
    return usg/u

def younger(x, p, ip):
    a,b,c,d=p
    u = x/ip
    return (a*(float(1)-float(1)/u)+b*(float(1)-float(1)/u)**float(2) + \
         c*np.log(u) + d*np.log(u)/u)/(u*ip**float(2))

def read_xsec(file):
    xsecfile = open(file)
    xsecfile.readline()
    energies = []
    xsec = []
    for i in xsecfile:
        energies.append(float(i.split()[0]))
        xsec.append(float(i.split()[1]))
    return energies, xsec


class XSEC:
    def __init__(self, file = None, energies = None, xsec = None, colstrength = False, units= 13.6058, \
                 ion = None, ip = None):
        if file is not None:
            self.energies, self.xsec = read_xsec(file)
        if energies is not None:
            self.energies = np.array(energies)
        if xsec is not None:
            self.xsec = np.array(xsec)
        if ion is not None:
            self.ion = ion   #ion state
        if ip is not None:
            self.ip = ip   #ionization potential
        self.colstrength = colstrength   #is this a collision strength or cross section
        self.units = units
        return

    def get_xsec_from_file(self, file):
        energies, xsec = read_xsec(file)
        self.energies = np.array(energies)
        self.xsec = np.array(xsec)
        print('Cross Section Read from File')
        return

    def get_fits_from_file(self, file, updateip = False):
        """
        Retrieve Rost-Pittard and Younger fits from file. If updateip == True,
        change self.ip to be the value in the file.
        """
        fitfile = open(file)
        ip = fitfile.readline().split()[-1]
        if updateip == True:
            self.ip = float(ip)
        self.fitstart, self.fitend = [self.ip*float(i) for i in fitfile.readline().split()] #the energy limits of the fits...
        lowfit = fitfile.readline()
        lowfit = [float(i) for i in lowfit.replace('D', 'E').split()[:4]] 
        if len(lowfit) == 6:
            print('CONTAINS A CHEBYCHEV FIT')
            self.rostorcheb = 'cheb'
        else:
            print('CONTAINS A ROST PITTARD FIT')
            self.rostorcheb = 'rost'
        youngfit = fitfile.readline()
        youngfit = [float(i) for i in youngfit.replace('D', 'E').split()]
        self.lowfit = lowfit
        self.youngfit = youngfit
        print('The fits are ', lowfit, youngfit)
        return

    def plot_raw(self, xsecorcolstrength = 'xsec'):
        """
        Plot the cross section. Default assumes the data already is a cross section
        and performs no conversions. If the data is a col. strength, the 'xsec' option will convert
        the plot and vice versa for setting to 'colstrength.'  <----NOT IMPLEMENTED YET!
        """
        import matplotlib.pyplot as plt
        fig = plt.figure()
        plt.plot(self.energies[:-1], self.xsec[:-1])
        fig.suptitle('Cross Section')
        plt.xlabel('Incident Energy (eV)')
        plt.ylabel('Cross Section')
        return

    def ionsn_plot_from_fit(self, lowfit = "rost", highfit = "younger", energyrange = None, fitrange = None):
        """
        Currently set to plot a Rost fit, lower energies, and a Younger fit.
        """
        if energyrange is None:
            energyrange = [self.ip, 5*self.ip]

        if fitrange is not None:
            self.fitstart = float(fitrange[0])
            self.fitend = float(fitrange[1])
        else:
            self.fitstart = 1.02*self.ip
            self.fitend = 2.00*self.ip
        fitstart = self.fitstart
        fitend = self.fitend
        
        import matplotlib.pyplot as plt
        ax=plt.figure().add_subplot(1,1,1)
        morexlow=np.linspace(fitstart, fitend, 100)
        if lowfit == "rost":
            ax.plot(morexlow,rost(morexlow,self.lowfit[:4], self.ip),color="red",label="Fit")
        elif lowfit == "cheb":
            ax.plot(morexlow,cheb(morexlow,self.fitstart, self.fitend, 6, \
                                  self.lowfit, self.ip),color="red",label="Fit")

        morexhigh=np.linspace(fitend, energyrange[1], 100)
        ax.plot(morexhigh,younger(morexhigh,self.youngerfit, self.ip),color="red",label="Fit")
        ax.legend(loc=2)
        plt.show()
    
    def rost_fit(self, maxenergy=None, showfitcomparison = True):
        """
        Fits lower energies (below 2 * ip) to a Rost-Pittard fit.
        """
        import scipy.optimize as scimin
        if self.ip is None:
            raise ValueError('MUST PROVIDE AN IONIZATION POTENTIAL TO FIT.')

        if maxenergy is not None:
            fitrange = [list(self.energies).index(i) for i in self.energies if \
                        (i <= maxenergy) and (i >= self.ip)]
        else:
            #fit on points from ip to 2*ip
            fitrange = [list(self.energies).index(i) for i in self.energies if \
                        (i <= 2*self.ip) and (i >= self.ip)]
        
        minindex = min(fitrange)
        maxindex = max(fitrange)
        fitxsec = self.xsec[minindex:maxindex+1]
        fitenergies = self.energies[minindex:maxindex+1]

        def residuals(p): # array of residuals
            return fitxsec-rost(fitenergies,p,self.ip)
        def sum_residuals(p): # the function we want to minimize
            return sum(residuals(p)**2)          
        p0=[1,1,1,1] # initial parameters guess
        p,cov,infodict,mesg,ier=scimin.leastsq(residuals, p0,full_output=True) #traditional least squares fit
#        pwith=scimin.fmin_slsqp(sum_residuals,pwithout,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
#        p=scimin.fmin_slsqp(fitfunc,p0,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip

        # plotting
        if showfitcomparison == True:
            import matplotlib.pyplot as plt
            ax=plt.figure().add_subplot(1,1,1)
            ax.plot(self.energies,self.xsec,ls="",marker="x",color="blue",mew=2.0,label="Raw")
            morex=np.linspace(self.ip,2*self.ip,100)
            ax.plot(morex,rost(morex,p, self.ip),color="red",label="Fit")
            ax.legend(loc=2)
            plt.show()
        np.insert(p, 4, self.ip)
        self.lowfit = p
        return p

    def younger_fit(self, bethepoint, energyrange = None, showfitcomparison = True):
        """
        Fits higher energies (above 2 * ip) to a Younger fit. The Bethe point 
        corresponds to the third parameter and can be obtained from a DW calculation.
        """
        import scipy.optimize as scimin
        if self.ip is None:
            raise ValueError('MUST PROVIDE AN IONIZATION POTENTIAL TO FIT.')

        if energyrange is not None:
            fitrange = [list(self.energies).index(i) for i in self.energies if \
                        (i <= energyrange[1]) and (i >= energyrange[0])]
        else:
            #fit on points from ip to 2*ip
            fitrange = [list(self.energies).index(i) for i in self.energies if \
                        (i <= 10*self.ip) and (i >= 2*self.ip)]
        
        minindex = min(fitrange)
        maxindex = max(fitrange)
        fitxsec = self.xsec[minindex:maxindex+1]
        fitenergies = self.energies[minindex:maxindex+1]

        def youngerbethe(energies,pnobethe,ionp):
            p = pnobethe.insert(2, bethepoint)
            return younger(energies, p, ionp)           
        def residuals(p): # array of residuals
            return fitxsec-youngerbethe(fitenergies,p,self.ip)
        def sum_residuals(p): # the function we want to minimize
            return sum(residuals(p)**2)          
        
        p0=[1,1,1] # initial parameters guess
        p,cov,infodict,mesg,ier=scimin.leastsq(residuals, p0,full_output=True) #traditional least squares fit
#        pwith=scimin.fmin_slsqp(sum_residuals,pwithout,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
#        p=scimin.fmin_slsqp(fitfunc,p0,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
        np.insert(p, 2, bethepoint)
        # plotting
        if showfitcomparison == True:
            import matplotlib.pyplot as plt
            ax=plt.figure().add_subplot(1,1,1)
            ax.plot(self.energies,self.xsec,ls="",marker="x",color="blue",mew=2.0,label="Raw")
            morex=np.linspace(self.ip,2*self.ip,100)
            ax.plot(morex,younger(morex,p, self.ip),color="red",label="Fit")
            ax.legend(loc=2)
            plt.show()
        self.youngerfit = p
        return p

    def write_fit_parameters(self, outfile, fitrange=None, bethe=None):
        fitfile = open(outfile, 'a')
        if self.lowfit is None:
            self.rost_fit()
        rost_fit = self.lowfit
        if self.youngerfit is None:
            self.younger_fit(bethe)
        younger_fit = self.youngerfit
        fitfile.write('    1    1    '+str(len(self.lowfit))+'   '+str(self.ip)+'\n')
        if fitrange is not None:
            fitfile.write('     '+str(np.round(fitrange[0],2))+'    '+ \
                                               str(np.round(fitrange[1],2))+'\n')
        else:
            fitfile.write('     1.02    2.0\n')
        def fitline(fit):
            thefit = '  '+'  '.join(["{:.5E}".format(i) for i in fit])+'\n'
            thefit = thefit.replace('  -', ' -')
            return thefit
        fitfile.write(fitline(rost_fit))
        fitfile.write(fitline(younger_fit))
        fitfile.close()
        print("Fits written to "+outfile)
        return

    def write_xsec(self,file):
        xsecfile = open(file, 'a')
        xsecfile.write('#Energy,  XSec\n')
        for i,j in zip(self.energies, self.xsec):
            xsecfile.write(str(i)+"\t"+str(j)+"\n")
        print("Cross Section written to",str(file))
        xsecfile.close()
        return
    

def scale_by_n4(xsec, ni, nf, ipi, ipf):
    """
    Xsec is a class XSEC object, and ni the original cross section's n number,
    nf the final cross section's n number. Returns another XSEC class with the new
    nf cross section. 
    """
    energies = (ipf/ipi)*xsec.energies
    newxsec = (nf**4/ni**4)*xsec.xsec
    return XSEC(energies=energies, xsec=newxsec, ip=ipf)



omega = OMEGA(omegafile)        
incidentenergies, xsec = omega.xtrctxsec([lowertrans, uppertrans], lowerterm, useXsec = xsecorcol)

#Write the raw data to a file
data = np.concatenate((np.expand_dims(np.array(incidentenergies), axis=1), \
                       np.expand_dims(np.array(xsec), axis=1)), axis=1)
file = open('sg.dat', 'w')
file.write('#Energy    Xsec')
for i in data:
    file.write('\t'.join([str(j) for j in i])+'\n')
file.close()

#import matplotlib.pyplot as plt
#fig = plt.figure()
#plt.plot(energies[:-1], xsec[:-1])
#fig.suptitle('Cross Section')
#plt.xlabel('Incident Energy (eV)')
#plt.ylabel('Cross Section')
    
    
