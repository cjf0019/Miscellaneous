# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 15:09:24 2018
Inputs:
1) cross section file
2) ionization potential
3) max energy for Rost fit
4) Younger fit range (lower bound then upper bound)
5) bethe point
6) output file for fit

@author: InfiniteJest
"""
import numpy as np
import sys

xsecfile = sys.argv[1]
ionpot = float(sys.argv[2])
maxrostfitenergy = float(sys.argv[3])
youngerfitrange = [float(sys.argv[4]), float(sys.argv[5])]
bethe = float(sys.argv[6])
fitfile = sys.argv[7]


def xtrctoneenergycolstrength(file, transrange, ntr, units=13.6058, useXsec=True, lowerenergy=0.0, ion=1, w=1):
    colstrength = 0
    begin = transrange[0]
    end = transrange[1]
    line = file.readline()
    energy = float(line.split()[0])
    line = ' '.join(line.split()[1:])
    transcount = 6  #from first line
    numlinestobegin = int(np.ceil(begin/6)) - 1  #how many lines to skip before the beginning trans
    beginshift = int(begin % 6) - 1    #indicates where inside a line indexwise to start counting the transition
    numlinestoend = int(np.ceil(end/6))
    endshift = int(end % 6)
    numtrlines = int(np.ceil(ntr/6))
    if begin <= 6:
        colstrength += sum([float(i) for i in line.split()[beginshift:endshift]])
        addone = 0
    else:    
        for i in range(numlinestobegin):
            line = file.readline()
            transcount += len(line.split())
        colstrength += sum([float(i) for i in line.split()[beginshift:]])
        for i in range(numlinestoend - numlinestobegin -1):
            line = file.readline()
            transcount += len(line.split())
            colstrength += sum([float(i) for i in line.split()])
        line = file.readline()
        colstrength += sum([float(i) for i in line.split()[:endshift]])
        addone = 1
    for i in range(numtrlines - numlinestoend - addone):  #finish reading through entire energy point
        line = file.readline()
        transcount += len(line.split())
    if energy == 0:
        pass
    else:
        energy = (units*ion**2)*(energy-lowerenergy)   #returns incident energy relative to term threshold
        if useXsec == True:
            print(colstrength)
            colstrength = colstrength*(8.79735*10**(-17)/(w*energy/units))    #convert back to a cross section from collision strength
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
        self.numtermlines = int(np.ceil(self.terms/5))
        self.energies = []
        for line in range(self.numtermlines):
            self.energies = self.energies + [float(i) for i in omega.readline().split()]
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
                                  w = (2*abs(self.qnums[lowerterm-1][0])+1)*abs(self.qnums[lowerterm-1][1]))
            incidentenergies.append(energy)
            xsec.append(onexsec)
        return incidentenergies, xsec
        

def rost(x,p, ip): # the rost pittard function
    a,b,c,d=p
    return b*((float(1)/(((x-ip)/a)+0.88731))*((x-ip)/a/(((x-ip)/a)+ \
           0.88731))**1.127)/0.25899 +((a/c)*d*((float(1)/(((a/c)*(x-ip)/c)+ \
            0.29577))*((a/c)*(x-ip)/c/(((a/c)*(x-ip)/c)+0.29577))**3.381))/ \
            ((float(1)/((a/c)+0.29577))*((a/c)/((a/c)+0.29577))**3.381)

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
    def __init__(self, energies = None, xsec = None, colstrength = False, units= 13.6058, \
                 ion = None, ip = None):
        if energies is not None:
            self.energies = np.array(energies)
        if xsec is not None:
            self.xsec = np.array(xsec)
        if ion is not None:
            self.ion = ion   #ion state
        if ip is not None:
            self.ip = float(ip)   #ionization potential
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
	self.bethe = lowfit[2]
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

    def ionsn_plot_from_fit(self, fit = "both", energyrange = None, fitrange = None):
        """
        Currently set to plot a Rost fit, lower energies, and a Younger fit.
        """
        if energyrange is None:
            energyrange = [self.ip, 8*self.ip]
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
        ax.plot(self.energies,self.xsec,ls="",marker="x",color="blue",mew=2.0,label="Raw")
        morexlow=np.linspace(1.02*self.ip, 10*self.ip, 100)
        ax.plot(morexlow,rost(morexlow,self.lowfit[:4], self.ip),color="red",label="Fit")
        morexhigh=np.linspace(1.3*self.ip, energyrange[1], 100)
        ax.plot(morexhigh,younger(morexhigh,self.youngerfit, self.ip),color="red",label="Fit")
        ax.legend(loc=4)
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
        p = np.insert(p, 4, self.ip)
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
            p = np.insert(pnobethe, 2, bethepoint)
            return younger(energies, p, ionp)           
        def residuals(p): # array of residuals
            return fitxsec-youngerbethe(fitenergies,p,self.ip)
        def sum_residuals(p): # the function we want to minimize
            return sum(residuals(p)**2)          
        
        p0=[1,1,1] # initial parameters guess
        p,cov,infodict,mesg,ier=scimin.leastsq(residuals, p0,full_output=True) #traditional least squares fit
#        pwith=scimin.fmin_slsqp(sum_residuals,pwithout,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
#        p=scimin.fmin_slsqp(fitfunc,p0,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
        p = np.insert(p, 2, bethepoint)

        # plotting
        if showfitcomparison == True:
            import matplotlib.pyplot as plt
            ax=plt.figure().add_subplot(1,1,1)
            ax.plot(self.energies,self.xsec,ls="",marker="x",color="blue",mew=2.0,label="Raw")
            morex=np.linspace(1.3*self.ip,7*self.ip,100)
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
            thefit = '  '+'  '.join(["{:.5e}".format(i) for i in fit])+'\n'
            thefit = thefit.replace('  -', ' -')
            return thefit
	print("Rost Fit:",rost_fit)
	print("Younger Fit:",younger_fit)
        fitfile.write(fitline(rost_fit))
        fitfile.write(fitline(younger_fit))
        fitfile.close()
        print("Fits written to "+outfile)
        return

energies, rawxsec = read_xsec(xsecfile)
xsec = XSEC(energies=energies, xsec=rawxsec, ip=ionpot)
xsec.rost_fit(maxenergy=maxrostfitenergy)
xsec.younger_fit(bethe,energyrange=youngerfitrange)
xsec.ionsn_plot_from_fit()
xsec.write_fit_parameters(fitfile)
