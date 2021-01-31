import numpy as np


def read_xsec(file):
    xsecfile = open(file)
    xsecfile.readline()
    energies = []
    xsec = []
    for i in xsecfile:
        energies.append(float(i.split()[0]))
        xsec.append(float(i.split()[1]))
    return energies, xsec


class XSEC(object):
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
