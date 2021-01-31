import numpy as np
import ecip
from xsec import read_xsec, XSEC
import sys

"""
Author: Connor Favreau
The base code for generating scaled ECIP fits to raw cross sections. A least squares fit is applied to a raw cross section to determine a multiplying constant 'a' for 'a * ecip_equation'. 

The code requires an f2py compiling of 'ecip_fxn.for' into ecip.so for use of the ECIP function.

Contains a child class of XSEC for ECIP cross sections. 
"""

def ecipscaled(egpt,a,ip,w):
    ecipscaled =  a*ecip.ecip(0,0,0,ip,w,egpt)
    return ecipscaled

def ecipegyshft(egpt,a,ip,w):
    engy = egpt-float(ip)
    if engy < 0:
        engy = 0.0
    return ecipscaled(engy,a,ip,w)

def calc_cont_egpts(ip,npts):
    emin = 0.001
    emax = 50.0 * ip
    de = (emax - emin)/float(npts)
    x = []
    energies = []
    for i in range(npts):
        egpt = emin + float(i)*de
        x.append(egpt)
        energies.append(x[-1]+ip)*de
    return x, energies

def write_ecip_xsec(energies,xsec,outfile):
    print("HERE!!!")
    print("OUTFILE",np.shape(xsec))
    sqxsec = np.squeeze(xsec)
    file = open(outfile, 'w')
    for i in range(len(energies)):
        file.write(str(np.round(energies[i],3)) + '\t' + str(np.round(sqxsec[i],3)) + '\n')
    file.close()

class ECIP(XSEC):
    def __init__(self, a=1.0, ip = None, w = 1, energies = None, xsec = None, colstrength = False, units=13.6058, ion=None):
        super(ECIP,self).__init__(ip=ip, energies=energies,xsec=xsec,ion=ion)
        self.a = float(a)
        self.w = int(w)

    def scaled_fit(self, label = None, energyrange = None, showfitcomparison = True, ecipprintoff=None):
        import scipy.optimize as scimin
        if self.ip is None or self.w is None:
            raise ValueError('MUST PROVIDE AN IONIZATION POTENTIAL AND NUMBER OF VALENCE ELECTRONS TO FIT.')
        a = self.a
#        if energyrange is not None:
#            fitrange = [list(self.energies).index(i) for i in self.energies if \
#                        (i <= energyrange[1]) and (i >= energyrange[0])]
#        else:
#            #fit on points from ip to 2*ip
#            fitrange = [list(self.energies).index(i) for i in self.energies if \
#                        (i <= 10*self.ip) and (i >= 2*self.ip)]
        
#        minindex = min(fitrange)
#        maxindex = max(fitrange)
        fitxsec = self.xsec[:]
        fitenergies = [i/self.units for i in self.energies]
        print("FITENERGIES",fitenergies)

        def residuals(a): # array of residuals
            ecips = np.zeros(np.shape(fitenergies))
            for i in range(len(fitenergies)):
                ecips[i] = ecipegyshft(fitenergies[i],a,self.ip,self.w)
            residuals = fitxsec-ecips
            return residuals

        def sum_residuals(p): # the function we want to minimize
            return sum(residuals(p)**2)          
        
        a0=a # initial parameter guess
        a,cov,infodict,mesg,ier=scimin.leastsq(residuals, a0,full_output=True) #traditional least squares fit
#        pwith=scimin.fmin_slsqp(sum_residuals,pwithout,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
#        p=scimin.fmin_slsqp(fitfunc,p0,bounds=[(self.ip, 2*self.ip)]) #add bounds between ip and 2*ip
        print("a",a,cov,mesg)
        print("IP",self.ip)

        # plotting
        if showfitcomparison == True:
            import matplotlib.pyplot as plt
            ax=plt.figure().add_subplot(1,1,1)
            ax.plot([i for i in self.energies],self.xsec,ls="",marker="x",color="blue",mew=2.0,label="Raw")
#            morex=np.linspace(0.0*self.ip,7*self.ip,100)
            morex = np.linspace(0.0,6,300) #400 eV
            morexipshift = (morex + self.ip)*self.units
            yreg = np.array([ecipscaled(i,1.0,self.ip,self.w) for i in morex])
            y = np.array([ecipscaled(i,a,self.ip,self.w) for i in morex])
            ax.plot(morexipshift,y,color="red",label="ECIP scaled")
            ax.plot(morexipshift,yreg,color="green",label='ECIP not scaled')
            ax.legend(loc=4)
            plt.xlabel("Incident Energy (eV)")
            plt.ylabel("Cross Section (Mb)")
            plt.savefig(label+".png")
            plt.show()

        if ecipprintoff is not None:
            write_ecip_xsec(morexipshift, y, ecipprintoff) 
            write_ecip_xsec(morexipshift, yreg, "noscale"+ecipprintoff)
            self.a = a
        return a

    def write_ecip_fit(self, outfile,label):
        fitfile = open(outfile, 'a')
        fitfile.write('\t'.join([str(self.ip),str(self.w),str(np.round(self.a[0],3)),label])+'\n')
        fitfile.close()
        print("Fit written to "+outfile)

def scale_ecip(ip,w,xsec_file,outfile,label,ecipprintoff=None):
    energies, rawxsec = read_xsec(xsec_file)
    xsec = ECIP(w = w, energies=energies, xsec=rawxsec, ip=ip)
    xsec.scaled_fit(label=label,ecipprintoff=ecipprintoff)
    xsec.write_ecip_fit(outfile,label)
    print("Fit of ",label,"finished!")

def fit_ecip_winfile(infile):
    file = open(infile)
    for line in file.readlines():
        if len(line.split()) == 6:
            ip,w,xsec_file,outfile,label,printfile = line.split()
            scale_ecip(ip,w,xsec_file,outfile,label,ecipprintoff=printfile)
        else:
            ip,w,xsec_file,outfile,label = line.split()
            scale_ecip(ip,w,xsec_file,outfile,label)

