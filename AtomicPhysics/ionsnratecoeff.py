import numpy as np
import scipy
from scipy import integrate
from omegautility import XSEC

rtcnst = 6.692353*10**(-11)
tconv = 0.00008617346  #K to eV
econv = 8065.479  #ev to cm

def rate_calc(temp, xsec, interpspace=0.005):
    maxenergy = np.max(xsec.energies)
    minenergy = np.min(xsec.energies)
    x = np.linspace(minenergy,maxenergy,num=interpspace)
    y = np.interp(x, xsec.energies, xsec.xsec)
    tmp = tconv * temp
    const = rtcnst /(tmp**1.5)
    def to_integrate(e):
        fxn = e * np.interp(e,xsec.energies,xsec.xsec) * np.exp(-e/tmp)
        print("FXN",fxn)
        return fxn
    rate = integrate.quad(to_integrate, minenergy, maxenergy)
    print("RATE, CONST", rate, const)
    rate = rate[0] * const
    return rate


class RATE:
    def __init__(self, tempgrid = None, xsec = None, xsecfile = None, ip = None):
        if tempgrid is not None:
            self.tempgrid = tempgrid
        if self.xsec is not None:
            self.xsec = xsec
        if self.xsecfile is not None:
            self.xsecfile = xsecfile
        if self.ip is not None:
            self.ip = ip

    def generate_rates(self):
        return

xsec = XSEC('eciptest')
print(rate_calc(5800,xsec))
print(rate_calc(11600,xsec))
print(rate_calc(29000,xsec))
