import adf04
import ecip
import numpy as np
from convertatomicunits import cmtoryd

wadf04 = adf04.adf04('adf04RyanFINAL')
d3ip = 0.0794
meta6sip = 0.0676



def ecipscalefxn(x, a):
    if isinstance(a, list):
        if len(a) == 3:
            return a[0]*x**2 + a[1]*x + a[2]
        elif len(a) == 2:
            return a[0]*x + a[1]
    else:
        return a*x
    return

ecipscales = []
outfile = open('ecipfactors', 'w')

#oddpar = [18.702,0.47084]
#evenpar = [23.106,-2.4895]

#BEFORE P ORB TEST
oddpar = [21.26,1.4295]
evenpar = [21.427,-1.922]

#AFTER P ORB TEST
oddpar = [21.427,-4.722]

for energy in wadf04.energies.split('\n')[1:]:
   print("ENERGYSPLIT",energy.split()[-1])
   ip = 0.577996 - cmtoryd(float(energy.split()[-1]))
   print("IP",ip) 
   conf = energy.split()[1][:-4]
   if conf in ['5d46s6p','5d36s26p','5p65d56p']:
       ecipsc = ecipscalefxn(ip, oddpar)
   elif conf in ['5d46s7s','5d46s6d','5p65d56s','5p65d56d','5s25p65d6','5p65d46s2']:
       ecipsc = ecipscalefxn(ip, evenpar)
   else:
       print("CONF OF TYPE ",conf," NEEDS DECISION")
       raise Exception

   outfile.write(str(np.round(ecipsc, 3)) + "\n")
   ecipscales.append(ecipsc)

outfile.close()
