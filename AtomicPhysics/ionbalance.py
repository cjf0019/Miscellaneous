# python3 function to read adf11 file
# returns numpy array with shape (z_max, ntemps, ndens)
# note that array is zero indexed, so z_max=1 is index 0
import os 
os.chdir("C:\\Users\\InfiniteJest\\Documents\\Python_Scripts\\Physics")
import numpy as np


def read_adf11(filename, z_max, ndens, ntemps):
    import numpy as np  
    with open(filename,'r') as f:  # open the file
        
        # determine the number lines in each density block
        dens_quot, dens_rem =  divmod(ndens,8)
        dens_lines = dens_quot
        if dens_rem > 0:
          dens_lines = dens_lines +1
        
        #determine the number of lines in each temperature block
        temp_quot, temp_rem =  divmod(ntemps,8)
        temp_lines = temp_quot
        if temp_rem > 0:
          temp_lines = temp_lines + 1
          
        line_dat = f.readline() # step forward one line
    
        while line_dat[1] != '-':  # skip past first line that begins with '-'
          line_dat = f.readline()
        
        # create an array with the density values
        dens_array = np.array(0)
        for i in range(dens_lines):
            dat = str.strip(f.readline()).split()
            dat = np.array(dat)
            dens_array = np.hstack((dens_array,dat))
            dens_array =  dens_array.astype(np.float)
            dens = 10 ** (dens_array[1:]) # standard ADAS format uses log10
        
        # create an array with the temperature values
        temp_array = np.array(0)
        for i in range(temp_lines):
            dat = str.strip(f.readline()).split()
            dat = np.array(dat)
            temp_array = np.hstack((temp_array,dat))
            temp_array =  temp_array.astype(np.float)
            temps = 10 ** (temp_array[1:]) # convert from log10
        
        # create a data array with shape ((z_max, ntemps, ndens))
        # and populate with the GCR coefficients from adf11 file
        adf11_dat = np.zeros((z_max, ntemps, ndens))
        for k in range(z_max):
            f.readline()
            qcd_array = dens_array
            for i in range(ntemps):
              ldat = qcd_array
              cdat = np.array(0)
              for j in range(dens_lines):
                dat = str.strip(f.readline()).replace('D','E').split()
                dat = np.array(dat)
                cdat = np.hstack((cdat,dat))
              qcd_array = np.vstack((ldat,cdat))
              qcd = (qcd_array.astype(np.float)[1:,1:])
            adf11_dat[k, :, :] = 10 ** qcd
        return(dens, temps, adf11_dat)


# %%
# call the function for both remobination (acd) and ionization (scd) files
dens, temps, acd_dat = read_adf11('acd89_fe.dat', 26, 26, 48)
dens, temps, scd_dat = read_adf11('scd89_fe.dat', 26, 26, 48)
    
import numpy as np
import matplotlib.pyplot as plt
#fig = plt.figure()
#plt.loglog(temps, acd_dat[0,:,0])
#fig.suptitle('ACD vs. Temperature')
#plt.xlabel('Temperature (K)')

#import matplotlib.pyplot as plt
#fig = plt.figure()
#plt.loglog(temps, scd_dat[0,:,0])
#fig.suptitle('SCD vs. Temperature')
#plt.xlabel('Temperature (K)')
   
def matrixline(A, state, tempind, densind):
    """
    Sets up a line of the A matrix given a charge state (0 onwards), temperature index,
    and density index (from adf11 file)
    """
    line = np.zeros(np.shape(A)[0])
    if state != 0:
        line[state-1] += scd_dat[state-1,tempind,densind]
        line[state] -= acd_dat[state-1,tempind,densind]
    
    if state != np.shape(A)[0]-1:
        line[state] -= scd_dat[state,tempind,densind]
        line[state+1] += acd_dat[state,tempind,densind]
    return line
                  

solutions = []
for temp in range(np.shape(temps)[0]):
    A = np.zeros((np.shape(acd_dat)[0]+1, np.shape(acd_dat)[0]+1)) #set up A matrix
    for state in range(len(A)):
        A[state] = matrixline(A, state, temp, 0)

    A[0] = np.ones((np.shape(acd_dat)[0])+1)

    Ainv = np.linalg.inv(A)
    steady = np.zeros(np.shape(A)[0])
    steady[0] = float(1)

    answer = np.matmul(Ainv, steady)
    solutions.append(answer)
    
solutions = np.array(solutions)
for i in solutions.T:
    plt.semilogx(temps,i)
plt.show()

