# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:27:17 2018
Using Metropolis Monte Carlo sampling, calculates the average energy and magnetization,
as well as their fluctuations (<E^2> and <mu^2>) for a 32 x 32 Ising model.
@author: InfiniteJest
"""

import numpy as np


lowinit = np.full((32, 32), -1)     #for close to T = 0... assumes all spins are at -1 (lowest energy)
highinit = np.random.uniform(low=-1, high=1, size=(32, 32))    #for high T... assumes spins are all random

def oneornegone(value):    #rounds random numbers to either 1 or -1
    if value < 0:
        newvalue = -1
    else:
        newvalue = 1
    return newvalue
    
oneornegone = np.vectorize(oneornegone)
highinit = oneornegone(highinit)     #convert highinit random numbers to 1 or -1

def one_spin_energy(xind, yind, matrix):
    """
    Calculates the energy at one location given its spin and the adjacent spins.
    xind and yind are the location in the matrix (a 32 x 32 configuration of spins).
    """
    spin = matrix[xind][yind]
    if xind == len(matrix[:][yind]) - 1:  #for periodic condition; for last index, go back to first index
        one = matrix[0][yind]
    else:
        one = matrix[xind+1][yind]
    two = matrix[xind-1][yind]
    if yind == len(matrix[xind][:]) - 1:
        three = matrix[xind][0]
    else:
        three = matrix[xind][yind+1]
    four = matrix[xind][yind-1]
    return -(spin*one + spin*two + spin*three + spin*four)


def energy_per_spin(matrix):    #calculate the energy at each position
    energymatrix = np.zeros(np.shape(matrix))
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            energymatrix[i][j] = one_spin_energy(i, j, matrix)
    return energymatrix

def calculate_avg_and_fluc(matrix):
    """
    Will calculate both <E^2> and <mu^2> depending on the matrix given.
    Backcalculates from the formula std = sqrt(<x^2> - <x>^2) . 
    """
    avg = np.mean(matrix)
    std = np.std(matrix)
    fluc = np.sqrt(std**2 + avg**2)
    return avg, fluc


def metropolisfxn(xn, xnp1, beta):   #calculate the probability of an energy at a given temperature
    return np.exp(-(xnp1 - xn)/beta)


def runmetropolisalgorithm(N, s, init, beta):
    """
    Performs the Metropolis algorithm to generate a Maxwellian sampling. 
    N represents the number of steps taken for sampling. At each step, a new (random)
    configuration is compared to the current configuration. If it is of less energy,
    it is automatically accepted. Otherwise, it is accepted with a Maxwellian probability
    based on the change in energy.
    """
    r = init
    initenergy = energy_per_spin(r)
    avginitenergy = calculate_avg_and_fluc(initenergy)[0]
    p = metropolisfxn(avginitenergy, avginitenergy, beta)

    samples = []
    for i in range(N):
        rn = oneornegone(np.random.uniform(low=-1, high=1, size=(32, 32)))  #try a new random configuration
        energyn = energy_per_spin(rn)
        avgenergyn = calculate_avg_and_fluc(energyn)[0]
        pn = metropolisfxn(avginitenergy, avgenergyn, beta)   #how likely the new configuration is from the older one at T=beta
        if avgenergyn <= avginitenergy:         #Accept all distributions with lower overall energy
            p = pn          #update probability
            r = rn          #update the configuration
            initenergy, avginitenergy = energyn, avgenergyn
        else:
            #accept a subset of higher energy configurations with Maxwellian probability
            u = np.random.rand()
            if u < pn:
                p = pn
                r = rn
        if i % s == 0:   #sample at only the 10th step
            samples.append(r)
    return samples
    

#Find <E>, <E**2>, <mu>, and <mu**2> at each temperature
temperatures = [0.01, 0.1, 1.0, 2.0, 3.0, 10]
avgenergies = []
energyfluc = []
avgmags = []
magfluc = []

for temp in temperatures:
    if temp != 10:
        sample = runmetropolisalgorithm(10000, 10, lowinit, temp)   #use low init for lower temperatures
    else:
        sample = runmetropolisalgorithm(10000, 10, highinit, temp)  #use high init for high temperature
    #energies
    sampleenergies = [calculate_avg_and_fluc(energy_per_spin(i))[0] for i in sample]
    avgenergy = np.mean([sampleenergies[i] for i in range(len(sampleenergies))])
    avgenergies.append(avgenergy)
    energyfluc.append(np.sqrt(np.std(sampleenergies)**2 + avgenergy**2))
    #magnetic moments
    samplemags = [calculate_avg_and_fluc(i)[0] for i in sample]
    avgmag = np.mean([samplemags[i] for i in range(len(samplemags))])
    avgmags.append(avgmag)
    magfluc.append(np.sqrt(np.std(sampleenergies)**2 + avgmag**2))
    

temperatures, avgenergies, energyfluc, avgmags, magfluc = np.array(temperatures), \
    np.array(avgenergies), np.array(energyfluc), np.array(avgmags), np.array(magfluc)

#Plot each as a function of temperature
import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot(temperatures, avgenergies)
fig.suptitle('Average Energy vs. Temperature')
plt.xlabel('Temperature (kT/J)')
plt.ylabel('Energy (J)')
fig.savefig('AvgEnergy.jpg')

fig2 = plt.figure()
plt.plot(temperatures, energyfluc)
fig2.suptitle('<E^2> vs. Temperature')
plt.xlabel('Temperature (kT/J)')
plt.ylabel('<E^2>')
fig2.savefig('EnergyFluc.jpg')

fig3 = plt.figure()
plt.plot(temperatures, avgmags)
fig.suptitle('Average Magnetization vs. Temperature')
plt.xlabel('Temperature (kT/J)')
plt.ylabel('Magnetization')
fig3.savefig('AvgMagnetization.jpg')

fig4 = plt.figure()
plt.plot(temperatures, energyfluc)
fig4.suptitle('<mu^2> vs. Temperature')
plt.xlabel('Temperature (kT/J)')
plt.ylabel('<mu^2>')
fig4.savefig('MagFluc.jpg')


#Write the raw data to a file
import pandas as pd
temperatures = np.expand_dims(temperatures, axis=1)
avgenergies = np.expand_dims(avgenergies, axis=1)
energyfluc = np.expand_dims(energyfluc, axis=1)
avgmags = np.expand_dims(avgmags, axis=1)
magfluc = np.expand_dims(magfluc, axis=1)

data = np.concatenate((temperatures, avgenergies, energyfluc, avgmags, magfluc), axis=1)
df = pd.DataFrame(data)
df.columns = ['Temperature', 'Average Energy', '<E^2>', 'Average Magnetization', '<mu^2>']

df.to_csv('MonteCarloIsingResults', sep='\t', index=False)

