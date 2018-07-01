# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 10:58:50 2018

@author: InfiniteJest
"""

from keras.models import Sequential
from keras.layers import Dense
from keras.callbacks import EarlyStopping
import numpy as np
import os
os.chdir('C:\\Users\\InfiniteJest\\Documents\\Python_Scripts\Physics')

file = open('orbenergycomp.dat')
file.readline()
inp = []
out = []
for line in file:
    inp.append([float(i) for i in line.split()[:2]])
    out.append([float(i) for i in line.split()[2:]])
 
inp = np.array(inp)
out = np.array(out)
    
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('Qt4Agg')
#import matplotlib.pyplot as plt


fig = matplotlib.pyplot.figure()
ax = fig.add_subplot(111, projection='3d')
x = inp.T[0]
y = inp.T[1]
X, Y = np.meshgrid(x, y)
zs = out.T[2]
#Z = zs.reshape(X.shape)

ax.plot_surface(X, Y, zs)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()
    


dnn = Sequential()
#early_stopping = EarlyStopping(monitor='val_loss', patience=3)
dnn.add(Dense(3, input_dim=2, init='uniform', activation='relu'))
#dnn.add(Dense(3, activation='sigmoid'))
#dnn.add(Dense(100))
dnn.compile(loss='mean_squared_error', optimizer='adam')
dnn.fit(inp, out, shuffle=True, validation_split=0.2, nb_epoch=100, batch_size=1)