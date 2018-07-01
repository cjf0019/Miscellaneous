# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 17:09:07 2017

@author: InfiniteJest
"""

import math
import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import SimpleRNN
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error

class MinimalRNNCell(keras.layers.Layer):

    def __init__(self, units, **kwargs):
        self.units = units
        self.state_size = units
        super(MinimalRNNCell, self).__init__(**kwargs)

    def build(self, input_shape):
        self.kernel = self.add_weight(shape=(input_shape[-1], self.units),
                  initializer='uniform',
                  name='kernel')
        self.recurrent_kernel = self.add_weight(
        shape=(self.units, self.units),
        initializer='uniform',
        name='recurrent_kernel')
        self.built = True

    def call(self, inputs, states):
        prev_output = states[0]
        h = K.dot(inputs, self.kernel)
        output = h + K.dot(prev_output, self.recurrent_kernel)
        return output, [output]

# Let's use this cell in a RNN layer:
#NOTE: Input must be of shape (batch_size,time_steps,input_dim)

def kinematic(a, t, v_0):
    return 0.5*(a)*t**2 + v_0*t 

def testkinematic(t):
    return kinematic(100, t, 0)

timeseq = list(range(101))
position = np.array(list(map(testkinematic, timeseq)))
#position = np.array(list(map(1, timeseq)))
maxpos = np.amax(position)
pixelscale = 500000

def getveclen(maxposition, scale):
    return int(math.ceil(maxposition/scale)) + 1

vectorlength = getveclen(maxpos, pixelscale)
print(maxpos)
print(vectorlength)


###CALCULATE FIRST VALUE SEPARATELY AND DON'T DIVIDE BY THE SPACE... ADD IT TO THE OTHERS
def convertfunctiontopointpixelvector(dist, fulldist, spacing):
    veclen = getveclen(fulldist, spacing)
    firstvecindex = int(dist // spacing)    
    probability = abs(((firstvecindex + spacing) - dist) / spacing)
    vector = np.zeros((veclen,))
    vector[firstvecindex] = probability
    try:
        vector[firstvecindex + 1] = abs(1 - probability)
    except:
        pass
    return vector

def testconvert(dist):
    return convertfunctiontopointpixelvector(dist, maxpos, pixelscale)

arrays = np.array(list(map(testconvert, position)))
fullarray = []
for array in arrays:
    fullarray.append(np.expand_dims(np.pad(array, 2, 'constant', constant_values=0), axis=0))
fullarray = np.array(fullarray)
traininp = fullarray[0:80]
trainout = fullarray[1:81]
testinp = fullarray[81:100]
testout = fullarray[82:]

#cell = MinimalRNNCell(vectorlength)
model = Sequential()
model.add(SimpleRNN(6, activation='relu', return_sequences=True, stateful=True,
               batch_input_shape=(1, 1, 6)))
model.add(SimpleRNN(6, activation='relu', return_sequences=True, stateful=True,
               batch_input_shape=(1, 1, 6)))

model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(traininp, trainout, shuffle=True, nb_epoch=200, batch_size=1)

weights = model.layers[0].get_weights()
model.save_weights('100a80pad1train.h5')

