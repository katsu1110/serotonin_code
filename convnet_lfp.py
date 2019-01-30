# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 14:07:15 2019

spike prediction from LFP by a simple 1D CNN
ref:
'https://blog.goodaudience.com/introduction-to-1d-convolutional-neural-networks-in-keras-for-time-sequences-3a7ff801a2cf'

@author: katsuhisa
"""
# libraries ====================
import glob
import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.layers import Conv1D, MaxPooling1D, GlobalAveragePooling1D

# datapath ===========================================
l = glob.glob(r'Z:/Katsuhisa/serotonin_project/LFP_project/Data/c2s/data/*/')

# 1D CNN model ===============================================
def oned_convnet(n_time):
    model = Sequential()
    model.add(Conv1D(64, 3, activation='relu', input_shape=(n_time, 1)))
    model.add(Conv1D(64, 3, activation='relu'))
    model.add(MaxPooling1D(3))
    model.add(Conv1D(128, 3, activation='relu'))
    model.add(Conv1D(128, 3, activation='relu'))
    model.add(GlobalAveragePooling1D())
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='sigmoid'))
#    print(model_m.summary())
    
    model.compile(loss='binary_crossentropy',
                optimizer='rmsprop', metrics=['accuracy'])
    return model


# fit and evaluate ======================
callbacks_list = [
    keras.callbacks.ModelCheckpoint(
        filepath='best_model.{epoch:02d}-{val_loss:.2f}.h5',
        monitor='val_loss', save_best_only=True),
    keras.callbacks.EarlyStopping(monitor='acc', patience=1)
]
BATCH_SIZE = 400
EPOCHS = 50
model = oned_convnet(141)
l = l[-2:]
print(l)
for c, fname in enumerate(l):
    # load data
    stlfp0 = sio.loadmat(l[c] + 'stlfp0.mat')
    stlfp1 = sio.loadmat(l[c] + 'stlfp1.mat')
    
    # baseline or drug
    for d in np.arange(2):
        # format data
        X = np.vstack((stlfp0['stlfp0'][0][d], stlfp1['stlfp1'][0][d]))
        X = np.expand_dims(X, axis=2)
        print(np.shape(X))
        len0 = np.shape(stlfp0['stlfp0'][0][d])[0]
        len1 = np.shape(stlfp1['stlfp1'][0][d])[0]
        y = np.concatenate((np.zeros(len0, dtype=int), np.ones(len1, dtype=int)))
        
        # fit the model
        history = model.fit(X, y,
                          batch_size=BATCH_SIZE,
                          epochs=EPOCHS,
                          callbacks=callbacks_list,
                          validation_split=0.2,
                          verbose=1)