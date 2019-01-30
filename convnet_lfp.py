# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 14:07:15 2019

spike prediction from LFP by a simple 1D CNN
ref:
'https://blog.goodaudience.com/introduction-to-1d-convolutional-neural-networks-in-keras-for-time-sequences-3a7ff801a2cf'

@author: katsuhisa
"""
# libraries ====================
import scipy.io as sio
import h5py
import numpy as np
import keras
from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.layers import Conv1D, MaxPooling1D, GlobalAveragePooling1D


# load mat file ===========================================
stlfp0 = sio.loadmat('Z:/Katsuhisa/serotonin_project/LFP_project/Data/c2s/data/stlfp0_rc.mat')
stlfp1 = sio.loadmat('Z:/Katsuhisa/serotonin_project/LFP_project/Data/c2s/data/stlfp1_rc.mat')

n_time = np.shape(stlfp1)[1]

# 1D CNN model ===============================================
def oned_convnet(n_time):
    model = Sequential()
    model.add(Conv1D(64, 3, activation='relu', input_shape=(n_time, )))
    model.add(Conv1D(64, 3, activation='relu'))
    model.add(MaxPooling1D(3))
    model.add(Conv1D(128, 3, activation='relu'))
    model.add(Conv1D(128, 3, activation='relu'))
    model.add(GlobalAveragePooling1D())
    model.add(Dropout(0.5))
    model.add(Dense(1, activation='sigmoid'))
    
    model.compile(loss='binary_crossentropy',
                optimizer='rmsprop', metrics=['accuracy'])
    return model
#    print(model_m.summary())

# fit and evaluate ======================
callbacks_list = [
    keras.callbacks.ModelCheckpoint(
        filepath='best_model.{epoch:02d}-{val_loss:.2f}.h5',
        monitor='val_loss', save_best_only=True),
    keras.callbacks.EarlyStopping(monitor='acc', patience=1)
]
BATCH_SIZE = 400
EPOCHS = 50
history = model.fit(x_train, y_train,
                      batch_size=BATCH_SIZE,
                      epochs=EPOCHS,
                      callbacks=callbacks_list,
                      validation_split=0.2,
                      verbose=1)