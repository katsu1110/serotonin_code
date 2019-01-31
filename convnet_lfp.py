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
import csv
import multiprocessing

#import matplotlib.pyplot as plt
#import seaborn as sns

from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.layers import Conv1D, MaxPooling1D, GlobalAveragePooling1D
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc

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


# fit and evaluate =====================
model = oned_convnet(141)
loo = LeaveOneOut()
def fit_session(i):    
    # load data
    stlfp0 = sio.loadmat(l[i] + 'stlfp0.mat')
    stlfp1 = sio.loadmat(l[i] + 'stlfp1.mat')
    
    # baseline or drug
    acc_ses = [0]*2
    auc_ses = [0]*2
    for d in np.arange(2):
        # format data
        X = np.vstack((stlfp0['stlfp0'][0][d], stlfp1['stlfp1'][0][d]))
        X = np.expand_dims(X, axis=2)
        len0 = np.shape(stlfp0['stlfp0'][0][d])[0]
        len1 = np.shape(stlfp1['stlfp1'][0][d])[0]
        y = np.concatenate((np.zeros(len0, dtype=int), np.ones(len1, dtype=int)))
        ypredc = np.zeros(len(y))
        ypredp = np.zeros(len(y))
        
        # fit the model with leave-one-out
        for train_idx, test_idx in loo.split(X):
            # train and test datasets
            X_train, X_test = X[train_idx], X[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            
            # model fitting
            model.fit(X_train, y_train,
                  batch_size=32,
                  epochs=50,
                  verbose=1)
            
            # prediction
            ypredc[test_idx] = model.predict_classes(X_test).ravel()
            ypredp[test_idx] = model.predict(X_test).ravel()
            
        # model scores
        acc_ses[d] = np.sum(y==ypredc)/len(y)
        fpr, tpr = roc_curve(y, ypredp)
        auc_ses[d] = auc(fpr, tpr)
    
    # outputs
    return [l[i][-8:-1], acc_ses[0], acc_ses[1], auc_ses[0], auc_ses[1]]
        
pool = multiprocessing.Pool(multiprocessing.cpu_count())
results = pool.map(fit_session, (i for i in range(len(l))))
#results = Parallel(n_jobs=num_cores)(delayed(fit_session)(i) for i in range(len(l)))   
            
# save matrices 
with open("Z:/Katsuhisa/serotonin_project/LFP_project/Data/c2s/results.csv", "w") as outfile:
   writer = csv.writer(outfile, quoting=csv.QUOTE_ALL)
   writer.writerows(results)