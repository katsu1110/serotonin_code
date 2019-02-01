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
from sklearn.model_selection import StratifiedKFold
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

# model performance metric ==============
def modelperf(y, ypredc, ypredp):
    fpr, tpr = roc_curve(y, ypredp)
    return np.sum(y==ypredc)/len(y), auc(fpr, tpr)

# fit and evaluate =====================
model = oned_convnet(141)
n_split = 10
cv = StratifiedKFold(n_splits=n_split, shuffle=True, random_state=1220)
def fit_session(i):    
    # load data
    stlfp0 = sio.loadmat(l[i] + 'stlfp0.mat')
    stlfp1 = sio.loadmat(l[i] + 'stlfp1.mat')
    
    # baseline or drug
    accuracy = np.zeros((n_split, 2))
    auroc = np.zeros((n_split, 2))
    for d in np.arange(2):
        # format data
        X = np.vstack((stlfp0['stlfp0'][0][d], stlfp1['stlfp1'][0][d]))
        X = np.expand_dims(X, axis=2)
        len0 = np.shape(stlfp0['stlfp0'][0][d])[0]
        len1 = np.shape(stlfp1['stlfp1'][0][d])[0]
        y = np.concatenate((np.zeros(len0, dtype=int), np.ones(len1, dtype=int)))
        
        # fit the model with leave-one-out
        c = 0
        for train_idx, test_idx in cv.split(X, y):            
            # model fitting
            model.fit(X[train_idx], y[train_idx],
                  batch_size=32,
                  epochs=50,
                  verbose=1)
            
            # prediction and evaluation
            ypredc = model.predict_classes(X[test_idx]).ravel()
            ypredp = model.predict(X[test_idx]).ravel()
            accuracy[c,d], auroc[c,d] = modelperf(y[test_idx], ypredc, ypredp)
            c += 1
    
    # outputs
    fname = l[i][l[i].find('data')+5:-1]
    print(fname + ' processed!')
    mean_acc = np.mean(accuracy, axis=0)
    mean_auc = np.mean(auroc, axis=0)
    return [fname, mean_acc[0], mean_acc[1], mean_auc[0], mean_auc[1]]
        
results = fit_session(0)
#pool = multiprocessing.Pool(multiprocessing.cpu_count())
#results = pool.map(fit_session, (i for i in range(len(l))))
#results = Parallel(n_jobs=num_cores)(delayed(fit_session)(i) for i in range(len(l)))   
            
# save matrices 
with open("Z:/Katsuhisa/serotonin_project/LFP_project/Data/c2s/results.csv", "w") as outfile:
   writer = csv.writer(outfile, quoting=csv.QUOTE_ALL)
   writer.writerows(results)