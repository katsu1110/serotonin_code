# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 12:14:59 2019

@author: katsuhisa
"""

# Save a dictionary into a pickle file.
import os
import platform
import numpy as np
from numpy.core import multiarray
from scipy.io import loadmat
import pickle

def load_data(filepath):
	"""
	Loads data in either pickle or MATLAB format.
	@type  filepath: string
	@param filepath: path to dataset
	@rtype: list
	@return: list of dictionaries containing the data
	"""

	if filepath.lower().endswith('.mat'):
		data = []
		data_mat = loadmat(filepath)

		if 'data' in data_mat:
			data_mat = data_mat['data'].ravel()

			for entry_mat in data_mat:
				entry = {}

				for key in entry_mat.dtype.names:
					entry[key] = entry_mat[key][0, 0]

				for key in ['calcium', 'spikes', 'spike_times']:
					if key in entry:
						entry[key] = entry[key].reshape(1, entry[key].size)
				if 'fps' in entry:
					entry['fps'] = float(entry['fps'])
				if 'cell_num' in entry:
					entry['cell_num'] = int(entry['cell_num'])

				data.append(entry)

		elif 'predictions' in data_mat:
			for predictions in data_mat['predictions'].ravel():
				data.append({'predictions': predictions.reshape(1, predictions.size)})

		return data

if 'Win' in platform.platform():
    mypath = '//172.25.250.112/nienborg_group/'
else:
    mypath = '/gpfs01/nienborg/group/'

datapath = mypath + '/Katsuhisa/serotonin_project/LFP_project/Data/c2s/data/'
fnames = ['data_train', 'data_test', 'data_drug']

for l in os.listdir(datapath):
    for f in fnames:
        p = datapath + l + '/' + f
        data = load_data(p + '.mat')
        with open(p + '.pck', 'wb') as f:
            pickle.dump(data, f, protocol=2)
        print(p + ' saved as pickle')