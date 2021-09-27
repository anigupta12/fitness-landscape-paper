#!/usr/bin/env python
# coding: utf-8

# In[3]:


"""To calculate frequency of OmpF evolution from landscape simulations in simulations_main.py """

import numpy as np
from numpy.random import binomial as nbinom
import random
import math

from timeit import default_timer as timer

# from numpy import nonzero as nnonzero
from six.moves import range as srange
import csv
import pandas as pd
from matplotlib import pyplot as plt

import plotly.plotly as py
import plotly.figure_factory as ff

from IPython.display import set_matplotlib_formats
# set_matplotlib_formats('pdf', 'svg')

time_step = 24*2
generations = 20*time_step
interval=10
switch_time_array = np.array(range(0,generations+int(generations/interval),int(generations/interval)))
# switch_generation = 864
switch_time_array = switch_time_array[1:-1]

total_number_runs = 300
threshold_OmpF = [5000]

batch_size = 30
batch_array = np.array(range(0,total_number_runs,batch_size))

output_data = {}
np.set_printoptions(precision=5)

for t in range(len(threshold_OmpF)):
    print('threshold {}'.format(threshold_OmpF[t]))
    OmpF_freq = np.zeros([len(batch_array),switch_time_array.size])
    batch_number = 0
    for batch in batch_array:
        count = 0
        for switch_generation in switch_time_array:
            y = np.zeros([batch_size,1])
            for n in range(0,batch_size,1):    
                filename = 'OmpF-abundance_Switch-generation_{}_run-number_{}.csv'.format(switch_generation,n+batch)
                data = pd.read_csv(filename)
                data = data.to_numpy() # np.array(data)
                data = data[:,1:] # to remove the first column of row numbers from the data frame
                OmpF_abundance = data[data>threshold_OmpF[t]]  # get the values in data array which is greater than threshold
                indices_OmpF_abundance = np.transpose(np.nonzero(data>threshold_OmpF[t]))  # get the indices in data array greater than a threshold. from https://docs.scipy.org/doc/numpy/reference/generated/numpy.nonzero.html
                y[n] = indices_OmpF_abundance.shape[0]
            total_OmpF_crossings = y[y>1] # total number of entries in the data which are above OmpF threshold for different runs in a switching simulations
#             print('threshold {}, switching time {}, total OmpF crossings- '.format(threshold_OmpF[t],switch_generation))
#             print(np.around(total_OmpF_crossings,3))
            OmpF_freq[batch_number, count] = len(total_OmpF_crossings)/batch_size
            count = count+1
        print('batch_number {}, OmpF_freq {}'.format(batch_number+1,OmpF_freq[batch_number,:]))
        batch_number = batch_number + 1
    #    output_data.update({str(switch_time_array) : OmpF_freq[t,:]})
    column_names = list(map(str,switch_time_array))
    row_index = list(range(1,len(batch_array)+1))
    df = pd.DataFrame(OmpF_freq, columns = column_names, index = row_index)
    filename = 'analysis_threshold-{}.csv'.format(threshold_OmpF[t])
    df.to_csv(filename)

# plt.bar(switch_time_array,OmpF_freq)
# plt.show()

