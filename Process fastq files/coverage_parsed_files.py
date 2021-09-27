#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np

folderpath = '27jun18_parsed'

# get all the filenames
filenames = []
for fname in os.listdir(folderpath):
    if fname.find('.txt')>0:
        filenames.append(fname) 
print(filenames)       

total_counts = []
for file in filenames:
    content = pd.read_csv('{folder_name}/{file_name}'.format(folder_name=folderpath, file_name=file), header = None)
    total_counts.append(sum(content[1]))
print(total_counts)  
average = np.mean(total_counts)
std_dev = np.std(total_counts)
print(average, std_dev)

outputfile = 'Coverage in parsed files.csv'
df = pd.DataFrame({'file': filenames,
                   'total_reads': total_counts})
df.to_csv(outputfile,index=False)


# In[ ]:




