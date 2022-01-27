
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')

import numpy as np
import matplotlib.pyplot as plt
import csv
import struct


# In[4]:


#File Converter#
file = open('p163241.03500_e5059', 'r')
pfile = file.readlines()
pfilecon = []
chunked_list = list()
chunk_size = 257
for i in range(len(pfile)):
    row = pfile[i].split()
    pfilecon.append(row)
for j in range(0,len(pfilecon),chunk_size):
    if int(float(pfilecon[j][0])) > 3:
        chunked_list.append(pfilecon[j:j+chunk_size])
    else:
        chunked_list.append(pfilecon[j:j+4])
        
#print(chunked_list[0][1][1])
print(chunked_list[1])


# In[ ]:


def pfileplot(l):
    for i in range(len(l)):
        

