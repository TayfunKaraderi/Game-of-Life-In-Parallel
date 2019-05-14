# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:20:25 2019

@author: karad
"""

import numpy as np
import matplotlib.pyplot as plt

# depending on the number of files change the below parameters
#-----
it=4 # number of iterations
nrows=3 # number of rows 
ncols=3 # number of columns (e.g. 4 processes -- 2row by 2column)
#-----

a = []
for i in range(0, it):
    for r in range(nrows):
        for c in range(ncols):
             globals()['data_%i_%d_%d' % (i,r,c)] = np.loadtxt('output_iteration_%d_processrow_%d_processcolumn_%d.txt'%(i,r,c),bool)
             a.append(np.loadtxt('output_iteration_%d_processrow_%d_processcolumn_%d.txt'%(i,r,c),bool))
d=[]          
c=[]
b = []
for iteration in range(it):
    for row in range(nrows):
        for col in range(ncols):
            b.append(a[iteration*nrows*ncols + row*ncols +col])
        c.append(b)
        b=[]
    d.append(c)
    c=[]



#row1 = np.concatenate(([a[0], a[1]]),axis=1) #column wise
f = []
e = []
for t in range(it):
    for col in range(nrows):
         globals()['row_%d_'%col] = np.concatenate(d[t][col],axis=1) # column wise merge
         e.append(np.concatenate(d[t][col],axis=1)) #column wise merge
    f.append(e)
    e = []

for t in range(it):
    matrix = np.concatenate((f[t]),axis=0) #row wise merge
    plt.imshow(matrix)
    plt.savefig('iteration_%d.jpg'%t)
    np.savetxt('iteration_%d.txt'%t,matrix,fmt='%d')
