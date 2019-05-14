# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 04:46:52 2019

@author: karad
"""
import numpy as np

image = np.loadtxt('iteration_9.txt',bool)

count=0
for i in range(image.shape[0]):
    for j in range(image.shape[1]):
        if (image[i][j]==1):
            count +=1

proportion_of_living_grids = count/(image.shape[0]*image.shape[0])
print(proportion_of_living_grids)