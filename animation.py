# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 08:54:58 2019

@author: karad
"""

import imageio

N_files = 1000 #number of files to read
filenames = [];
for i in range(N_files):
    filenames.append("iteration_%d.jpg"%i)

with imageio.get_writer('animation.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)