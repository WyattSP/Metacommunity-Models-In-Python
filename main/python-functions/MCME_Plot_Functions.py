#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 07:39:16 2024

@author: wyattpetryshen
"""

import matplotlib.pyplot as plt
import numpy as np

def mcme_plot(population_list, patch):
    # This is technically the incorrect number of generations
    if len(patch) == 1:
        # Function of selected save intervals 
        x = np.arange(0,len(population_list))
        # Stack save points into array
        pop_arr = np.dstack(population_list)
        # Plot dynamics for desired patch
        for i in np.arange(0, pop_arr.shape[1]):
            plt.plot(x, pop_arr[patch,i,:])
        plt.show()
    else:
        pop_arr = np.dstack(population_list)
        plot_rows = pop_arr.shape[0] 
        fig, axarr = plt.subplots(plot_rows, 1, sharex=True, figsize=(12,6))

        x = np.arange(0,len(population_list))
        
        count = 0
        for ax in axarr.ravel():
            for j in np.arange(0, pop_arr.shape[1]):
                ax.plot(x, pop_arr[count,j,:])
            count += 1
        fig.tight_layout()
