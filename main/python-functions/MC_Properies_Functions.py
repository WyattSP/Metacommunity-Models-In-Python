#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:41:23 2024

@author: wyattpetryshen
"""
import numpy as np


def gamma_diversity(population_list):
    iter_l = len(population_list)
    # Convert to array
    # Get time slice dimension
    t_gam = list()
    for t in range(iter_l):
        # list to matrix
        pop_mat = np.array(population_list[t])
        t_gam.append(len(set(np.where(pop_mat[:,:] > 0)[1])))
    return(t_gam)

def temporal_gamma_diversity(population_list):
    # Convert to array
    pop_arr = np.dstack(population_list)
    # Get time slice dimension
    m_slices = pop_arr.shape[0]
    m_gam = list()
    for m in range(m_slices):
        m_gam.append(len(set(np.where(pop_arr[m,:,:] > 0)[0])))
    return(m_gam)

def alpha_richness(population_list):
    iter_l = len(population_list)
    alpha_t = list()
    # Average number of species in each patch
    for t in range(iter_l):
        pop_mat = np.array(population_list[t])
        m_slices = pop_mat.shape[0]
        m_alpha = 0
        for m in range(m_slices):
            m_alpha += len(np.where(pop_mat[m,:] > 0)[0])
        alpha_t.append(m_alpha / m_slices)
    return(alpha_t)

