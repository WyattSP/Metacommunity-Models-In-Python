#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:41:23 2024

@author: wyattpetryshen
"""
import numpy as np
import matplotlib.pyplot as plt

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

def get_hill_number(patch_occupancy, q):
    # Reference to Hill Number Calculation from:
    #   Chao et al. 2014. doi: 10.1146/annurev-ecolsys-120213-091540
    # input must be an numpy array or matrix
    # get relative frequncies / normalized spcies richness
    # drop zero values and get sum
    patch_occupancy = patch_occupancy[patch_occupancy != 0]
    species_freq = patch_occupancy / sum(patch_occupancy)
    # check sum equals 1
    if round(sum(species_freq),5) != 1.0:
        return(print("Error in Species Frequnecy (Sum not equal to 1)"))
    # if q == 1 need to take limit 
    if q == 1:
        # for all other values of q
        hill_number = np.exp(-sum(species_freq * np.log(species_freq)))
    else:
        hill_number = sum(species_freq ** q) ** (1/(1-q))
    return(hill_number)

def compare_patch_by_hill(species_patch_matrix, q_sequence, species_axis, plot = True):
    if species_axis == 0:
        # transpose matrix so species are by column
        species_patch_matrix = np.transpose(species_patch_matrix)
        patch_index = 0 # species values are by rows
    elif species_axis == 1:
        patch_index = 0  # species values are by columns
    else:
        return(print("Error: species axis must be provided"))
    # create matrix to store values based on q_sequence and n_species
    n_patches =  species_patch_matrix.shape[patch_index]
    hill_matrix = np.empty((n_patches,len(q_sequence))) 
    # loop through q_sequence and patchs 
    for i, q_val in enumerate(q_sequence):
        for p in np.arange(n_patches):
            hill_matrix[p, i] = get_hill_number(species_patch_matrix[p,:], q_val)
    if plot == True:
        fig, axs = plt.subplots(1, 1, figsize=(12, 10))
        for rows in np.arange(hill_matrix.shape[0]):
            axs.scatter(q_sequence, hill_matrix[rows,:])
            axs.plot(q_sequence, hill_matrix[rows,:])
        axs.set_xlim([0,np.max(q_sequence)])
        axs.set_xlabel("Order q",  size = 15)
        axs.set_ylabel("Hill Number",  size = 15)
        return(hill_matrix)
    else:
        return(hill_matrix)
    
    
# Test of Hill Numbers
#test_patch_matrix = np.array([[10,20,0,5,6,2,5,10,11,33,99],[5,30,2,2,6,2,5,9,22,4,1]])
#compare_patch_by_hill(test_patch_matrix, np.arange(0,3.5,0.5), 1, False)
 
        
        
        
        
        