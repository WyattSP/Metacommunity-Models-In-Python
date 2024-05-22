#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 10:32:06 2024

@author: wyattpetryshen
"""

# Library imports
import numpy as np


def initalize_aij(end_member, unif_start, unif_end, S, percent_dominant):
    alpha_ij_vector = [unif_start, unif_end]
    d_species = int(round(S * percent_dominant))
    
    if str(end_member) == "equal":
        # Array of Ones
        alpha = np.ones((S,S)) * 1.0
        
    elif str(end_member) == "stabalizing":
        # Array of Low Values Non-Zero's
        alpha = np.random.uniform(alpha_ij_vector[0], alpha_ij_vector[1], (S, S))
        # Set aii
        np.fill_diagonal(alpha, 1.0)
        
    elif str(end_member) == "mixed":
        # Array
        alpha = np.random.uniform(alpha_ij_vector[0], alpha_ij_vector[1], (S, S)) 
        # Set aii
        np.fill_diagonal(alpha, 1.0)
        
    elif str(end_member) == "tradeoff":
        # Array Higher Value Non-Zero's
        alpha = np.random.uniform(0, 1, (S, S))
        alpha_hold = np.random.uniform(0, 1, (S, S))
        # Randomly create dominant species
        alpha[0:d_species, :] = np.random.uniform(1, 1.5, (d_species, S))
        # Set lower traingle of matrix to be lower competitorss 
        alpha[np.tril_indices(S)] = alpha_hold[np.tril_indices(S)]
        # Set aii
        np.fill_diagonal(alpha, 1.0)
    
    else:
        return(print("Please Specify Aij"))
    
    # Multiply all values by 0.05 to increase EQ abundances
    alpha *= 0.05
    
    return(alpha)


def initalize_aii(aii_vec):
    return(aii_vec)


def random_dispersal_vector(vector_length):
    # dispersal_rate = 1 No dispersal resistance; dispersal rate approaching 0 High dispersal resistence
    # Initial configuration gives equal dispersal ability for all s
    # Equally sample along a line 
    dispersal_rate = np.exp(np.linspace(np.log(1e-5), np.log(1), vector_length))
    # Only allow values less than 1
    dispersal_rate = dispersal_rate[dispersal_rate < 1]
    
    return(dispersal_rate)


def random_sigma_vector(vector_length):
    # Equally sample along a line 
    sigma_vector = np.exp(np.linspace(np.log(0.001), np.log(10), vector_length))
    
    return(sigma_vector)


# Need to check that this function works 
def value_per_species(trait_value, species_number, unique_vales):
    # Sample for multiple sigma or dispersal
    ret_trait_vec = np.repeat(np.random.choice(trait_value, size = unique_vales), species_number)
    
    return(ret_trait_vec)

def get_dispersal_probability(dispersal_dist_k, decay_constant):
    # Get probabilities
    temp_probs = np.exp(- dispersal_dist_k * decay_constant)
    # Set values equal to 1 to zero
    temp_probs[temp_probs == 1] = 0
    # Normalize array to sum to 1
    return_probs = temp_probs / sum(temp_probs)
    
    return(return_probs)
    
def get_emigration(N, dispersal_rate, alpha_matrix):
    # Find species that emigrate
    emigrants = [np.random.binomial(n, dispersal_rate) for n in N.flatten()]
    # Add additional emigrant manipulation for competition-colonization tradeoff 
    # Make this easy to remove; Just replace values with new values 
    dom_index = alpha_matrix > 0.05
    dom_index = [True in dom_index[i,:] for i in np.arange(0,dom_index.shape[0])]
    fac = [1 if a == False else 0.01 for a in dom_index]

    # Find Emigration
    emigrants = np.zeros(N.shape)
    for i in range(N.shape[0]): # M patchs
        for j in range(N.shape[1]): # S species
            emigrants[i,j] = np.random.binomial(N[i,j], dispersal_rate * fac[j])
            
    return(emigrants)

def get_immigration(M, S, distance_matrix, emigrants):
    immigrants = np.zeros((M, S))
    for k in range(immigrants.shape[0]):
        # find index of emigrants
        index_emig = np.where(emigrants[k,:] > 0)[0]
        # assign immigration
        for i in index_emig:
            # remove emigration patch from possible immigration patchs
            # potential_patch = np.delete(np.arange(0, M), k)
            dispersal_dist_k = distance_matrix[k,:]
            # Get dispersal probabilities based on exponential function
            dispersal_probability = get_dispersal_probability(dispersal_dist_k, 0.05)
            # Disperse based on distance
            draws = int(np.sum(emigrants[k,:], axis = 0))
            #np.random.binomial(np.arange(0,M), dispersal_probability)
            img = np.random.choice(np.arange(0, M), draws, p = dispersal_probability)
            # Assign immigrants to correct location in immigrants array
            img_ind, img_count = np.unique(img, return_counts = True)
            immigrants[k, img_ind] = img_count
            
    return(immigrants)







