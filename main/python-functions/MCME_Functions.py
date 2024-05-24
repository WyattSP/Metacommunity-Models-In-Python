#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 10:32:06 2024

@author: wyattpetryshen
"""

# Library imports
import numpy as np
import scipy.sparse as sparse


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
    for k in range(M):
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
            immigrants[img_ind, i] += img_count
            
    return(immigrants)


def get_speciation(N, speciation_rate):
    # An interesting addition would be to have a different speciation rate for each species
    # Or a speciation rate that is dependent upon environmental temperature
    # HOWEVER since speciation_rate is a function of population abundance, which is also
    # a function of environment, this is already included
    # get relative abundances of species i in metacommunity
    # column sums of species i 
    total_abundance_species = np.sum(N, axis = 0)
    # get ancestors
    ancestors = np.nonzero(total_abundance_species)[0]
    # Get rid of zeros
    total_abundance_species_nonzero = total_abundance_species[total_abundance_species != 0]
    # total sum of matrix
    total_abundance = np.sum(total_abundance_species_nonzero)
    # species occurrence frequency
    relative_species_frequency = total_abundance_species_nonzero/total_abundance  
    # probabilites
    probs = speciation_rate * relative_species_frequency
    probs_no = 1 - sum(probs)
    # potential ancestors
    # inputs for draws
    spec_in = np.insert(arr = ancestors, obj = 0, values = 0)
    prob_in = np.insert(arr = probs, obj = 0, values = probs_no)
    # draw new species
    new_species = np.random.default_rng().choice(a = spec_in, size = len(spec_in) - 1, p = prob_in, replace = True)
    # return array of anscetor species for new species
    anscestor_species = new_species[new_species != 0]
    return(anscestor_species)


def save_ancestory(N, anscestor_species, current_simulation_time_step):
    n_anc = len(anscestor_species)
    # Get current number of species
    current_species_number = N.shape[1]
    # Get ID's of new species
    new_species_ids = list(range(current_species_number, current_species_number + len(anscestor_species)))
    # Stack ancestor with its corresponding descendent
    anc_out = np.vstack((anscestor_species,new_species_ids))
    # Save divergence time
    previous_step = [(current_simulation_time_step - 1)] * n_anc
    current_step = [current_simulation_time_step] * n_anc
    dec_out = np.vstack((previous_step,current_step))
    return(anc_out, dec_out)

def add_new_species(N, anc):
    in_species = anc
    # Number of new species
    n_species = len(in_species[0,])
    # Construct new dense array
    rows = N.shape[0]
    columns = N.shape[1]
    new_N = np.zeros((rows, columns + n_species))    
    # Copy old dense array to fill top-left quadrant
    new_N[0:rows, 0:columns] = N
    # Re-distribute new species in new_N
    # Randomly select habitat patch with new species i
    # first in_species[0,] = anscestor; in_species[1,] = descendent 
    for i in range(n_species):
        # find indexes with species i and select a random patch for speciation 
        sample_in = np.nonzero(new_N[:,in_species[0,i]])[0]
        if len(sample_in) <= 0:
            print(sample_in, n_species)
        # for point mutation add 1 unit individual for descendent and remove 1 for ancestor
        patch_index = np.random.default_rng().choice(a = sample_in, size = 1)[0]
        # add new species
        new_N[patch_index, in_species[1,i]] += 1
        # subtract old species
        new_N[patch_index, in_species[0,i]] -= 1
    # Return new N matrix
    return(new_N)

def update_interactions(z, alpha_matrix, anc):
    in_species = anc
    # Old number of species
    old_n = len(z)
    # New number of species
    n_species = len(in_species[0,])
    
    # Get new species optimum for new species 
    # Create new z matrix
    new_z = np.zeros((1 , old_n + n_species))
    new_z[0, 0:old_n] = z
    # Fill in new z's
    # Get ancestor z value
    for i in range(n_species):
        # old z value for species i ancestor
        old_z_val = new_z[0,in_species[0,i]]
        # add random jitter from uniform 
        new_z[0,in_species[1,i]] = (old_z_val + np.random.uniform(-0.05, 0.05, 1)[0])
        
    # Expand Z to original 
    
    # Create new alpha matrix
    new_alpha_matrix = np.zeros((old_n + n_species , old_n + n_species))
    new_alpha_matrix[0:old_n,0:old_n] = alpha_matrix
    
    # diagonal
    diag = new_alpha_matrix[list(in_species[0,:]),list(in_species[0,:])]
    # portion of lower triangle missing bottom right
    triang = new_alpha_matrix[list(in_species[0,:]), 0:old_n]
    # rebuild off diagonals and diagonal except lower right corner
    # full diagonal
    new_alpha_matrix[list(in_species[1,:]),list(in_species[1,:])] = diag
    # lower triangle
    new_alpha_matrix[list(in_species[1,:]), 0:old_n] = triang
    # upper triangle
    new_alpha_matrix[0:old_n, list(in_species[1,:])] = np.rot90(triang)
    # add bottom right corner
    iter_a = [(x, y) for x in list(in_species[0,:]) for y in list(in_species[0,:])]
    iter_r = [(x, y) for x in list(in_species[1,:]) for y in list(in_species[1,:])]
    for i in range(len(iter_r)):
        new_alpha_matrix[iter_r[i]] = alpha_matrix[iter_a[i]]
    
    return(new_z[0], new_alpha_matrix)
    


