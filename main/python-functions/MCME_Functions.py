#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 10:32:06 2024

@author: wyattpetryshen
"""

# Library imports
import numpy as np
import scipy.sparse as sparse
from scipy.spatial.distance import pdist, squareform
import geopy.distance


############################
# Initialization Functions #
############################

def initalize_aij(end_member, unif_start, unif_end, S, percent_dominant):
    alpha_ij_vector = [unif_start, unif_end]
    d_species = int(round(S * percent_dominant))
    
    if str(end_member) == "equal":
        # Array of Ones
        alpha = np.ones((S,S)) * 1.0
        
    elif str(end_member) == "neutral":
        # Array of zeros
        alpha = np.zeros((S,S))
        
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
    # In the original Thompson paper this is toggled on
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

def initial_species_niche_optimum(S, M, initial_climate_slice):
    # This function produces random niche optimum for starting number of species
    # The niche optimum is constant for ALL habitat patches
    # The random value, z, is chosen from a standard deviation of the current time step
    # Could the sd be correlated to previous climate variability???
    
    # Get mean of first t-slice 
    time_zero_mean = np.mean(initial_climate_slice)
    sd_all = np.std(initial_climate_slice)
    return_vals = np.random.normal(time_zero_mean, sd_all, S)
    
    # Need to return a M by S matrix of starting niche optimums 
    niche_z = np.zeros(shape = (M,S))
    niche_z[:,:] = return_vals

    return(niche_z)

############################
#   Dispersal Functions    #
############################

def get_dispersal_probability(dispersal_dist_k, decay_constant):
    # Get probabilities
    # If you want to add other forms of dispersal resistence change dispersal_dist_k to dispersal_resist_k
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

def get_immigration(M, S, distance_matrix, emigrants, resistance):
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
            dispersal_probability = get_dispersal_probability(dispersal_dist_k, resistance)
            # Disperse based on distance
            draws = int(np.sum(emigrants[k,:], axis = 0))
            #np.random.binomial(np.arange(0,M), dispersal_probability)
            img = np.random.choice(np.arange(0, M), draws, p = dispersal_probability)
            # Assign immigrants to correct location in immigrants array
            img_ind, img_count = np.unique(img, return_counts = True)
            immigrants[img_ind, i] += img_count
            
    return(immigrants)


############################
#  Speciation Functions    #
############################

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
    # 
    spec_in = np.insert(arr = ancestors, obj = 0, values = -999)
    prob_in = np.insert(arr = probs, obj = 0, values = probs_no)
    # draw new species
    new_species = np.random.default_rng().choice(a = spec_in, size = len(spec_in) - 1, p = prob_in, replace = True)
    # return array of anscetor species for new species
    anscestor_species = new_species[new_species != -999]
    return(anscestor_species)


def save_ancestory(N, ancestor_species, current_simulation_time_step):
    n_anc = len(ancestor_species)
    # Get current number of species
    current_species_number = N.shape[1]
    # Get ID's of new species
    new_species_ids = list(range(current_species_number, current_species_number + len(ancestor_species)))
    # Stack ancestor with its corresponding descendent
    anc_out = np.vstack((ancestor_species,new_species_ids))
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
        if sample_in.size == 0:
            print(sample_in, n_species)
        else:
            # for point mutation add 1 unit individual for descendent and remove 1 for ancestor
            patch_index = np.random.default_rng().choice(a = sample_in, size = 1)[0]
            # add new species
            new_N[patch_index, in_species[1,i]] += 1
            # subtract old species
            new_N[patch_index, in_species[0,i]] -= 1
    # Return new N matrix
    return(new_N)


def random_endmember_interactions(end_member, iter_nan, min_aij, max_aij):

    if str(end_member) == "equal":
        # vector of all ones
        new_aij = [0.05 for i in range(len(iter_nan))]
        
    elif str(end_member) == "neutral":
        # vector of all zeros
        new_aij = [0.0 for i in range(len(iter_nan))]
        
    elif str(end_member) == "stabalizing":
        # Array of values in range
        new_aij = [np.random.uniform(min_aij, max_aij) for i in range(len(iter_nan))]
        
    elif str(end_member) == "mixed":
        # Array of values in range
        new_aij = [np.random.uniform(min_aij, max_aij) for i in range(len(iter_nan))]        
        
    return(new_aij)
    

def update_interactions(z, alpha_matrix, anc, end_member):
    # z = niche optimums of previous time step species
    # alpha_matrix = species interactions
    # anc = phylogeny of new species (ancestor, descendent)
    
    # Determine new species ancestory
    in_species = anc
    # Old number of species from niche optimum vector
    old_n = z.shape[1]
    # New number of species; doesn't actually matter which row you look at
    n_species = len(in_species[0,:]) # look at all possible columns
    # Total number of species for new Z array
    total_S = old_n + n_species
    # Get new species optimum for new species 
    # Create new z vector (1d-array)
    #new_z = np.copy(z)
    new_z = np.zeros((z.shape[0], total_S))
    # Fill in with existing z
    new_z[0:z.shape[0],0:z.shape[1]] = z
    
    # Fill in new z's
    # Get ancestor z value
    for i in range(n_species):
        # index old species
        old_j = in_species[0,i]
        # index new species
        new_j = in_species[1,i]
        # old z value for species i ancestor
        # in a 1 species system first value will be 0 for the index
        anc_z_val = z[:,old_j]
        # add random jitter from uniform 
        # append new value to new_z
        new_z[:,new_j] = anc_z_val + np.random.uniform(-0.05, 0.05, 1)
        
    if new_z.shape[1] != total_S:
        print("WARNING: Error in New Z Calculation")
    
    # Create new interaction matrix (alpha matrix)
    new_alpha_matrix = np.empty((old_n + n_species , old_n + n_species))
    new_alpha_matrix[:] = np.nan
    # Copy existing alpha matrix into new matrix
    new_alpha_matrix[0:alpha_matrix.shape[0], 0:alpha_matrix.shape[1]] = alpha_matrix
    # Set diagonal
    np.fill_diagonal(new_alpha_matrix, 0.05)
    
    # Set new interactions
    iter_nan = np.argwhere(np.isnan(new_alpha_matrix))
    # Get max and min interaction values from orginal interaction matrix
    max_aif = np.amax(alpha_matrix)
    min_aif = np.amin(alpha_matrix)
    # Get new interaction values... note that their are no tradeoffs here...
    new_aijs = random_endmember_interactions(end_member, iter_nan, min_aif, max_aif)
    # Fill in interaction matrix
    for i, val in enumerate(iter_nan):
        # row
        row = val[0]; col = val[1]
        new_alpha_matrix[row, col] = new_aijs[i]
    
    # return output
    return(new_z, new_alpha_matrix)
    
############################
#    Climate Functions     #
############################


def random_coordinates(N_coords, climate_type = "HadCM3"):
    if climate_type == "HadCM3":
        lat_centers = [90.0000,87.5000,85.0000,82.5000,80.0000,77.5000,75.0000,72.5000,
                     70.0000,67.5000,65.0000,62.5000,60.0000,57.5000,55.0000,52.5000,
                     50.0000,47.5000,45.0000,42.5000,40.0000,37.5000,35.0000,32.5000,
                     30.0000,27.5000,25.0000,22.5000,20.0000,17.5000,15.0000,12.5000,
                     10.0000,7.5000,5.0000,2.5000,0.0000,-2.5000,-5.0000,-7.5000,
                     -10.0000,-12.5000,-15.0000,-17.5000,-20.0000,-22.5000,-25.0000,-27.5000,
                     -30.0000,-32.5000,-35.0000,-37.5000,-40.0000,-42.5000,-45.0000,-47.5000,
                     -50.0000,-52.5000,-55.0000,-57.5000,-60.0000,-62.5000,-65.0000,-67.5000,
                     -70.0000,-72.5000,-75.0000,-77.5000,-80.0000,-82.5000,-85.0000,-87.5000,
                     -90.0000]
        lon_centers = [0.0000, 3.7500,7.5000,11.2500,15.0000,18.7500,22.5000,26.2500,
                     30.0000,33.7500,37.5000,41.2500,45.0000,48.7500,52.5000,56.2500,
                     60.0000,63.7500,67.5000,71.2500,75.0000,78.7500,82.5000,86.2500,
                     90.0000,93.7500,97.5000,101.2500,105.0000,108.7500,112.5000,116.2500,
                     120.0000,123.7500,127.5000,131.2500,135.0000,138.7500,142.5000,146.2500,
                     150.0000,153.7500,157.5000,161.2500,165.0000,168.7500,172.5000,176.2500,
                     180.0000,183.7500,187.5000,191.2500,195.0000,198.7500,202.5000,206.2500,
                     210.0000,213.7500,217.5000,221.2500,225.0000,228.7500,232.5000,236.2500,
                     240.0000,243.7500,247.5000,251.2500,255.0000,258.7500,262.5000,266.2500,
                     270.0000,273.7500,277.5000,281.2500,285.0000,288.7500,292.5000,296.2500,
                     300.0000,303.7500,307.5000,311.2500,315.0000,318.7500,322.5000,326.2500,
                     330.0000,333.7500,337.5000,341.2500,345.0000,348.7500,352.5000,356.2500]
        
        zipped_coords = zip(np.random.choice(lon_centers, N_coords),np.random.choice(lat_centers, N_coords))
    else:
        return(print("Only HadCM3 currently supported"))
    return(list(zipped_coords))


def get_graph_distance_matrix_HadCM3(coord_list):
    # pdist can calculate a variety of different distance metrics
    # found here: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
    # Get vector form of distances
    # first value in coordlist tupple should be lat, second value lon

    # Input coordinates
    # Number of points
    points = len(coord_list)
    # Create iterator
    iter_a =  np.array(np.meshgrid(range(points), range(points))).reshape(2, points * points).T
    # Create matrix to store results
    distance_vector = [geopy.distance.distance(geopy.distance.lonlat(*coord_list[i[0]]), geopy.distance.lonlat(*coord_list[i[1]])).km for i in iter_a]
    # Keep this as a vector for now, but you should only calculate the lower triangle of the matrix 
    # matrix[np.tril_indices(3, k = 1)] ... you would need the off diagonal from the above coordinate list
    # can use scipy.spatial.distance.squareform(t) to put off-diagonal back into square matrix

    # Put back into a square... might change this to a long form array later
    output_distance_matrix = np.zeros((points, points))
    for i in np.arange(0, len(iter_a)):
        row_i = iter_a[i][0]
        col_j = iter_a[i][1]
        output_distance_matrix[row_i, col_j] = distance_vector[i]

    return(output_distance_matrix)

############################
#   Evolution Functions    #
############################

def evolve_trait(z, S, sd):
    i = z.shape[0]
    # Get length of needed random numbers
    z_additions = np.random.normal(loc = 0, scale = sd, size = (i * S)).reshape(i, S)
    return(z_additions)

def trait_frequnecy_homonization(z, N):
    # Loop over species; homogenize species trait by frequnecy in each patch
    for j in range(N.shape[1]):
        N_in = N[:,j]
        
        if sum(N_in) == 0:
            N_in = [i + 1 for i in N_in]
        
        new_trait_val = sum(z[:,j][g] *  N_in[g] for g in range(len(z[:,j]))) / sum(N_in)
        z[:,j] = new_trait_val
        
    return(z)


############################
#      Save Functions      #
############################















