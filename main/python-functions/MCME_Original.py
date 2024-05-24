#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:03:28 2024

@author: wyattpetryshen
"""

# Main simulation loop
import numpy as np
#from scipy.stats import poisson, uniform, lognormal, normal
#from scipy.spatial.distance import cdist
import os

# Custom imports
import MCME_Functions

def MCME(time_in,
         species_number_ini,
         max_growth_r,
         alpha_matrix_ini,
         distance_matrix_ini,
         climate_input,
         niche_optimum_ini,
         dispersal_rate_ini,
         niche_breadth_ini):

    #########################
    #     Function inputs   #
    #########################

    # Import model parameters
    seed_time = time_in[0]
    seed_step = time_in[1]
    burn_in = time_in[2] # Burn in time
    simulation_time = time_in[3] # Simulation time

    S = species_number_ini # Starting number of species
    r = max_growth_r # r max

    # Import landscape
    #distance_matrix = landscape_configuration
    distance_matrix = distance_matrix_ini
    M = distance_matrix.shape[0]

    # Import climate
    # Climate dimensions should equal habitat patch number x simulation length (not including burn in)
    climate = climate_input

    z = niche_optimum_ini
    # Environmental optimum of species i
    dispersal_rate = dispersal_rate_ini # Initial dispersal ability
    sigma_vector = niche_breadth_ini # Initial niche breadth

    #########################
    # Intialize model start #
    #########################

    # Set alpha values
    alpha_matrix = alpha_matrix_ini

    # Set simulation times
    total_generation_time = int(seed_time + burn_in + simulation_time) # May need to add a round function...

    # Hashed out the bottom rows. Don't think I need this. Provide correct dim array to start
    # Environmental dimensions from matrix
    #sim_array = np.zeros((simulation_time, M))
    #for i in np.arange(0, simulation_time):
    #    # Adding climate data into an array that can be later referenced
    #    sim_array[i, :] = climate[i].to_numpy()

    sigma_niche = sigma_vector

    # Initial populations drawn from a poisson distribution
    N = np.random.poisson(0.5, (M, S)) * 1.0

    # Save interval
    samp_V = np.arange(0, simulation_time + 1, 20, dtype = int)

    # Intialize save lists
    N_save = list()
    lambda_save = list()
    env_save = list()
    den_save = list()
    niche_opt = list()

    # Step Counter
    counter = 0 # Remomve for Cluster

    #########################
    #     Seeding Model     #
    #########################
    # Seed local patches randomly every 10 steps
    if seed_time > 0:
        # Set initial static environment and growth rate
        x_ini = np.repeat(climate[:, 0], S).reshape(M, S)
        r_ini = r * np.exp(-((z - x_ini) / (2.0 * sigma_niche)) ** 2.0)

        seed_V = np.arange(0 + seed_step - 1, seed_time + seed_step - 1, seed_step, dtype = int)
        # Seed in loop
        for t in range(1, seed_time + 1):
            # Seed if t in seed_time
            if t in seed_V:
                N += np.random.poisson(0.5, (M, S)) * 1.0
            # Population Dynamics with static Env
            density = N @ alpha_matrix
            # Competition Model
            lambda_vector = r_ini * N * (1.0 / (1.0 + density))
            lambda_vector[lambda_vector < 0.0] = 0.0
            # Random Death
            N = [np.random.poisson(l) for l in lambda_vector.flatten()]
            N = np.array(N).reshape(M, S) * 1.0
            # Emigration
            emigrants = MCME_Functions.get_emigration(N, dispersal_rate, alpha_matrix)
            # Immigration
            immigrants = MCME_Functions.get_immigration(M, S, distance_matrix, emigrants)
            # Substract emigrants and add immigrants
            N -= emigrants
            N += immigrants

            counter += 1 # Remomve for Cluster
        print("Seed: %d" % counter) # Remomve for Cluster
    #########################
    #     Burn in Model     #
    #########################
    if burn_in > 0:
        # Burn in loop
        for b in range(0, burn_in):
            # Population Dynamics with static Env
            density = N @ alpha_matrix
            # Competition Model
            lambda_vector = r_ini * N * (1.0 / (1.0 + density))
            lambda_vector[lambda_vector < 0.0] = 0.0
            # Random Death
            N = [np.random.poisson(l) for l in lambda_vector.flatten()]
            N = np.array(N).reshape(M, S) * 1.0
            # Emigration
            emigrants = MCME_Functions.get_emigration(N, dispersal_rate, alpha_matrix)
            # Immigration
            immigrants = MCME_Functions.get_immigration(M, S, distance_matrix, emigrants)
            # Substract emigrants and add immigrants
            N -= emigrants
            N += immigrants

            counter += 1 # Remomve for Cluster
        print("Burn In: %d" % counter) # Remomve for Cluster
    # Save N after burn in
    N_save.append(N)
    niche_opt.append(z)

    #########################
    #    Main Simulation    #
    #########################

    # Simulate population dynamics for generation time
    for g in range(0, simulation_time):

        # Extract climate for each species in patches
        x = np.repeat(climate[:, g], S).reshape(M, S)

        # Density-independent Growth
        # z = species environmental optimum; x = environment
        r_ix = r * np.exp(-((z - x) / (2.0 * sigma_niche)) ** 2.0)

        # Competition
        density = N @ alpha_matrix

        # Species density change
        lambda_vector = r_ix * N * (1.0 / (1.0 + density))
        lambda_vector[lambda_vector < 0.0] = 0.0

        # Random value pulled from poisson distribution
        N = [np.random.poisson(l) for l in lambda_vector.flatten()]
        N = np.array(N).reshape(M, S) * 1.0

        # Get emigrants
        emigrants = MCME_Functions.get_emigration(N, dispersal_rate, alpha_matrix)

        # Get immigrants
        immigrants = MCME_Functions.get_immigration(M, S, distance_matrix, emigrants)

        # Substract emigrants and add immigrants
        N -= emigrants
        N += immigrants

        # For population values below 0 set to 0
        N[N < 0] = 0

        # Set up save pipeline
        if g in samp_V:
            # You are only ever iterating N(t) during simulation
            N_save.append(N) # Save populations
            lambda_save.append(lambda_vector) # Save Lambda output
            env_save.append(x) # Save environment at time-step
            den_save.append(distance_matrix)  # Save distance between patchs
            niche_opt.append(z) # Save niche optimum

        counter += 1 # Remomve for Cluster

    print("Simulation End: %d" % counter) # Remomve for Cluster
    print("Total Simulation Generations: %d" % total_generation_time)
    # Write to file once function is complete
    # https://docs.python.org/3/library/multiprocessing.shared_memory.html
    # Using shared memory for multiple threads might be the best idea
    return(N_save, lambda_save, env_save, den_save, niche_opt, alpha_matrix)
