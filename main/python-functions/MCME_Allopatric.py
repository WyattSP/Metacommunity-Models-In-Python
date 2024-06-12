#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:04:08 2024

@author: wyattpetryshen
"""

# Main simulation loop
import numpy as np
#from scipy.stats import poisson, uniform, lognormal, normal
#from scipy.spatial.distance import cdist

import MCME_Functions
# Custom imports
#import importlib.machinery
#MCME_Functions = importlib.machinery.SourceFileLoader('MCME_Functions','/home/wp288/project/MCME/python-functions/MCME_Functions.py').load_module()

def MCME(time_in,
         species_number_ini,
         max_growth_r,
         max_population,
         speciation_rate,
         alpha_matrix_ini,
         distance_matrix,
         climate_input,
         niche_optimum_ini,
         dispersal_rate_ini,
         niche_breadth_ini,
         end_member):

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
    n_max = max_population

    # Import landscape
    #distance_matrix = landscape_configuration
    distance_matrix = distance_matrix
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
    ####### Can set save interval per simulation here
    samp_V = np.arange(0, simulation_time + 1, 1, dtype = int)

    # Intialize save lists
    N_save = list()
    lambda_save = list()
    env_save = list()
    den_save = list()
    niche_opt = list()

    # Matrix storing evolutionary history
    ancestory = list()
    divergence_time = list()

    # Step Counter
    counter = 0 # Remomve for Cluster

    #########################
    #     Seeding Model     #
    #########################

    # Seed local patches randomly every 10 steps
    # No chance for speciation or evolution during seeding
    if seed_time > 0:
        # Set initial static environment and growth rate
        x_ini = np.repeat(climate[0, :], S).reshape(M, S)
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
            # Set max N
            #N[N > n_max] = n_max
            # Emigration
            emigrants = MCME_Functions.get_emigration(N, dispersal_rate, alpha_matrix)
            # Immigration
            immigrants, trait_flow = MCME_Functions.get_immigration(M, S, distance_matrix, emigrants, 0.00025, z)
            # Substract emigrants and add immigrants
            N -= emigrants
            N += immigrants
            # Trait frequency homogenization
            z = MCME_Functions.trait_frequnecy_homonization(z, N, trait_flow)

            counter += 1 # Remomve for Cluster
        print("Seed: %d" % counter) # Remomve for Cluster

    #########################
    #     Burn in Model     #
    #########################

    if burn_in > 0:
        # Burn in loop
        for b in range(0, burn_in):
            # Find r
            r_burnin = r * np.exp(-((z - x_ini) / (2.0 * sigma_niche)) ** 2.0)
            # Population Dynamics with static Env
            density = N @ alpha_matrix
            # Competition Model
            lambda_vector = r_burnin * N * (1.0 / (1.0 + density))
            lambda_vector[lambda_vector < 0.0] = 0.0
            # Random Death
            N = [np.random.poisson(l) for l in lambda_vector.flatten()]
            N = np.array(N).reshape(M, S) * 1.0
            # Set max N
            #N[N > n_max] = n_max
            # Emigration
            emigrants = MCME_Functions.get_emigration(N, dispersal_rate, alpha_matrix)
            # Immigration
            immigrants, trait_flow = MCME_Functions.get_immigration(M, S, distance_matrix, emigrants, 0.00025, z)
            # Substract emigrants and add immigrants
            N -= emigrants
            N += immigrants
            # Trait frequency homogenization
            z = MCME_Functions.trait_frequnecy_homonization(z, N, trait_flow)
            # Random trait evolution of Z
            ### May want to reduce the amount for a random trait to evolve
            z += MCME_Functions.evolve_trait(z, S, 0.01)
            
            counter += 1 # Remomve for Cluster
        print("Burn In: %d" % counter) # Remomve for Cluster
    # Save N after burn in
    N_save.append(N)
    niche_opt.append(z)

    #########################
    #    Main Simulation    #
    #########################
    # Currently no real function calls that redefine M and S; shouldn't be a problem but could pose issues for adding M's

    # Simulate population dynamics for generation time
    for g in range(0, simulation_time):
        ###################
        # Abiotic Forcing #
        ###################

        # Extract climate for each species in patches
        x = np.repeat(climate[g, :], S).reshape(M, S)

        #####################
        #     Selection     #
        #####################
        # i.e Ecology

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

        # Set metacommunity trait
        # Note that the order of selection on the trait will be significant.
        # if you want to add a trade-off do so here...

        ###################
        #    Dispersal    #
        ###################

        # Get emigrants
        emigrants = MCME_Functions.get_emigration(N, dispersal_rate, alpha_matrix)

        # Get immigrants
        # This needs to return an identical matrix for trait values
        immigrants, trait_flow = MCME_Functions.get_immigration(M, S, distance_matrix, emigrants, 0.00025, z)

        # Substract emigrants and add immigrants
        N -= emigrants
        N += immigrants

        # Trait frequency homogenization
        z = MCME_Functions.trait_frequnecy_homonization(z, N, trait_flow)

        # For population values below 0 set to 0
        N[N < 0] = 0
        
        N[N > n_max] = n_max

        # Set up save pipeline
        if g in samp_V:
            # You are only ever iterating N(t) during simulation
            N_save.append(N) # Save populations
            lambda_save.append(lambda_vector) # Save Lambda output
            env_save.append(x) # Save environment at time-step
            den_save.append(distance_matrix)  # Save distance between patchs
            niche_opt.append(z) # Save niche optimum

        ###################
        #    Speciation   #
        ###################


        ###################
        #    Sympatric   #
        ###################        
        ancestor_species = MCME_Functions.get_speciation(N, speciation_rate)

        ###################
        #    Allopatric   #
        ###################  
        #ancestor_species = MCME_Functions.get_allopatric_speciation(N, z)

        # Save Ancestory
        if len(ancestor_species) > 0:
            anc, dec = MCME_Functions.save_ancestory(N, ancestor_species, g)
            ancestory.append(anc); divergence_time.append(dec)

            # Adding new species to N and updating interaction coefficents after generation save
            # Add new species
            N = MCME_Functions.add_new_species(N, anc, g)
            N[N < 0] = 0

            # Update interaction coefficients
            z, alpha_matrix = MCME_Functions.update_interactions(z, alpha_matrix, anc, end_member)

            # Redefine S
            S = N.shape[1]

        ###################
        #   Evolutuion    #
        ###################
        # Random trait evolution of Z
        New_z =  MCME_Functions.evolve_trait(z, S, 0.01)
        z += New_z

        ###################
        #      Save       #
        ###################
        counter += 1 # Remomve for Cluster
        # Set up save pipeline
        if g in samp_V:
            # You are only ever iterating N(t) during simulation
            N_save.append(N) # Save populations
            lambda_save.append(lambda_vector) # Save Lambda output
            env_save.append(x) # Save environment at time-step
            den_save.append(distance_matrix)  # Save distance between patchs
            niche_opt.append(z) # Save niche optimum


    #print("Simulation End: %d" % counter) # Remomve for Cluster
    print("Total Simulation Generations: %d" % total_generation_time) # Remomve for Cluster
    # Write to file once function is complete
    # https://docs.python.org/3/library/multiprocessing.shared_memory.html
    # Using shared memory for multiple threads might be the best idea
    # Removed den_save from output
    # Removed lambda_save from output
    # Removed env_save from outuput
    return(N_save, niche_opt, alpha_matrix, ancestory, divergence_time)
    #return() # Return for profiling; Delete for complete model
