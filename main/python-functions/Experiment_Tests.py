#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:49:00 2024

@author: wyattpetryshen
"""

# Import functions
from MCME import MCME
import MCME_Functions
import numpy as np
from MCME_Plot_Functions import mcme_plot

seed_time = 200
burn_in = 800 # Burn in time
simulation_time = 1200 # Simulation time
seed_step = 10
    
species_number_ini = 50 # Starting number of species
max_growth_r = 5 # r max
M = 5

# Set interaction matrix
alpha_matrix_ini = MCME_Functions.initalize_aij("stabalizing", 0.0, 1.0, species_number_ini, 0.3)

# Distance matrix
distance_matrix_ini = np.ones([M,M])
np.fill_diagonal(distance_matrix_ini, 0.0)

# Provide Climate
climate_input = (np.random.rand(M*simulation_time) * 1.0).reshape(M,simulation_time) 

# Provide initial niche optimum
niche_optimum_ini = np.repeat(np.random.rand(species_number_ini), M).reshape(species_number_ini, M).T 

# Provide Initial dispersal ability
# dispersal_rate_ini = np.random.choice(MCME_Functions.random_dispersal_vector(16))
dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16)[10]

# Provide Initial niche breadth
# niche_breadth_ini = np.random.choice(MCME_Functions.random_sigma_vector(13))
niche_breadth_ini = MCME_Functions.random_sigma_vector(13)[10]

N_save, lambda_save, env_save, den_save, niche_opt, alpha_matrix = MCME(seed_time, 
                                                                        burn_in, 
                                                                        simulation_time, 
                                                                        seed_step, 
                                                                        species_number_ini, 
                                                                        max_growth_r, 
                                                                        alpha_matrix_ini,
                                                                        distance_matrix_ini, 
                                                                        climate_input, 
                                                                        niche_optimum_ini, 
                                                                        dispersal_rate_ini, 
                                                                        niche_breadth_ini)

mcme_plot(N_save, [0,1,2,3,4])
