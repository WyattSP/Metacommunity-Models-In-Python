# Main simulation loop
import numpy as np
import pandas as pd
from scipy.stats import poisson, uniform, lognormal, normal
from scipy.spatial.distance import cdist
import os
# Need to find some JIT libraries that are fast...

# For initalization might have to pull variables / simulation commands from an initalization file.
# Still use an argpaser argument but have its sole input be the inialization file.

# Function definitions

# User Inputs
# Number of patchs within the metacommunity
# Starting number of species
# Dominant metacommunity model or relative strength of such models

# Global imports
climate = NULL
distance_matrix = NULL
landscape_configuration = NULL
simulation_time = NULL
burn_in = NULL

# Starting number of species
S = 50
# Environmental dimensions from matrix
M = distance_matrix.shape[0] # This may become a matrix from the landscape input
# Environmental landscape to sca results
# Extend matrix by the time dimnsion to save results
# Need to see what this actually does... appears to store results or feed input climate per time step
sim_array = np.zeros((max(simulation_time), M)) # Probably replace env.time with time
for i in range(1, max(simulation_time)+1):
    sim_array[i-1, :] = climate[i].to_numpy()

# Set simulation times
generations = int(burn_in + simulation_time) # May need to add a round function...

# Set initial values or max values
d_species = int(round(S * 0.3)) # Number of super species ... might remove this
r = 5 # Max growth rate
z = np.repeat(np.random.rand(S), M).reshape(S, M).T # Environmental optimum of species i
a_vector = np.exp(np.linespace(np.log(1e-5), np.log(1), 16)) # Dispersal Rate
a_vector = a_vector[a_vector < 1] # Dispersal Rate
sigma_vector = np.exp(np.linspace(np.log(0.001), np.log(10), 13)) # Niche breadth/width
alpha_ij_vector = [0.5, 1.5] # Interspecific competitive effects (aij) drawn from a uniform distribution
# alpha_ii_vector = [0.5, 1.5] # Intraspecific competitive effects (aii) drawn from a uniform distribution

# The entire below for loop can be seperated into a function call for intial parameter set up
# They may iterate to explore all 4 scenarios which you can alter...
# Set up interaction coefficients for the simulation
for k in np.arange(1, len(alpha_ij_vector) + 2):
    if k <= len(alpha_ij_vector):
        if alpha_ij_vector[k - 1] == 0:
            # Equal
            # If no competitive interactions aii = aij
            alpha = np.zeros((S,S)) * 1.0
            alpha_write = str(alpha_ij_vector[k - 1])
        else:
            # Mixed
            # Competitive interaction randomly drawn from a uniform distribution
            alpha = np.random.uniform(0, alpha_ij_vector[k - 1], (S, S))
            alpha_write = str(alpha_ij_vector[k - 1])
    elif k == len(alpha_ij_vector) + 1:
        # Equal
        # aii = aij
        alpha = np.ones((S, S)) * 1.0
        alpha_write = "equal"
    else:
        # Competition-Colonization Tradeoff
        # sets values of aij including for dominant species
        # I may want to remove the "dominant" species
        alpha = np.random.uniform(0, 1, (S, S))
        alpha_hold = np.random.uniform(0, 1, (S, S))
        alpha[0:d_species, :] = np.random.unifrom(1, 1.5, (d_species, S))
        alpha[np.tril_indices(S)] = alpha_hold[np.tril_indices(S)]
        a_write = "patch_dynamics"

    np.fill_diagonal(alpha, 1.0)
    alpha *= 0.05

    # Iterating through each landscape
    for i in range(len(a_vector)):
        a = a_vector[i]

        # Iterating through each species
        for j in range(len(sigma_vector)):
            sigma_niche = sigma_vector[j]

            # Initial populations drawn from a poisson distribution
            N = np.random.poisson(0.5, (M, S)) * 1.0

            # Might want to adjust this seed setting?
            seed_V = np.arange(burn_in//10, burn_in//2+1, burn_in//10, dtype = int)
            # Actual simulation length
            samp_V = np.arange(burn_in+800, generations + 1, 20, dtype = int)
            # Save population values at each step
            N_save = N.copy()
            lambda_save = np.zeros((M, S))
            env_save = z.copy()
            den_save = np.zeros((M, S))
            env_match_save = np.zeros((M, S))

            # Simulate population dynamics for generation time
            for g in range(1, generations + 1):
                if g in seed_V:
                    N += np.random.poisson(0.5, (M, S)) * 1.0

                x = np.repeat(sim_array[g - 1, :], S).reshape(M, S)
                env = np.exp(-((x - z) / (2.0 * sigma_niche)) ** 2.0)
                density = N * alpha
                lambda_vector = r * N * (1.0 / (1.0 + density)) * env
                lambda_vector[lambda_vector < 0.0] = 0.0

                N = [np.random.poisson(l) for l in lambda_vector.flatten()]
                N = np.array(N).reshape(M, S)

                # Immigration/Emmigration
                # If immigration is possible...
                if k < len(alpha_ij_vector) + 2:
                    emigrants = [np.random.binomial(n, alpha) for n in N.flatten()]
                    emigrants = np.array(emigrants).reshape(M, S)
                    immigrants_exp = distance_matrix @ emigrants
                    immigrants_S = np.sum(emigrants, axis = 0)
                    immigrants = np.zeros((M, S))
                    for l in range(S):
                        immigrants[:, l] = np.random.poisson(immigrants_exp[:, l])
                    N += immigrants
                else:
                    N += np.random.poisson(distance_matrix @ (N * alpha), (M, S))

                N[N < 0] = 0

                # Set up save pipeline

                # Write to file once function is complete

# https://docs.python.org/3/library/multiprocessing.shared_memory.html
# Using shared memory for multiple threads might be the best idea
