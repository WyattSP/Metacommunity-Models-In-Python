# Set working directory
import os
os.chdir('/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/python-functions')

# Import functions
#from MCME_Speciation import MCME
from MCME_Allopatric import MCME
import MCME_Functions
import MC_Properies_Functions
import numpy as np
from MCME_Plot_Functions import mcme_plot, species_richness_plot, plot_origination_extinction, plot_richness_map, get_richness
import MCME_Plot_Functions
import butterworth_filter
import matplotlib.pyplot as plt
import sys
import Climate_Import
import pandas as pd
import random
import math
import matplotlib

#%% Parameters
### General improvements for speed
# Convert distance matrix to off-diagonal (upper or lower triangle) only, or sparse-matrix
# Having a static-global distance matrix with a resistence matrix may be better
# Resistence matrix defined by habitat niche
### You need to work through an example with 1 starting species; can't just copy old interaction matrix's

#########################
#  Set Function Inputs  #
#########################
Tmax_4kyr = Climate_Import.import_climate_data("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Proc_Bricks/Tmax.mon.4kyr.tif", "Geotif", "rasterio")

#########################
#   Global Parameters   #
#########################

S = 25 # set sarting species number
max_growth_r = 5 # set max growth rate
M = 165 # set patch number

###################################
#  Import Continental Coordinates #
###################################

# Metacommunity on North America
Cont_Coords = pd.read_csv("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/World_Continent_Coords.csv")

# Extract only NA coords
# Keep these below lines for full dataset
mask = Cont_Coords.Cont == "South America"
Clim_array = Cont_Coords[mask].reset_index()

# Convert Data Frame into array of zip (lon,lat)
Clim_Coords = list()
for i in range(len(Clim_array)):
    Clim_Coords.append((Clim_array.Lon_d[i], Clim_array.Lat_d[i]))
    
# Randomly sample points in North America
coord_index = np.random.choice(range(len(Clim_array)), M, replace = False)
coords = [Clim_Coords[i] for i in coord_index]

# Define distance matrix
distance_matrix = MCME_Functions.get_graph_distance_matrix_HadCM3(coords)

###################################
#       Interpolate Climate       #
###################################
# Intervals: 4 kya, 50 kya, 100 kya
# Native resolution climate for desired coordinates
climate_4kya = Climate_Import.return_climate_array(coords, 400, 901, Tmax_4kyr) - 273.5

###################################
#    Set up Model Parameters      #
###################################
# Max Simulation time
sim_time = climate_4kya.shape[0]
# Simulation times
time_in = [200, 10, 100, sim_time] # seed time, seed inteval, burn-in time, simulation time
# Set intial species niche optimums; this will stay static between all simulations
niche_optimum_ini_4 = MCME_Functions.initial_species_niche_optimum(S, M, climate_4kya[0,:])
###################################
#     Iteratable Parameters       #
###################################
# Dispersal ability
#dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16, 1e-5, 1)
dispersal_rate_ini = [0.01]
# Niche Breadth
niche_breadth = [2.5]
# narrow = 0.5
# wide = 2.5
# Specify the speciation rate
speciation_threshold = [0.25]
# Species interaction end-members
end_member = ["equal", "stabalizing", "mixed", "neutral"]
# Equal species interactions
alpha_matrix_equal = MCME_Functions.initalize_aij("equal", 0.0, 1.0, S, 0.3)
# Stabalizing species interactions
alpha_matrix_stabalizing = MCME_Functions.initalize_aij("stabalizing", 0.0, 0.5, S, 0.3)
# Mixed species interactions
alpha_matrix_mixed = MCME_Functions.initalize_aij("mixed", 0.0, 1.5, S, 0.3)
# Neutral species interactions
alpha_matrix_neutral = MCME_Functions.initalize_aij("neutral", 0.0, 1.0, S, 0.3)
alpha_matrix_list = list([alpha_matrix_equal, alpha_matrix_stabalizing, alpha_matrix_mixed, alpha_matrix_neutral])

#%% Simulation
#########################
#   Run Simulation Exp  #
#########################

# 4 kya forcings wide
results, niche_opt, alpha_matrix, phylogeny, divergence_time, patch_origin = MCME(time_in =  [100, 10, 100, 250],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = climate_4kya,
                                                                   niche_optimum_ini = niche_optimum_ini_4,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 2.5,
                                                                   end_member = "stabalizing")


#%% Diagnositcs
##########################
# Calculate Comm Metrics #
##########################
#%%
n_gamma = MC_Properies_Functions.gamma_diversity(results)
n_alpha = MC_Properies_Functions.alpha_richness(results)
# Hill Numbers
# Get mean hill numbers between all patchs? 
# Let's try for the last time step
#%%
ls_16 = MC_Properies_Functions.compare_patch_by_hill(results[-1], np.arange(0,4,1), 1, True)
delta_ls_16 = MC_Properies_Functions.compare_hill_difference(results[-1], [0,2], 1, True)
sum(delta_ls_16[~np.isnan(delta_ls_16)])
#%%
# Plot Hill differences
MCME_Plot_Functions.plot_hill_difference_map(delta_ls_16, Clim_array, coord_index, False)
#########################
#          Plot         # 
#########################
#%%
# Plot of gamma and alpha diversity
x = np.arange(0,len(n_alpha))
plt.plot(x, n_gamma, label = "gamma")
plt.plot(x, n_alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")
#%%
# Plot of species richness change
orgination, extinction = plot_origination_extinction(results, True)
plt.plot(range(len(np.cumsum(orgination[:,1]))), np.cumsum(orgination[:,1]), c = "blue")
plt.plot(range(len(np.cumsum(extinction[:,1]))), np.cumsum(extinction[:,1]), c = "red")
# Plot map of geographic points with colours illustrating community richness
#%%
plot_richness_map(results[-1], Clim_array, coord_index, False)

# %%
