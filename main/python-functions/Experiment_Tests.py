#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:49:00 2024

@author: wyattpetryshen
"""
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
import matplotlib.pyplot as plt
import sys
import Climate_Import
import pandas as pd
import random
import math

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
######################################
#       Fourier Smooth Climate       #
######################################
# Smooth out variance at: 4 kya, 40 kya, 100 kya
# lfreq: 
#tseries_freq = Climate_Import.get_timeseries_frequency_range(climate_4kya[:,1])
sample_rate = 4 # t series sample rate
# double sampling 4 * 3 = 12 kya
lfreq_1 = round(1 / (16 / sample_rate), 3)
# 10 times sampling 4 * 10 = 40 kya
lfreq_2 = round(1 / (40 / sample_rate), 3)
# 25 times sampling 4 * 25 = 100 kya
lfreq_3 = round(1 / (100 / sample_rate), 3)
# Get new climates
fClimate_16kya = Climate_Import.get_fourier_smooth_climate(climate_4kya[1:], lfreq_1, 0.5)
fClimate_40kya = Climate_Import.get_fourier_smooth_climate(climate_4kya[1:], lfreq_2, 0.5)
fClimate_100kya = Climate_Import.get_fourier_smooth_climate(climate_4kya[1:], lfreq_3, 0.5)

###################################
#    Set up Model Parameters      #
###################################
# Max Simulation time
sim_time = climate_4kya.shape[0]

# Simulation times
time_in = [200, 10, 100, sim_time] # seed time, seed inteval, burn-in time, simulation time

# Set intial species niche optimums; this will stay static between all simulations
niche_optimum_ini_4 = MCME_Functions.initial_species_niche_optimum(S, M, climate_4kya[0,:])
niche_optimum_ini_100 = MCME_Functions.initial_species_niche_optimum(S, M, fClimate_100kya[0,:])

###################################
#     Iteratable Parameters       #
###################################

# Dispersal ability
#dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16, 1e-5, 1)
dispersal_rate_ini = [0.01]

# Niche Breadth
niche_breadth = [0.5]

# Specify the speciation rate
speciation_threshold = [0]

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


#########################
#   Run Simulation Exp  #
#########################

# 4 kya forcings
results, niche_opt, alpha_matrix, phylogeny, divergence_time, patch_origin = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = climate_4kya,
                                                                   niche_optimum_ini = niche_optimum_ini_4,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 0.5,
                                                                   end_member = "stabalizing")


# 100 kya forcings
results_100, niche_opt_100, alpha_matrix_100, phylogeny_100, divergence_time_100, patch_origin_100 = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = fClimate_100kya,
                                                                   niche_optimum_ini = niche_optimum_ini_100,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 0.5,
                                                                   end_member = "stabalizing")

#########################
#      Profile          #
#########################
from cProfile import Profile
from pstats import SortKey, Stats

# When profiling set the return function to nothing so you don't overflow your console

with Profile() as profile:
    print(f"{ MCME(time_in =  [100, 10, 100, 500],species_number_ini = 25,max_growth_r = 5,speciation_threshold = 0.25,alpha_matrix_ini = alpha_matrix_stabalizing,distance_matrix = distance_matrix, climate_input = fClimate_100kya,niche_optimum_ini = niche_optimum_ini_100,dispersal_rate_ini = 0.1,niche_breadth_ini = 0.5,end_member = 'stabalizing') = }")
    (
     Stats(profile)
     .strip_dirs()
     .sort_stats(SortKey.CALLS)
     .print_stats()
     )

sys.getsizeof(results)

######################################
# Notes for Performance Improvements #
######################################
import time
start = time.time()
np.random.poisson(10000,int(1e7))
end = time.time()
print(end-start)


##########################
# Calculate Comm Metrics #
##########################

gamma = MC_Properies_Functions.gamma_diversity(results_100)
alpha = MC_Properies_Functions.alpha_richness(results_100)

# t_gamma = MC_Properies_Functions.temporal_gamma_diversity(N_save)

#########################
#          Plot         # 
#########################

# Plot of gamma and alpha diversity
x = np.arange(0,len(alpha))
plt.plot(x, gamma, label = "gamma")
plt.plot(x, alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")


# Plot of species richness change
orgination, extinction = plot_origination_extinction(results, True)
plt.plot(range(len(np.cumsum(orgination[:,1]))), np.cumsum(orgination[:,1]), c = "blue")
plt.plot(range(len(np.cumsum(extinction[:,1]))), np.cumsum(extinction[:,1]), c = "red")

# Plot map of geographic points with colours illustrating community richness
plot_richness_map(results[-1], Clim_array, coord_index, False)

#np.savetxt("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/test_out_raster.csv", t, delimiter=",")

###########################
#         WBF Plots       #
###########################

def return_richness_pivot_table(species_patch_array, Climate_input, Coord_index):
    lon_in = [Climate_input.Lon[i] for i in Coord_index]
    lat_in = [Climate_input.Lat[i] for i in Coord_index]
    r_in = MCME_Plot_Functions.get_richness(species_patch_array)
    r_in = [i for i in r_in]
        
    in_data = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "R": r_in})
    in_data = in_data.pivot_table(index='Lat', columns='Lon', values='R')
    return(in_data)
 
# Get richness pivot table for last time step 
data_4 = return_richness_pivot_table(results[-1], Clim_array, coord_index)
data_100 = return_richness_pivot_table(results_100[-1], Clim_array, coord_index)

# For smoothing if you want it
methods = [None, 'none', 'nearest', 'bilinear', 'bicubic', 'spline16',
           'spline36', 'hanning', 'hamming', 'hermite', 'kaiser', 'quadric',
           'catrom', 'gaussian', 'bessel', 'mitchell', 'sinc', 'lanczos']

### 4 kya plot of spatial richness
hmap = plt.imshow(data_4.values, cmap='viridis', interpolation='nearest', vmin=0)
#plt.colorbar(hmap, ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
plt.axis('off')
plt.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/s05_4_stab_heat.png', format='png', dpi=300, transparent=True)

### 100 kya plot of spatial richness
htmap = plt.imshow(data_100.values, cmap='viridis', interpolation='nearest', vmin=0)
#plt.colorbar(htmap,ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
plt.axis('off')
plt.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/s05_100_stab_heat.png', format='png', dpi=300, transparent=True)

### Anomaly plot of 4 - 100 stabilizing
htmap = plt.imshow(data_4.values-data_100.values, cmap='RdBu', interpolation='nearest', vmin=-12, vmax=12)
plt.title("Difference Map of Species Ricnhess (Low - High)")
plt.colorbar(htmap,ticks=[-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12])
plt.axis('off')
plt.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/4-100_stab_heat.png', format='png', dpi=300, transparent=True)

####### Alpha Gamma for WBF Narrow Niche
gamma_4 = MC_Properies_Functions.gamma_diversity(results)
alpha_4 = MC_Properies_Functions.alpha_richness(results)
gamma_100 = MC_Properies_Functions.gamma_diversity(results_100)
alpha_100 = MC_Properies_Functions.alpha_richness(results_100)

fig, axs = plt.subplots(1, 1)
#fig.tight_layout(pad=0.5)
axs.plot(range(len(gamma_4)), gamma_4, color='red', linestyle='-')
axs.plot(range(len(gamma_100)), gamma_100, color='black', linestyle='-')
axs.plot(range(len(alpha_4)), alpha_4, color='red', linestyle=':')
axs.plot(range(len(alpha_100)), alpha_100, color='black', linestyle=':')
axs.grid(True, which='major', axis='both')
axs.grid(True, which='major', axis='both')
axs.xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs.xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs.set_xlabel("Thousand Years",  size = 15)
axs.tick_params(axis='x', colors='black')
axs.tick_params(axis='y', colors='black')
axs.spines['bottom'].set_color('black')
axs.spines['top'].set_color('black') 
axs.spines['right'].set_color('black')
axs.spines['left'].set_color('black')
axs.yaxis.label.set_color('black')
axs.xaxis.label.set_color('black')
axs.set_ylabel("Diversity", size = 15)
plt.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/Richness_sig05_sims.png', format='png', dpi=300, transparent=True)






hmap = plt.imshow(data_100.values, cmap='YlGn', interpolation='bilinear', vmin=0, vmax=20)
plt.colorbar(hmap, ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

hmap = plt.imshow(data_4.values - data_100.values, cmap='seismic', interpolation='none', vmin=-10, vmax=10)
plt.colorbar(hmap, ticks=[-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])

sns.set(rc={'axes.facecolor':'transparent', 'figure.facecolor':'transparent'})
sns.heatmap(data_4.values - data_100.values, annot=True, fmt="g", cmap='viridis', annot_kws={"size": 5}, cbar = True, vmin=-8, vmax=8, square=True)


#### Anomaly count
p_anomlist = list()
n_anomlist = list()
z_anomlist = list()

for i in np.arange(0,500):
    data_4 = return_richness_pivot_table(results[i], Clim_array, coord_index)
    data_100 = return_richness_pivot_table(results_100[i], Clim_array, coord_index)
    x_in = data_4.values-data_100.values
    x_out = x_in[~np.isnan(x_in)]
    m_pos = np.count_nonzero(x_out > 0)
    m_neg = np.count_nonzero(x_out < 0)
    m_zero = np.count_nonzero(x_out == 0)
    p_anomlist.append(round(m_pos/len(x_out)*100,0))
    n_anomlist.append(round(m_neg/len(x_out)*100,0))
    z_anomlist.append(round(m_zero/len(x_out)*100,0))
    
print("Postive Anomaly: ", round(m_pos/len(x_out)*100,0), "%")
print("Negative Anomaly: ", round(m_neg/len(x_out)*100,0), "%")
print("Zero Anomaly: ", round(m_zero/len(x_out)*100,0), "%")


plt.plot(range(len(p_anomlist)), p_anomlist, c = "blue")
plt.plot(range(len(p_anomlist)), n_anomlist, c = "red")
plt.plot(range(len(p_anomlist)), z_anomlist, c = "grey")

#### Where are species originating / going extinct?

all_origins = [x[0] for xs in patch_origin for x in xs]
values, counts = np.unique(all_origins, return_counts=True)

all_origins_100 = [x[0] for xs in patch_origin_100 for x in xs]
values_100, counts_100 = np.unique(all_origins_100, return_counts=True)

# Order counts by the 
save_patch_index = np.zeros(coord_index.shape)
for i in range(len(values)):
    for j in range(len(coord_index)):
        if values[i] == coord_index[j]:
           save_patch_index[j] = counts[i]
        else:
            continue
save_patch_index = [i for i in save_patch_index]

patch_in_data = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "P": save_patch_index})
patch_in_data = patch_in_data.pivot_table(index='Lat', columns='Lon', values='P')

hmap = plt.imshow(patch_in_data, cmap='viridis', interpolation='none')
plt.colorbar(hmap, ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])


##### 100 kya
# Order counts by the 
save_patch_index_100 = np.zeros(coord_index.shape)
for i in range(len(values_100)):
    for j in range(len(coord_index)):
        if values[i] == coord_index[j]:
           save_patch_index_100[j] = counts_100[i]
        else:
            continue
save_patch_index_100 = [i for i in save_patch_index_100]

patch_in_data_100 = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "P": save_patch_index_100})
patch_in_data_100 = patch_in_data_100.pivot_table(index='Lat', columns='Lon', values = "P")

hmap = plt.imshow(patch_in_data_100, cmap='viridis', interpolation='none')
plt.colorbar(hmap, ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])



t_patch = patch_in_data - patch_in_data_100
np.count_nonzero(t_patch > 0)
np.count_nonzero(t_patch < 0)
np.count_nonzero(t_patch == 0)
plt.imshow(patch_in_data - patch_in_data_100, cmap = "RdBu")

#### How variable is the climate?
std_4clim = list()
std_100clim = list()
for i in range(climate_4kya[1:].shape[1]):
    local_std_4 = np.std(climate_4kya[1:,i])
    std_4clim.append(local_std_4)
    local_std_100 = np.std(fClimate_100kya[1:,i])
    std_100clim.append(local_std_100)
    
std_4clim = [i for i in std_4clim]
std_100clim = [i for i in std_100clim]


clim_in_initial = range(0, 165)
clim_index_4 = np.zeros(coord_index.shape)
clim_index_100 = np.zeros(coord_index.shape)
for i in range(len(std_4clim)):
    for j in range(len(std_4clim)):
        if clim_in_initial[i] == coord_index[j]:
           clim_index_4[j] = std_4clim[i]
           clim_index_100[j] = std_100clim[i]
        else:
            continue
        
clim_index_4 = [i for i in clim_index_4]
clim_index_100 = [i for i in clim_index_100]

clim_4_std_in = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "STD": clim_index_4})
clim_4_std_in = clim_4_std_in.pivot_table(index='Lat', columns='Lon', values = "STD")

clim_100_std_in = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "STD": clim_index_100})
clim_100_std_in = clim_100_std_in.pivot_table(index='Lat', columns='Lon', values = "STD")


hmap1 = plt.imshow(clim_4_std_in, cmap='viridis', interpolation='none')
plt.colorbar(hmap1, ticks=[0, 0.5, 1, 1.5, 2, 2.5])

hmap2 = plt.imshow(clim_100_std_in, cmap='viridis', interpolation='none')
plt.colorbar(hmap2, ticks=[0, 0.5, 1, 1.5, 2, 2.5])

##############################################
#      Plot of species traits at single site #
##############################################
### Climate 4
clime_s = 0
clime_e =  500
niche_breadth = 0.5
for s in range(0,500):
    climate_index = s
    nt = niche_opt[s]
    t = results[s]
    
    # Select a single patch
    collists = [np.nonzero(i)[0] for i in t]
    patches_with_species = [len(i) for i in collists]
    
    # Going to use patch 9
    n_index = np.argwhere(t[29,:])
    patch_pop = t[29, n_index]
    pop_nitch = nt[29, n_index]
    patch_clim = climate_4kya[:,29][climate_index]
    
    # Patch optimal growth
    env_start = patch_clim - 7.5
    env_end = patch_clim + 7.5
    env = np.arange(env_start,env_end,0.1)
    # Get smooth curve
    yx = [100 * np.exp(-((patch_clim - e)/(2*niche_breadth))**2) for e in env]
        
    
    # Plot of Niche optimum
    plt.ioff()
    fig, axs = plt.subplots(2, 1, figsize=(12, 10))
    for i in range(len(n_index)):
        # Get population nitche
        pop_niche_opt = pop_nitch[i]
        patch_pop_n = patch_pop[i]
        # Plot 
        axs[0].scatter(pop_niche_opt, patch_pop_n, c = "black")
        axs[0].plot(env, yx)
        axs[0].set_ylim([0,100])
        axs[0].set_xlim([15,45])
        axs[0].text(15, 100,'Species Richness: %d' % len(n_index), ha='left', va='bottom', color = 'black',  size = 15)
        
    # Plot of current time step iin climate as well
    axs[1].plot(range(clime_s, clime_e), climate_4kya[clime_s:clime_e,29], c = "black")
    axs[1].axvline(climate_index, c = "red")
    axs[0].axvline(patch_clim, c = "red")
    fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/kya_4_gif/%d_step.png' % s)


# Create Gif
import os
import re
path = "/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/kya_4_gif/"
# Import files
filenames = os.listdir(path)
# Remove the DS Store
if '.DS_Store' in filenames:
  filenames.remove('.DS_Store')
# Sort
numeric_files = [int(re.match("[0-9]+", x)[0]) for x in filenames]
file_index = np.argsort(numeric_files)
file_index = [i for i in file_index]
# Organize files
filnames_sorted = [path + filenames[i] for i in file_index]

import imageio
images = []
for filename in filnames_sorted:
    try:
        images.append(imageio.imread(filename))
    except:
        continue
imageio.mimsave("/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/Narrow_Cell_9_4kya_Trait_Dynamics.gif", images)

##################
# 100 kya climate
##################

### Climate 100kya
clime_s = 0
clime_e =  500
niche_breadth = 0.5
for s in range(0,500):
    climate_index = s
    nt = niche_opt_100[s]
    t = results_100[s]
    
    # Select a single patch
    collists = [np.nonzero(i)[0] for i in t]
    patches_with_species = [len(i) for i in collists]
    
    # Going to use patch 9
    n_index = np.argwhere(t[9,:])
    patch_pop = t[29, n_index]
    pop_nitch = nt[29, n_index]
    patch_clim = fClimate_100kya[:,29][climate_index]
    
    # Patch optimal growth
    env_start = patch_clim - 7.5
    env_end = patch_clim + 7.5
    env = np.arange(env_start,env_end,0.1)
    # Get smooth curve
    yx = [100 * np.exp(-((patch_clim - e)/(2*niche_breadth))**2) for e in env]
        
    # Plot of Niche optimum
    plt.ioff()
    fig, axs = plt.subplots(2, 1, figsize=(12, 10))
    for i in range(len(n_index)):
        # Get population nitche
        pop_niche_opt = pop_nitch[i]
        patch_pop_n = patch_pop[i]
        # Plot 
        axs[0].scatter(pop_niche_opt, patch_pop_n, c = "black")
        axs[0].plot(env, yx)
        axs[0].set_ylim([0,100])
        axs[0].set_xlim([15,45])
        axs[0].text(15, 100,'Species Richness: %d' % len(n_index), ha='left', va='bottom', color = 'black',  size = 15)
        
    # Plot of current time step iin climate as well
    axs[1].plot(range(clime_s, clime_e), fClimate_100kya[clime_s:clime_e,29], c = "black")
    axs[1].axvline(climate_index, c = "red")
    axs[0].axvline(patch_clim, c = "red")
    fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/kya_100_gif/%d_step.png' % s)


# Create Gif
import os
import re
path = "/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/kya_100_gif/"
# Import files
filenames = os.listdir(path)
# Remove the DS Store
if '.DS_Store' in filenames:
  filenames.remove('.DS_Store')
# Sort
numeric_files = [int(re.match("[0-9]+", x)[0]) for x in filenames]
file_index = np.argsort(numeric_files)
file_index = [i for i in file_index]
# Organize files
filnames_sorted = [path + filenames[i] for i in file_index]

import imageio
images = []
for filename in filnames_sorted:
    try:
        images.append(imageio.imread(filename))
    except:
        continue
imageio.mimsave("/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/Narrow_Cell_9_100kya_Trait_Dynamics.gif", images)


###########
# GIF MAP #
###########
##################

### Climate 100kya
clime_s = 0
clime_e =  500
niche_breadth = 0.5
for s in range(0,500):
    t = results_100[s]
    map_in = return_richness_pivot_table(t, Clim_array, coord_index)
    ### 100 kya plot of spatial richness
    plt.ioff()
    htmap = plt.imshow(map_in.values, cmap='viridis', interpolation='nearest', vmin=0)
    #plt.colorbar(htmap,ticks=[0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    plt.axis('off')
    plt.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/map_gif_test/%d_step.png' % s)


# Create Gif
import os
import re
path = "/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/map_gif_test/"
# Import files
filenames = os.listdir(path)
# Remove the DS Store
if '.DS_Store' in filenames:
  filenames.remove('.DS_Store')
# Sort
numeric_files = [int(re.match("[0-9]+", x)[0]) for x in filenames]
file_index = np.argsort(numeric_files)
file_index = [i for i in file_index]
# Organize files
filnames_sorted = [path + filenames[i] for i in file_index]

import imageio
images = []
for filename in filnames_sorted:
    try:
        images.append(imageio.imread(filename))
    except:
        continue
imageio.mimsave("/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/SA_100kya_Map.gif", images)




