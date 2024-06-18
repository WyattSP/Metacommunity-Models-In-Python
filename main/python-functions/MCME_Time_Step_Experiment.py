#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 18:53:37 2024

@author: wyattpetryshen
"""

# Set working directory
import os
os.chdir('/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/python-functions')

# Import functions
from MCME_Speciation import MCME
import MCME_Functions
import MC_Properies_Functions
import numpy as np
from MCME_Plot_Functions import mcme_plot, species_richness_plot, plot_origination_extinction, plot_richness_map
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
import Climate_Import
import pandas as pd

###################################
#          Import climate         #
###################################

Tmax_4kyr = Climate_Import.import_climate_data("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Proc_Bricks/Tmax.mon.4kyr.tif", "Geotif", "rasterio")

#########################
#   Global Parameters   #
#########################

S = 5  # set sarting species number
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
# Intervals: 4 kya, , 8kya, 40 kya, 100 kya
# Native resolution climate for desired coordinates
climate_4kya = Climate_Import.return_climate_array(coords, 400, 901, Tmax_4kyr) - 273.5
# If you need to check output climate lengts use below function
Ct = Climate_Import.linear_climate_interpolation(climate_4kya[:,0], 25, True)
 
# Resample at new intervals
climate_16kya = Climate_Import.interpolate_climate_array(climate_4kya, 4) # skip every second points
climate_40kya = Climate_Import.interpolate_climate_array(climate_4kya, 10) # skip every fifteen points
climate_100kya = Climate_Import.interpolate_climate_array(climate_4kya, 25) # skip every fifty points
# Put climates into a giant array
climate_input_list = list([climate_4kya, climate_16kya, climate_40kya, climate_100kya])
climate_input_names = ["4kya_step", "12kya_step", "40kya_step", "100kya_step"]

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

# Put climates into a giant array
climate_input_list = list([climate_4kya, fClimate_16kya, fClimate_40kya, fClimate_100kya])
climate_input_names = ["4kya_step", "fClimate_16kya", "fClimate_40kya", "fClimate_100kya"]

###################################
#    Set up Model Parameters      #
###################################
# Max Simulation time
sim_time = climate_4kya.shape[0] - 1

# Simulation times
time_in = [200, 10, 200, sim_time] # seed time, seed inteval, burn-in time, simulation time

# Set intial species niche optimums; this will stay static between all simulations
niche_optimum_ini = MCME_Functions.initial_species_niche_optimum(S, M, climate_4kya[0,:])

###################################
#     Iteratable Parameters       #
###################################

# Dispersal ability
#dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16, 1e-5, 1)
dispersal_rate_ini = [0.001, 0.01, 0.1]

# Niche Breadth
niche_breadth = [1, 10]

# Specify the speciation rate
speciation_rate = [0.01]

# Species interaction end-members
#end_member = ["equal", "stabalizing", "mixed", "neutral"]
end_member = ["stabalizing"]

# Equal species interactions
alpha_matrix_equal = MCME_Functions.initalize_aij("equal", 0.0, 1.0, S, 0.3)
# Stabalizing species interactions
alpha_matrix_stabalizing = MCME_Functions.initalize_aij("stabalizing", 0.0, 1.0, S, 0.3)
# Mixed species interactions
alpha_matrix_mixed = MCME_Functions.initalize_aij("mixed", 0.0, 1.5, S, 0.3)
# Neutral species interactions
#alpha_matrix_neutral = MCME_Functions.initalize_aij("neutral", 0.0, 1.0, S, 0.3)

alpha_matrix_list = list([alpha_matrix_equal, alpha_matrix_stabalizing, alpha_matrix_mixed])

###################################
#        Simulation Loops         #
###################################

# Save Path 
save_path = "/Volumes/Wyatt_SSD_A/MCME_Results/"

for reps in range(10):
    simulation_results = list()
    simulation_phylogeny = list()
    s_names = list()
    run_num = 1
    for t_i, t in enumerate(climate_input_names):
        for n in niche_breadth:
            for sp in speciation_rate:
                for end_i, end_m in enumerate(end_member):
                    for disp in dispersal_rate_ini:
                        num_sims = len(climate_input_names) * len(niche_breadth) * len(speciation_rate) * len(end_member) * len(dispersal_rate_ini)
                        # Set save name
                        name = "Rep:" + str(reps) + "_Clim:" + str(t) + "_n:" + str(n) + "_sp:" + str(sp) + "_" + str(end_m) + "_disp:" + str(round(disp,3))
                        s_names.append(name)
                        
                        # Run simulation
                        results, niche_opt, alpha_matrix, phylogeny, divergence_time = MCME(time_in,
                                                                                           S,
                                                                                           max_growth_r,
                                                                                           sp,
                                                                                           alpha_matrix_stabalizing,
                                                                                           distance_matrix,
                                                                                           climate_input_list[t_i],
                                                                                           niche_optimum_ini,
                                                                                           disp,
                                                                                           n,
                                                                                           end_m)
                        
                        # Place into save format
                        # Set full save path
                        full_save_path = save_path + name + ".csv"
                        # Save to disk as csv
                        #np.savetxt(full_save_path, results, delimiter=",")
                        simulation_results.append(results)
                        simulation_phylogeny.append(phylogeny)
                        
                        print("%d" % run_num, " of %d" % num_sims)
                        run_num += 1
                        

###########################
#     Calculate Gammas    #
###########################

# Calculate gamma and alpha

stabalizing_gamma = list()
sstabalizing_alpha = list()

end_slice = len(simulation_results[0]) - 1

for i, output in enumerate(simulation_results):
    t_gamma = MC_Properies_Functions.gamma_diversity(output)[-1]
    t_alpha = MC_Properies_Functions.alpha_richness(output)[-1]
    stabalizing_gamma.append(t_gamma)
    sstabalizing_alpha.append(t_alpha)
    

color_scatter = ["blue","blue","blue","green","green","green","red","red","red","orange","orange","orange"]
plt.scatter(range(12), stabalizing_gamma, label = '_nolegend_' , color = color_scatter, marker = "+")
plt.scatter(range(12), sstabalizing_alpha, label = '_nolegend_', color = color_scatter, marker = "x")
plt.plot(range(12), stabalizing_gamma, c = "blue", label = "gamma")
plt.plot(range(12), sstabalizing_alpha, c = "green", label = "alpha")
plt.xlabel("Experiment Combinations")
plt.ylabel("Richness")
plt.xticks(range(12), s_names, rotation='vertical')
plt.title("")
plt.legend(loc="upper left")

plot_richness_map(results[-1], Clim_array, coord_index, False)
plot_richness_map(results[-1], Clim_array, coord_index, False)

plot_richness_map(simulation_results[0][-1], Clim_array, coord_index, False)
plot_richness_map(simulation_results[3][-1], Clim_array, coord_index, False)
plot_richness_map(simulation_results[6][-1], Clim_array, coord_index, False)
plot_richness_map(simulation_results[9][-1], Clim_array, coord_index, False)

orgination, extinction = plot_origination_extinction(simulation_results[0], True)
orgination, extinction = plot_origination_extinction(simulation_results[3], True)
orgination, extinction = plot_origination_extinction(simulation_results[6], True)
orgination, extinction = plot_origination_extinction(simulation_results[9], True)

###########################
#     Individual Plots    #
###########################

gamma_all = list()
alpha_all = list()

end_slice = len(simulation_results[0]) - 1

for i, output in enumerate(simulation_results):
    t_gamma = MC_Properies_Functions.gamma_diversity(output)
    t_alpha = MC_Properies_Functions.alpha_richness(output)
    gamma_all.append(t_gamma)
    alpha_all.append(t_alpha)
    
    
x = np.arange(0,len(alpha_all[0]))
plt.plot(x, gamma_all[0], c = "blue", linestyle = "solid", label = "4 kya gamma")
plt.plot(x, alpha_all[0], c = "blue", linestyle = "dashed", label = "4 kya alpha")
plt.plot(x, gamma_all[3], c = 'green',linestyle = "solid", label = "8 kya gamma")
plt.plot(x, alpha_all[3], c = "green", linestyle = "dashed", label = "8 kya alpha")
plt.plot(x, gamma_all[6], c = 'red', linestyle = "solid", label = "40 kya gamma")
plt.plot(x, alpha_all[6], c = "red", linestyle = "dashed", label = "40 kya alpha")
plt.plot(x, gamma_all[9], c = 'orange', linestyle = "solid", label = "100 kya gamma")
plt.plot(x, alpha_all[9], c = "orange", linestyle = "dashed", label = "100 kya alpha")
plt.title("Low Dispersal")
plt.legend(loc="upper left")

x = np.arange(0,len(alpha_all[0]))
plt.plot(x, gamma_all[2], c = "blue", linestyle = "solid", label = "4 kya gamma")
plt.plot(x, alpha_all[2], c = "blue", linestyle = "dashed", label = "4 kya alpha")
plt.plot(x, gamma_all[5], c = 'green', linestyle = "solid", label = "8 kya gamma")
plt.plot(x, alpha_all[5], c = "green", linestyle = "dashed", label = "8 kya alpha")
plt.plot(x, gamma_all[8], c = 'red', linestyle = "solid", label = "40 kya gamma")
plt.plot(x, alpha_all[8], c = "red", linestyle = "dashed", label = "40 kya alpha")
plt.plot(x, gamma_all[11], c = 'orange', linestyle = "solid", label = "100 kya gamma")
plt.plot(x, alpha_all[11], c = "orange", linestyle = "dashed", label = "100 kya alpha")
plt.title("High Dispersal")
plt.legend(loc="upper left")

###########################
#  Climate Plot Example   #
###########################
from scipy import signal
freq4, spec4 = signal.periodogram(climate_4kya[:,1], detrend="constant")

fig, axs = plt.subplots(3, 1, figsize=(12, 8))
axs[0].plot(range(len(climate_4kya[:,1])),climate_4kya[:,1])
axs[1].plot(freq4,spec4)
axs[2].plot(np.log(freq4), np.log(spec4))

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(12, 16))
fig.tight_layout(pad=2)
axs[0].plot(range(len(climate_4kya[:,1])),climate_4kya[:,1], color = 'black')
axs[0].set_xlabel("Thousand Years")
axs[0].set_ylabel("T ($^\circ$C)")
axs[0].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[0].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[0].tick_params(axis='x', colors='darkgrey')
axs[0].tick_params(axis='y', colors='darkgrey')
axs[0].spines['bottom'].set_color('darkgrey')
axs[0].spines['top'].set_color('darkgrey') 
axs[0].spines['right'].set_color('darkgrey')
axs[0].spines['left'].set_color('darkgrey')
axs[0].yaxis.label.set_color('darkgrey')
axs[0].xaxis.label.set_color('darkgrey')

axs[1].plot(freq4,spec4, color = 'black')
axs[1].set_xlabel("Frequency")
axs[1].set_ylabel("Power")
axs[1].tick_params(axis='x', colors='darkgrey')
axs[1].tick_params(axis='y', colors='darkgrey')
axs[1].spines['bottom'].set_color('darkgrey')
axs[1].spines['top'].set_color('darkgrey') 
axs[1].spines['right'].set_color('darkgrey')
axs[1].spines['left'].set_color('darkgrey')
axs[1].yaxis.label.set_color('darkgrey')
axs[1].xaxis.label.set_color('darkgrey')

axs[2].plot(np.log(freq4), np.log(spec4), color = 'black')
axs[2].set_xlabel("Log Frequency")
axs[2].set_ylabel("Log Power")
axs[2].tick_params(axis='x', colors='darkgrey')
axs[2].tick_params(axis='y', colors='darkgrey')
axs[2].spines['bottom'].set_color('darkgrey')
axs[2].spines['top'].set_color('darkgrey') 
axs[2].spines['right'].set_color('darkgrey')
axs[2].spines['left'].set_color('darkgrey')
axs[2].yaxis.label.set_color('darkgrey')
axs[2].xaxis.label.set_color('darkgrey')

fig.patch.set_facecolor('mintcream')
axs[0].set_facecolor('mintcream')
axs[1].set_facecolor('mintcream')
axs[2].set_facecolor('mintcream')


fig, axs = plt.subplots(figsize=(12, 4))
fig.tight_layout(pad=2)
fig.patch.set_facecolor('mintcream')
axs.plot(range(len(climate_4kya[:,1])),climate_4kya[:,1], color = 'black')
axs.set_xlabel("Thousand Years")
axs.set_ylabel("T ($^\circ$C)")
axs.xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs.xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs.tick_params(axis='x', colors='darkgrey')
axs.tick_params(axis='y', colors='darkgrey')
axs.spines['bottom'].set_color('darkgrey')
axs.spines['top'].set_color('darkgrey') 
axs.spines['right'].set_color('darkgrey')
axs.spines['left'].set_color('darkgrey')
axs.yaxis.label.set_color('darkgrey')
axs.xaxis.label.set_color('darkgrey')
axs.set_facecolor('mintcream')

fig, axs = plt.subplots(figsize=(12, 4))
fig.tight_layout(pad=2)
fig.patch.set_facecolor('mintcream')
axs.plot(freq4,spec4, color = 'black')
axs.set_xlabel("Frequency")
axs.set_ylabel("Power")
axs.tick_params(axis='x', colors='darkgrey')
axs.tick_params(axis='y', colors='darkgrey')
axs.spines['bottom'].set_color('darkgrey')
axs.spines['top'].set_color('darkgrey') 
axs.spines['right'].set_color('darkgrey')
axs.spines['left'].set_color('darkgrey')
axs.yaxis.label.set_color('darkgrey')
axs.xaxis.label.set_color('darkgrey')
axs.set_facecolor('mintcream')

fig, axs = plt.subplots(figsize=(12, 4))
fig.tight_layout(pad=2)
fig.patch.set_facecolor('mintcream')
axs.plot(np.log(freq4), np.log(spec4), color = 'black')
axs.set_xlabel("Log Frequency")
axs.set_ylabel("Log Power")
axs.tick_params(axis='x', colors='darkgrey')
axs.tick_params(axis='y', colors='darkgrey')
axs.spines['bottom'].set_color('darkgrey')
axs.spines['top'].set_color('darkgrey') 
axs.spines['right'].set_color('darkgrey')
axs.spines['left'].set_color('darkgrey')
axs.yaxis.label.set_color('darkgrey')
axs.xaxis.label.set_color('darkgrey')
axs.set_facecolor('mintcream')


###########################
#  Line Plot of Climate   #
###########################

fig, axs = plt.subplots(7, 1, figsize=(13, 10))
#fig.tight_layout(pad=0.5)
fig.patch.set_facecolor('mintcream')
axs[0].plot(range(len(climate_4kya[:,1])), climate_4kya[:,1], color='black', linestyle='-')
axs[0].text(0, 25,'4 kya Native Resolution', ha='left', va='bottom', color = 'darkgrey')
axs[0].set_xticks([])
axs[0].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[0].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[0].axhline(y = 40,lw=0.5, color='black', alpha=0.5)

axs[1].plot(range(len(climate_16kya[:,1])), climate_16kya[:,1], color='red', linestyle='-')
axs[1].set_xticks([])
axs[1].text(0, 25,'16 kya Smoothed', ha='left', va='bottom', color = 'darkgrey')
axs[1].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[1].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[1].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
#axs[1].set_facecolor('xkcd:light moss green')

axs[2].plot(range(len(climate_40kya[:,1])), climate_40kya[:,1], color='blue', linestyle='-')
axs[2].set_xticks([])
axs[2].text(0, 25,'40 kya Smoothed', ha='left', va='bottom', color = 'darkgrey')
axs[2].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[2].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[2].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
#axs[2].set_facecolor('xkcd:light moss green')

axs[3].plot(range(len(climate_100kya[:,1])), climate_100kya[:,1], color='green', linestyle='-')
axs[3].set_xticks([])
axs[3].text(0, 25,'100 kya Smoothed', ha='left', va='bottom', color = 'darkgrey')
axs[3].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[3].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[3].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
axs[3].set_ylabel("T ($^\circ$C)")
#axs[3].set_facecolor('xkcd:light moss green')

axs[4].plot(range(len(fClimate_16kya[:,1])), fClimate_16kya[:,1], color='red', linestyle='-')
axs[4].set_xticks([])
axs[4].text(0, 25,'16 kya Fourier Transformed', ha='left', va='bottom', color = 'darkgrey')
axs[4].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[4].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[4].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
#axs[4].set_facecolor('xkcd:cloudy blue')

axs[5].plot(range(len(fClimate_40kya[:,1])), fClimate_40kya[:,1], color='blue', linestyle='-')
axs[5].set_xticks([])
axs[5].text(0, 25,'40 kya Fourier Transformed' , ha='left', va='bottom', color = 'darkgrey')
axs[5].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[5].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[5].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
#axs[5].set_facecolor('xkcd:cloudy blue')

axs[6].plot(range(len(fClimate_100kya[:,1])), fClimate_100kya[:,1], color='green', linestyle='-')
axs[6].text(0, 25,'100 kya Fourier Transformed', ha='left', va='bottom', color = 'darkgrey')
axs[6].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[6].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[6].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
axs[6].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[6].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[6].set_xlabel("Thousand Years")
#axs[6].set_facecolor('xkcd:cloudy blue')

for i in range(7):
    axs[i].tick_params(axis='x', colors='darkgrey')
    axs[i].tick_params(axis='y', colors='darkgrey')
    axs[i].spines['bottom'].set_color('darkgrey')
    axs[i].spines['top'].set_color('darkgrey') 
    axs[i].spines['right'].set_color('darkgrey')
    axs[i].spines['left'].set_color('darkgrey')
    axs[i].yaxis.label.set_color('darkgrey')
    axs[i].xaxis.label.set_color('darkgrey')
    axs[i].set_facecolor('mintcream')
    axs[i].set_ylim([25,35])

#########################################
# Figure of Climate Fourier Transformed #
#########################################
from scipy import signal
# Original
freq4, spec4 = signal.periodogram(climate_4kya[:,1], detrend="constant")
# Smoothed
freq12, spec12 = signal.periodogram(climate_16kya[:,1], detrend="constant")
freq40, spec40 = signal.periodogram(climate_40kya[:,1], detrend="constant")
freq100, spec100 = signal.periodogram(climate_100kya[:,1], detrend="constant")
# Fourier transformed
freq_f12, spec_f12 = signal.periodogram(fClimate_16kya[:,1], detrend="constant")
freq_f40, spec_f40 = signal.periodogram(fClimate_40kya[:,1], detrend="constant")
freq_f100, spec_f100 = signal.periodogram(fClimate_100kya[:,1], detrend="constant")

fig, axs = plt.subplots(7, 1, figsize=(12, 8))
fig.patch.set_facecolor('mintcream')
axs[0].plot(np.log(freq4), np.log(spec4), color='black', linestyle='-')
axs[0].text(0-6, -8,'4 kyr Native Resolution', ha='left', va='bottom', color = 'darkgrey')
axs[0].set_xticks([])
axs[0].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[0].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[0].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[0].axhline(y = 0,lw=0.5, color='black', alpha=0.5)


axs[1].plot(np.log(freq12), np.log(spec12), color='red', linestyle='-')
axs[1].set_xticks([])
axs[1].text(-6, -8,'16 kyr Smoothed', ha='left', va='bottom', color = 'darkgrey')
axs[1].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[1].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[1].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[1].axhline(y = 0,lw=0.5, color='black', alpha=0.5)

axs[2].plot(np.log(freq40), np.log(spec40), color='blue', linestyle='-')
axs[2].set_xticks([])
axs[2].text(-6, -8,'40 kyr Smoothed', ha='left', va='bottom', color = 'darkgrey')
axs[2].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[2].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[2].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[2].axhline(y = 0,lw=0.5, color='black', alpha=0.5)

axs[3].plot(np.log(freq100), np.log(spec100), color='green', linestyle='-')
axs[3].set_xticks([])
axs[3].text(-6, -8,'100 kyr Smoothed', ha='left', va='bottom', color = 'darkgrey')
axs[3].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[3].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[3].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[3].axhline(y = 0,lw=0.5, color='black', alpha=0.5)
axs[3].set_ylabel("Power")

axs[4].plot(np.log(freq_f12), np.log(spec_f12), color='red', linestyle='-')
axs[4].set_xticks([])
axs[4].text(-6, -8, '16 kyr Fourier Transformed', ha='left', va='bottom', color = 'darkgrey')
axs[4].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[4].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[4].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[4].axhline(y = 0,lw=0.5, color='black', alpha=0.5)

axs[5].plot(np.log(freq_f40), np.log(spec_f40), color='blue', linestyle='-')
axs[5].set_xticks([])
axs[5].text(-6, -8,'40 kyr Fourier Transformed' , ha='left', va='bottom', color = 'darkgrey')
axs[5].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[5].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[5].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[5].axhline(y = 0,lw=0.5, color='black', alpha=0.5)

axs[6].plot(np.log(freq_f100), np.log(spec_f100), color='green', linestyle='-')
axs[6].text(-6, -8,'100 kyr Fourier Transformed', ha='left', va='bottom', color = 'darkgrey')
axs[6].axvspan(np.log(Climate_Import.year_to_frequency(25,4)), np.log(Climate_Import.year_to_frequency(18,4)), alpha=0.5, color='grey')
axs[6].axvspan(np.log(Climate_Import.year_to_frequency(56,4)), np.log(Climate_Import.year_to_frequency(53,4)), alpha=0.5, color='grey')
axs[6].axvspan(np.log(Climate_Import.year_to_frequency(126,4)), np.log(Climate_Import.year_to_frequency(122,4)), alpha=0.5, color='grey')
axs[6].axhline(y = 0,lw=0.5, color='black', alpha=0.5)
axs[6].set_xlabel("Frequency")

for i in range(7):
    axs[i].set_ylim([-15,10])
    axs[i].tick_params(axis='x', colors='darkgrey')
    axs[i].tick_params(axis='y', colors='darkgrey')
    axs[i].spines['bottom'].set_color('darkgrey')
    axs[i].spines['top'].set_color('darkgrey') 
    axs[i].spines['right'].set_color('darkgrey')
    axs[i].spines['left'].set_color('darkgrey')
    axs[i].yaxis.label.set_color('darkgrey')
    axs[i].xaxis.label.set_color('darkgrey')
    axs[i].set_facecolor('mintcream')
    axs[i].set_ylim([-8,8])
    axs[i].set_xlim([-6.5,0])
    
    

# Estimate slopes between a range
# Get slopes of fit
low = Climate_Import.year_to_frequency(126,4)
high = Climate_Import.year_to_frequency(18,4)

a4 = np.where(freq4 >= low)
b4 = np.where(freq4 <= high)
index_int_4 = np.intersect1d(a4,b4)

freq_filt_4 = freq4[index_int_4]
spec_filt_4 = freq4[index_int_4]

freq_filt_16 = freq12[index_int_4]
spec_filt_16 = spec12[index_int_4]

freq_filt_40 = freq40[index_int_4]
spec_filt_40 = spec40[index_int_4]

freq_filt_100 = freq100[index_int_4]
spec_filt_100 = spec100[index_int_4]

freq_filt_f16 = freq_f12[index_int_4]
spec_filt_f16 = spec_f12[index_int_4]


trend_4 = np.polyfit(np.log(freq_filt_4), np.log(spec_filt_4), deg = 1)
trend_line_4 = np.poly1d(trend_4)

trend_16 = np.polyfit(np.log(freq_filt_16), np.log(spec_filt_16), deg = 1)
trend_line_16 = np.poly1d(trend_16)

trend_40 = np.polyfit(np.log(freq_filt_40), np.log(spec_filt_40), deg = 1)
trend_line_40 = np.poly1d(trend_40)

trend_100 = np.polyfit(np.log(freq_filt_100), np.log(spec_filt_100), deg = 1)
trend_line_100 = np.poly1d(trend_100)

trend_f16 = np.polyfit(np.log(freq_filt_f16), np.log(spec_filt_f16), deg = 1)
trend_line_f16 = np.poly1d(trend_f16)

trend_4[0]
trend_16[0] # -3.6819645143338438
trend_40[0]
trend_100[0]
trend_f16[0] # 0.02206042894767503

fig, axs = plt.subplots(2, 1,figsize=(12, 8))
axs[0].plot(np.log(freq40), np.log(spec40), color='blue', linestyle='-')
axs[0].set_ylim([-20,10])
axs[0].set_xlim([-6.5,0])
axs[0].scatter(np.log(freq_filt_40),np.log(spec_filt_40))
axs[0].plot(np.log(freq_filt_40), trend_line_40(np.log(freq_filt_40)), c = "black")

axs[1].plot(np.log(freq12), np.log(spec12), color='black', linestyle='-')
axs[1].set_ylim([-20,10])
axs[1].set_xlim([-6.5,0])
axs[1].scatter(np.log(freq_filt_16),np.log(spec_filt_16))
axs[1].plot(np.log(freq_filt_16), trend_line_16(np.log(freq_filt_16)), c = "black")


axs.set_ylim([-20,10])
axs.set_xlim([-6.5,0])
axs.scatter(np.log(freq_filt_40),np.log(spec_filt_40))
axs.plot(np.log(freq_filt_40), trend_line_40(np.log(freq_filt_40)), c = "black")


fig, axs = plt.subplots(figsize=(12, 8))
axs.plot(np.log(freq4), np.log(spec4), color='black', linestyle=(0, (5, 5)), alpha=0.5, label = "4 kya")
axs.plot(np.log(freq12), np.log(spec12), color='blue', linestyle=(0, (5, 5)), alpha=0.5, label = "16 kya")
axs.plot(np.log(freq40), np.log(spec40), color='red', linestyle=(0, (5, 5)), alpha=0.5, label = "40 kya")
axs.plot(np.log(freq100), np.log(spec100), color='green', linestyle=(0, (5, 5)), alpha=0.5, label = "100 kya")

axs.plot(np.log(freq_filt_4), trend_line_4(np.log(freq_filt_4)), c = "black")
axs.plot(np.log(freq_filt_16), trend_line_16(np.log(freq_filt_16)), c = "blue")
axs.plot(np.log(freq_filt_40), trend_line_40(np.log(freq_filt_40)), c = "red")
axs.plot(np.log(freq_filt_100), trend_line_100(np.log(freq_filt_100)), c = "green")



#### Fourier Only

###########################
#  Line Plot of Climate   #
###########################

fig, axs = plt.subplots(4, 1, figsize=(13, 10))
#fig.tight_layout(pad=0.5)
fig.patch.set_facecolor('mintcream')
axs[0].plot(range(len(climate_4kya[1:,1])), climate_4kya[1:,1], color='black', linestyle='-', zorder = 1)
axs[0].scatter(range(len(climate_4kya[1:,1])), climate_4kya[1:,1], color = 'black', s = 0.5, zorder = 2)
axs[0].text(0, 25,'4 kya Native Resolution', ha='left', va='bottom', color = 'black',  size = 15)
axs[0].set_xticks([])
axs[0].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[0].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[0].axhline(y = 40,lw=0.5, color='black', alpha=0.5)

axs[1].plot(range(len(fClimate_16kya[:,1])), fClimate_16kya[:,1], color='green', linestyle='-',zorder = 1)
axs[1].scatter(range(len(fClimate_16kya[:,1])), fClimate_16kya[:,1], color = 'black', s = 0.5, zorder = 2)
axs[1].set_xticks([])
axs[1].text(0, 25,'16 kya Fourier Transformed', ha='left', va='bottom', color = 'black',  size = 15)
axs[1].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[1].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[1].axhline(y = 40,lw=0.5, color='black', alpha=0.5)

axs[2].plot(range(len(fClimate_40kya[:,1])), fClimate_40kya[:,1], color='blue', linestyle='-', zorder = 1)
axs[2].scatter(range(len(fClimate_40kya[:,1])), fClimate_40kya[:,1], color = 'black', s = .5, zorder = 2)
axs[2].set_xticks([])
axs[2].text(0, 25,'40 kya Fourier Transformed' , ha='left', va='bottom', color = 'black',  size = 15)
axs[2].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[2].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[2].axhline(y = 40,lw=0.5, color='black', alpha=0.5)

axs[3].plot(range(len(fClimate_100kya[:,1])), fClimate_100kya[:,1], color='red', linestyle='-', zorder = 1)
axs[3].scatter(range(len(fClimate_100kya[:,1])), fClimate_100kya[:,1], color = 'black', s = 0.5, zorder = 2)
axs[3].text(0, 25,'100 kya Fourier Transformed', ha='left', va='bottom', color = 'black',  size = 15)
axs[3].axhline(y = 20,lw=0.5, color='black', alpha=0.5)
axs[3].axhline(y = 30,lw=0.5, color='black', alpha=0.5)
axs[3].axhline(y = 40,lw=0.5, color='black', alpha=0.5)
axs[3].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[3].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[3].set_xlabel("Thousand Years",  size = 15)

for i in range(4):
    axs[i].tick_params(axis='x', colors='black')
    axs[i].tick_params(axis='y', colors='black')
    axs[i].spines['bottom'].set_color('black')
    axs[i].spines['top'].set_color('black') 
    axs[i].spines['right'].set_color('black')
    axs[i].spines['left'].set_color('black')
    axs[i].yaxis.label.set_color('black')
    axs[i].xaxis.label.set_color('black')
    axs[i].set_facecolor('mintcream')
    axs[i].set_ylim([25,35])
    axs[i].text(475, 33,'n = 500', ha='left', va='bottom', color = 'black',  size = 15)
    axs[i].set_ylabel("T ($^\circ$C)", size = 15)


fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/fourier_steps_v2.png', format='png', dpi=300, transparent=True)



#### Line plot of smoothed climate


##############################################
#  Line Plot Example of Different Time Steps #
##############################################
fig, axs = plt.subplots(3, 1, figsize=(13, 10))
#fig.tight_layout(pad=0.5)
axs[0].plot(range(len(climate_4kya[:,43])), climate_4kya[:,43], color='black', linestyle='-')
axs[0].text(0, 26,'Climate Forcings ≥ 4000-Years', ha='left', va='bottom', color = 'black', size = 18)
axs[0].set_xticks([])
#axs[0].axhline(y = 29,lw=0.5, color='black', alpha=0.5)

axs[1].plot(range(len(climate_40kya[:,43])), climate_40kya[:,43], color='blue', linestyle='-')
axs[1].set_xticks([])
axs[1].text(0, 26,'Climate Forcings ≥ 40,000-Years', ha='left', va='bottom', color = 'black', size = 18)
#axs[1].axhline(y = 29,lw=0.5, color='black', alpha=0.5)

axs[2].plot(range(len(climate_100kya[:,43])), climate_100kya[:,43], color='green', linestyle='-')
axs[2].set_xticks([])
axs[2].text(0, 26,'Climate Forcings ≥ 100,000-Years', ha='left', va='bottom', color = 'black', size = 18)
#axs[2].axhline(y = 29,lw=0.5, color='black', alpha=0.5)
axs[2].set_ylabel("T ($^\circ$C)")

for i in range(3):
    axs[i].tick_params(axis='x', colors='black')
    #axs[i].tick_params(axis='y', colors='black')
    axs[i].spines['bottom'].set_color('black')
    axs[i].spines['top'].set_color('black') 
    axs[i].spines['right'].set_color('black')
    axs[i].spines['left'].set_color('black')
    axs[i].yaxis.label.set_color('black')
    axs[i].xaxis.label.set_color('black')
    axs[i].set_ylim([26,31])
    axs[i].set_ylabel("T ($^\circ$C)", size = 15)
    axs[i].grid(False)

fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/skip_time_ex.png', format='png', dpi=300, transparent=True)




