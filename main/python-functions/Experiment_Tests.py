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
import butterworth_filter
import matplotlib.pyplot as plt
import sys
import Climate_Import
import pandas as pd
import random
import math
import matplotlib

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

# Interpolated Climate
# Resample at new intervals
climate_16kya = Climate_Import.interpolate_climate_array(climate_4kya, 4) # skip every second points
climate_40kya = Climate_Import.interpolate_climate_array(climate_4kya, 10) # skip every fifteen points
climate_100kya = Climate_Import.interpolate_climate_array(climate_4kya, 25) # skip every fifty points

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

##########################################
#      Butterworth Low-Pass Filter       #
##########################################
sample_rate = 4 # t series sample rate
# double sampling 4 * 3 = 16 kya
lfreq_1 = round(1 / (16 / sample_rate), 3)
# 10 times sampling 4 * 10 = 40 kya
lfreq_2 = round(1 / (40 / sample_rate), 3)
# 25 times sampling 4 * 25 = 100 kya
lfreq_3 = round(1 / (100 / sample_rate), 3)

# Note that for evenly spaced data you will use analog = False in the butter function 
lowpass_freq = list([lfreq_1, lfreq_2, lfreq_3])

BFc16 = butterworth_filter.get_the_butter(climate_array = climate_4kya[1:], lpass = lfreq_1, order = 9)
BFc40 = butterworth_filter.get_the_butter(climate_array = climate_4kya[1:], lpass = lfreq_2, order = 9)
BFc100 = butterworth_filter.get_the_butter(climate_array = climate_4kya[1:],lpass = lfreq_3, order = 9)

###################################
#    Set up Model Parameters      #
###################################
# Max Simulation time
sim_time = climate_4kya.shape[0]

# Simulation times
time_in = [200, 10, 100, sim_time] # seed time, seed inteval, burn-in time, simulation time

# Set intial species niche optimums; this will stay static between all simulations
niche_optimum_ini_4 = MCME_Functions.initial_species_niche_optimum(S, M, climate_4kya[0,:])
niche_optimum_ini_100 = MCME_Functions.initial_species_niche_optimum(S, M, BFc100[0,:])

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

# 4 kya forcings wide
results, niche_opt, alpha_matrix, phylogeny, divergence_time, patch_origin = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_equal,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = climate_4kya,
                                                                   niche_optimum_ini = niche_optimum_ini_4,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 2.5,
                                                                   end_member = "equal")


# 100 kya forcings wide
results_100, niche_opt_100, alpha_matrix_100, phylogeny_100, divergence_time_100, patch_origin_100 = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = BFc100,
                                                                   niche_optimum_ini = niche_optimum_ini_100,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 2.5,
                                                                   end_member = "stabalizing")

# 4 kya forcings narrow
n_results, n_niche_opt, n_alpha_matrix, n_phylogeny, n_divergence_time, n_patch_origin = MCME(time_in =  [100, 10, 100, 500],
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


# 100 kya forcings narrow
n_results_100, n_niche_opt_100, n_alpha_matrix_100, n_phylogeny_100, n_divergence_time_100, n_patch_origin_100 = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = BFc100,
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

n_gamma = MC_Properies_Functions.gamma_diversity(results)
n_alpha = MC_Properies_Functions.alpha_richness(results)
n_gamma_100 = MC_Properies_Functions.gamma_diversity(n_results_100)
n_alpha_100 = MC_Properies_Functions.alpha_richness(n_results_100)

# Hill Numbers
# Get mean hill numbers between all patchs? 
# Let's try for the last time step
ls_16 = MC_Properies_Functions.compare_patch_by_hill(results[-1], np.arange(0,4,1), 1, True)
ls_100 = MC_Properies_Functions.compare_patch_by_hill(results_100[-1], np.arange(0,4,1), 1, True)
# Plot of q orders
MCME_Plot_Functions.plot_hill_number_maps(ls_100 - ls_16, Clim_array, coord_index, False)

delta_ls_16 = MC_Properies_Functions.compare_hill_difference(results[-1], [0,2], 1, True)
delta_ls_100 = MC_Properies_Functions.compare_hill_difference(results_100[-1], [0,2], 1, True)

sum(delta_ls_16[~np.isnan(delta_ls_16)])
sum(delta_ls_100[~np.isnan(delta_ls_100)])

# Plot Hill differences
MCME_Plot_Functions.plot_hill_difference_map(delta_ls_16, Clim_array, coord_index, False)
MCME_Plot_Functions.plot_hill_difference_map(delta_ls_100, Clim_array, coord_index, False)
#########################
#          Plot         # 
#########################

# Plot of gamma and alpha diversity
x = np.arange(0,len(n_alpha))
plt.plot(x, n_gamma, label = "gamma")
plt.plot(x, n_alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")


# Plot of species richness change
orgination, extinction = plot_origination_extinction(results, True)
plt.plot(range(len(np.cumsum(orgination[:,1]))), np.cumsum(orgination[:,1]), c = "blue")
plt.plot(range(len(np.cumsum(extinction[:,1]))), np.cumsum(extinction[:,1]), c = "red")

# Plot map of geographic points with colours illustrating community richness
plot_richness_map(results_100[-1], Clim_array, coord_index, False)

#np.savetxt("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/test_out_raster.csv", t, delimiter=",")

##########################################
#        Plots of Climate Inputs         #
##########################################
# Mean Climate Untransformed
climate_mean = np.mean(climate_4kya[1:],1)
# Mean Butterworths
Bfc16_mean = np.mean(BFc16,1)
Bfc40_mean = np.mean(BFc40,1)
Bfc100_mean = np.mean(BFc100,1)
# Mean Fourier
fC16_mean = np.mean(fClimate_16kya, 1)
fC40_mean = np.mean(fClimate_40kya, 1)
fC100_mean = np.mean(fClimate_100kya, 1)
# Mean Interpolations
Sc16_mean = np.mean(climate_16kya, 1)
Sc40_mean = np.mean(climate_40kya, 1)
Sc100_mean = np.mean(climate_100kya, 1)

##########################
#### Plot Time Domain ####
##########################
fig, axs = plt.subplots(4, 1, figsize=(12, 8))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
# Plot of Original Data
axs[0].plot(range(len(climate_mean)), climate_mean, c = "black")
axs[0].set_ylim([23, 28])
axs[0].text(0, 23.25,'Mean Untransformed Climate')
axs[0].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[0].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[0].set_ylabel('Temperature (\u00B0C)')

# Butterworth Plot
axs[1].plot(range(len(Bfc16_mean)), Bfc16_mean)
axs[1].plot(range(len(Bfc40_mean)), Bfc40_mean)
axs[1].plot(range(len(Bfc100_mean)), Bfc100_mean)
axs[1].set_ylim([23, 28])
axs[1].text(0, 23.25,'Mean Butterworth Climate')
axs[1].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[1].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[1].set_ylabel('Temperature (\u00B0C)')
# Fourier Plot
axs[2].plot(range(len(fC16_mean)), fC16_mean)
axs[2].plot(range(len(fC40_mean)), fC40_mean)
axs[2].plot(range(len(fC100_mean)), fC100_mean)
axs[2].set_ylim([23, 28])
axs[2].text(0, 23.25,'Mean Fourier Climate')
axs[2].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[2].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[2].set_ylabel('Temperature (\u00B0C)')
# Smooth Interpolated Plot
axs[3].plot(range(len(Sc16_mean)), Sc16_mean)
axs[3].plot(range(len(Sc40_mean)), Sc40_mean)
axs[3].plot(range(len(Sc100_mean)), Sc100_mean)
axs[3].set_ylim([23, 28])
axs[3].text(0, 23.25,'Mean Smooth Interpolated Climate')
axs[3].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[3].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[3].set_xlabel("Thousand Years")
axs[3].set_ylabel('Temperature (\u00B0C)')
# Legend
fig.legend(labels=['4 kya', '16 kya', '40 kya', '100 kya'], loc='lower left', 
           shadow=False, ncol=4, bbox_to_anchor=(0.075, 0.05))

fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/git_time_domain-env.png', format='png', dpi=300, transparent=False)

#########################################
# Plot of Climate Fourier Transformed #
#########################################
from scipy import signal
# Original
Or_fr, Sp_fr = signal.periodogram(climate_mean, detrend="constant")
# Butterworth
Bfc16_fr, Bfc16_sp = signal.periodogram(Bfc16_mean, detrend="constant")
Bfc40_fr, Bfc40_sp = signal.periodogram(Bfc40_mean, detrend="constant")
Bfc100_fr, Bfc100_sp = signal.periodogram(Bfc100_mean, detrend="constant")
# Fourier
freq_f12, spec_f12 = signal.periodogram(fC16_mean, detrend="constant")
freq_f40, spec_f40 = signal.periodogram(fC40_mean, detrend="constant")
freq_f100, spec_f100 = signal.periodogram(fC100_mean, detrend="constant")
# Smoothed
Sif12_fr, Sif12_sp = signal.periodogram(Sc16_mean, detrend="constant")
Sif40_fr, Sif40_sp = signal.periodogram(Sc40_mean, detrend="constant")
Sif100_fr, Sif100_sp = signal.periodogram(Sc100_mean, detrend="constant")

# The acutal plot
fig, axs = plt.subplots(4, 1, figsize=(12, 8))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
# Plot of Original Data
axs[0].plot(np.log(Or_fr), np.log(Sp_fr), c = "black")
axs[0].set_ylim([-14, 6])
axs[0].text(-6, -13,'PSD Mean Untransformed Climate')
axs[0].set_ylabel('Log Power')
# Butterworth Plot
axs[1].plot(np.log(Bfc16_fr), np.log(Bfc16_sp))
axs[1].plot(np.log(Bfc40_fr), np.log(Bfc40_sp))
axs[1].plot(np.log(Bfc100_fr), np.log(Bfc100_sp))
axs[1].set_ylim([-14, 6])
axs[1].text(-6, -13,'PSD Mean Butterworth Climate')
axs[1].set_ylabel('Log Power')
# Fourier Plot
axs[2].plot(np.log(freq_f12), np.log(spec_f12))
axs[2].plot(np.log(freq_f40), np.log(spec_f40))
axs[2].plot(np.log(freq_f100), np.log(spec_f100))
axs[2].text(-6, -72,'PSD Mean Fourier Climate')
axs[2].set_ylabel('Log Power')
# Smooth Interpolated Plot
axs[3].plot(np.log(Sif12_fr), np.log(Sif12_sp))
axs[3].plot(np.log(Sif40_fr), np.log(Sif40_sp))
axs[3].plot(np.log(Sif100_fr), np.log(Sif100_sp))
axs[3].set_ylim([-14, 6])
axs[3].text(-6, -13,'PSD Mean Smooth Interpolated Climate')
axs[3].set_xlabel("Log Frequency ")
axs[3].set_ylabel('Log Power')
# Legend
fig.legend(labels=['4 kya', '16 kya', '40 kya', '100 kya'], loc='lower left', 
           shadow=False, ncol=4, bbox_to_anchor=(0.075, 0.05))

fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/git_frequency_domain_env.png', format='png', dpi=300, transparent=False)

#########################################
# Line Plot of Summed Hill Differences  #
#########################################
dsum_hv_wide = np.zeros(len(results))
dsum_lv_wide = np.zeros(len(results))
porp = np.zeros(len(results))
for i in range(0, len(results)):
    a = MC_Properies_Functions.compare_hill_difference(results[i], [0,2], 1, False)
    b = MC_Properies_Functions.compare_hill_difference(results_100[i], [0,2], 1, False)
    c = a-b
    dsum_hv_wide[i] = sum(a[~np.isnan(a)])
    dsum_lv_wide[i] = sum(b[~np.isnan(b)])
    porp[i] = len(c[c > 0])
# Lower y-values indicates increased evenness
# 100 kya variance is consistenly more even
##########################
fig, axs = plt.subplots(1, 1, figsize=(12, 8))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
axs.plot(range(0, len(results)), dsum_hv_wide/165, c = "orange")
axs.plot(range(0, len(results_100)), dsum_lv_wide/165, c = "blue")
axs.xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs.xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs.text(0, 0.445,'Relatively More Even')
axs.text(0, 0.625,'Relatively More Uneven')
axs.set_xlabel("Thousand Years")
axs.set_ylabel('$\Delta D^{0} - D^{2}$')
axs.set_title("$\Delta D^{0} - D^{2}$ under Wide Niche and Stabalizing Interactions")
fig.legend(labels=['4 kya Variance', '100 kya Butterworth Filtered Variance'], loc='lower right', 
           shadow=False, ncol=1, bbox_to_anchor=(0.89, 0.14))
fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/delta_hill_wide.png', format='png', dpi=300, transparent=False)

fig, axs = plt.subplots(1, 1, figsize=(12, 8))
axs.plot(range(0, len(results_100)), porp/165, c = "black")
axs.set_ylim(0.25,0.75)
axs.axhline(0.5,linestyle=':', c = 'r')
axs.fill_between(range(0, len(results_100)),0.5, 1, where = porp/165 > 0.5,fc='blue', alpha=0.1)
axs.fill_between(range(0, len(results_100)),0, 0.5, where = porp/165 < 0.5,fc='green', alpha=0.1)
axs.xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs.xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs.text(0, 0.7,'Low Variance is More Even')
axs.text(390, 0.3,'High Variance is More Even')
axs.set_title("Normalized Proportion of $\Delta D^{0} - D^{2}$ between High - Low Variance")
fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/wide_d_hill_line_plot.png', format='png', dpi=300, transparent=False)


#########################################
#            Gif of anomaly             #
#########################################
from matplotlib.colors import CenteredNorm
# Coordinates for plotting
lon_in = [Clim_array.Lon[i] for i in coord_index]
lat_in = [Clim_array.Lat[i] for i in coord_index]
for s in range(0,500):
    a = MC_Properies_Functions.compare_hill_difference(results[s], [0,2], 1, False)
    b = MC_Properies_Functions.compare_hill_difference(results_100[s], [0,2], 1, False)
    fig, axs = plt.subplots(1, 3, figsize=(12, 10))
    # Loop through columns for plotting
    in_data_a = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": a})
    in_data_a = in_data_a.pivot_table(index='Lat', columns='Lon', values='dH')
    in_data_b = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": b})
    in_data_b = in_data_b.pivot_table(index='Lat', columns='Lon', values='dH')
    in_data_c = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": -(a-b)})
    in_data_c = in_data_c.pivot_table(index='Lat', columns='Lon', values='dH')
    axs[0].imshow(in_data_a, cmap='viridis', interpolation='nearest')
    axs[1].imshow(in_data_b, cmap='viridis', interpolation='nearest')
    axs[2].imshow(in_data_c, cmap='seismic', interpolation='nearest', norm=CenteredNorm())
    axs[0].set_title("Q Order Difference High Variance")
    axs[1].set_title("Q Order Difference Low Variance")
    axs[2].set_title("Q Order Difference High - Low Variance")
    fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/wide_d_hill_maps/%d_step.png' % s)

# Create Gif
import os
import re
path = "/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/wide_d_hill_maps/"
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
imageio.mimsave("/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/wide_d_hill_maps.gif", images)


#################
#     NARROW    #
#################
dsum_hv_narrow = np.zeros(len(n_results))
hv_non_z_count = np.zeros(len(n_results))
dsum_lv_narrow = np.zeros(len(n_results))
lv_non_z_count = np.zeros(len(n_results))
porp_narrow = np.zeros(len(n_results))
for i in range(0, len(n_results)):
    a = MC_Properies_Functions.compare_hill_difference(n_results[i], [0,2], 1, False)
    b = MC_Properies_Functions.compare_hill_difference(n_results_100[i], [0,2], 1, False)
    a_c = a[~np.isnan(a)]
    b_c = b[~np.isnan(b)]
    dsum_hv_narrow[i] = sum(a_c)
    dsum_lv_narrow[i] = sum(b_c)
    hv_non_z_count[i] = len(a_c)
    lv_non_z_count[i] = len(b_c)
    c = a-b
    porp_narrow[i] = len(c[c > 0])
# Lower y-values indicates increased evenness
# 100 kya variance is consistenly more even
##########################
fig, axs = plt.subplots(1, 1, figsize=(12, 8))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
axs.plot(range(0, len(n_results)), dsum_hv_narrow/hv_non_z_count, c = "orange")
axs.plot(range(0, len(n_results_100)), dsum_lv_narrow/lv_non_z_count, c = "blue")
axs.xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs.xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs.text(0, 0.245,'Relatively More Even')
axs.text(0, 0.675,'Relatively More Uneven')
axs.set_xlabel("Thousand Years")
axs.set_ylabel('$\Delta D^{0} - D^{2}$')
axs.set_title("$\Delta D^{0} - D^{2}$ under Wide Niche and Stabalizing Interactions")
fig.legend(labels=['4 kya Variance', '100 kya Butterworth Filtered Variance'], loc='lower right', 
           shadow=False, ncol=1, bbox_to_anchor=(0.89, 0.14))
fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/delta_hill_narrow_2lineplot.png', format='png', dpi=300, transparent=False)

fig, axs = plt.subplots(1, 1, figsize=(12, 8))
axs.plot(range(0, len(n_results_100)), porp_narrow/165, c = "black")
axs.set_ylim(0.0,1.0)
axs.axhline(0.5,linestyle=':', c = 'r')
axs.fill_between(range(0, len(n_results_100)),0.5, 1, where = porp_narrow/165 > 0.5,fc='blue', alpha=0.1)
axs.fill_between(range(0, len(n_results_100)),0, 0.5, where = porp_narrow/165 < 0.5,fc='green', alpha=0.1)
axs.xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs.xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs.text(0, 0.55,'Low Variance is More Even')
axs.text(390, 0.45,'High Variance is More Even')
axs.set_title("Normalized Proportion of $\Delta D^{0} - D^{2}$ between High - Low Variance")
fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/narrow_d_hill_line_plot.png', format='png', dpi=300, transparent=False)

#########################################
#            Gif of anomaly             #
#########################################
from matplotlib.colors import CenteredNorm
# Coordinates for plotting
lon_in = [Clim_array.Lon[i] for i in coord_index]
lat_in = [Clim_array.Lat[i] for i in coord_index]
for s in range(0,500):
    a = MC_Properies_Functions.compare_hill_difference(n_results[s], [0,2], 1, False)
    b = MC_Properies_Functions.compare_hill_difference(n_results_100[s], [0,2], 1, False)
    fig, axs = plt.subplots(1, 3, figsize=(12, 10))
    # Loop through columns for plotting
    in_data_a = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": a})
    in_data_a = in_data_a.pivot_table(index='Lat', columns='Lon', values='dH')
    in_data_b = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": b})
    in_data_b = in_data_b.pivot_table(index='Lat', columns='Lon', values='dH')
    in_data_c = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": -(a-b)})
    in_data_c = in_data_c.pivot_table(index='Lat', columns='Lon', values='dH')
    axs[0].imshow(in_data_a, cmap='viridis', interpolation='nearest')
    axs[1].imshow(in_data_b, cmap='viridis', interpolation='nearest')
    axs[2].imshow(in_data_c, cmap='seismic', interpolation='nearest', norm=CenteredNorm())
    axs[0].set_title("Q Order Difference High Variance")
    axs[1].set_title("Q Order Difference Low Variance")
    axs[2].set_title("Q Order Difference High - Low Variance")
    axs[0].set_ylim(26.5,-0.5)
    axs[1].set_ylim(26.5,-0.5)
    axs[2].set_ylim(26.5,-0.5)
    axs[0].set_xlim(-0.5,12.5)
    axs[1].set_xlim(-0.5,12.5)
    axs[2].set_xlim(-0.5,12.5)
    fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/narrow_d_hill_maps/%d_step.png' % s)

# Create Gif
import os
import re
path = "/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/narrow_d_hill_maps/"
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
imageio.mimsave("/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/narrow_d_hill_maps.gif", images)





