#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 22:32:57 2024

@author: wyattpetryshen
"""
import os
os.chdir('/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/python-functions')

import numpy as np
import glob
from itertools import compress
import MC_Properies_Functions
import matplotlib.pyplot as plt
import scipy.stats as st
import re

# Import files and put into large list
# String match by Experiment 
# niche = 1, 10
# speciation = 0.01
# dispersal = 0.001, 0.01, 0.1
# end member = stabalizing, equal 
# climate: 4, 16, 40 , 100 
# Glob will be better to filter of phylogeny files
files = glob.glob("/Volumes/Wyatt_SSD_A/time-experiment-results/time-experiment/fourier-climate/*.npz")

# get unique combinations without rep number
unique_files_no_rep = glob.glob("/Volumes/Wyatt_SSD_A/time-experiment-results/time-experiment/fourier-climate/Rep:0*.npz")

# strip the rep so can iterate through unique values
unique_experiments = [i.split("Rep:0_")[1] for i in unique_files_no_rep]

# Master to store alpha mean, gamma mean, alpha UC, alpha LC, gamma UC, alpha LC
master_gamma_mean = np.zeros((len(unique_experiments), 1001))
master_gamma_uc = np.zeros((len(unique_experiments), 1001))
master_gamma_lc = np.zeros((len(unique_experiments), 1001))
master_alpha_mean = np.zeros((len(unique_experiments), 1001))
master_alpha_uc = np.zeros((len(unique_experiments), 1001))
master_alpha_lc = np.zeros((len(unique_experiments), 1001))

for e in range(len(unique_experiments)):
    # string match first unique
    unique_exp_reps = [i for i in files if unique_experiments[e] in i]
    # load files
    rep_files = [np.load(i) for i in unique_exp_reps]
    # Convert npz to list of arrays
    rep_list = list()
    for i in range(len(rep_files)):
        rep_list.append([rep_files[i][k] for k in rep_files[i]])
    # mean gamma and alpha
    end_slice = len(rep_list[0]) - 1
    # create a matrix to store values
    rep_gamma = np.zeros((len(rep_list), end_slice+1))
    rep_alpha = np.zeros((len(rep_list), end_slice+1))
    
    for i in range(len(rep_list)):
        rep_gamma[i,:] = MC_Properies_Functions.gamma_diversity(rep_list[i])
        rep_alpha[i,:] = MC_Properies_Functions.alpha_richness(rep_list[i])
        
    # mean
    rep_mean_gamma = [np.mean(rep_gamma[:,cols]) for cols in range(end_slice+1)]
    rep_mean_alpha =  [np.mean(rep_alpha[:,cols]) for cols in range(end_slice+1)]
    
    # calculate confidence intervals
    # standard deviation 
    # https://rowannicholls.github.io/python/statistics/confidence_intervals.html
    rep_std_gamma = [np.std(rep_gamma[:,cols], ddof = 1) for cols in range(end_slice+1)] # ddof = 1 is sample standard deviation
    rep_std_alpha = [np.std(rep_alpha[:,cols], ddof = 1) for cols in range(end_slice+1)] # ddof = 1 is sample standard deviation
    
    alpha = 1 - 0.95
    tails = 2
    
    t_star_gamma = st.t.ppf( 1 - (alpha / tails), len(rep_gamma[:,0]) - 1)
    t_star_alpha = st.t.ppf( 1 - (alpha / tails), len(rep_alpha[:,0]) - 1)
    
    rep_gamma_ci_upper = [rep_mean_gamma[i] + t_star_gamma * rep_std_gamma[i] / np.sqrt(len(rep_gamma[:,0])) for i in range(end_slice+1)]
    rep_gamma_ci_lower = [rep_mean_gamma[i] - t_star_gamma * rep_std_gamma[i] / np.sqrt(len(rep_gamma[:,0])) for i in range(end_slice+1)]
    
    rep_alpha_ci_upper = [rep_mean_alpha[i] + t_star_alpha * rep_std_alpha[i] / np.sqrt(len(rep_alpha[:,0])) for i in range(end_slice+1)]
    rep_alpha_ci_lower = [rep_mean_alpha[i] - t_star_alpha * rep_std_alpha[i] / np.sqrt(len(rep_alpha[:,0])) for i in range(end_slice+1)]
    
    # need a master list to store alpha and gamma means and upper confidence interval and lower confidence interval 
    master_gamma_mean[e,:] = rep_mean_gamma
    master_gamma_uc[e,:] = rep_gamma_ci_upper
    master_gamma_lc[e,:] = rep_gamma_ci_lower
    
    master_alpha_mean[e,:] = rep_mean_alpha
    master_alpha_uc[e,:] = rep_alpha_ci_upper
    master_alpha_lc[e,:] = rep_alpha_ci_lower

# Get index for color/step
col_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    for i in spltName:
        if re.search(r'\b' + "4kya" + r'\b', i):
            col_index.append(4)
        if re.search(r'\b' + "16kya" + r'\b', i):
            col_index.append(16)
        if re.search(r'\b' + "40kya" + r'\b', i):
            col_index.append(40)
        if re.search(r'\b' + "100kya" + r'\b', i):
            col_index.append(100)

index_4kya = np.where(np.isin(col_index,4))
index_16kya = np.where(np.isin(col_index,16))
index_40kya = np.where(np.isin(col_index,40))
index_100kya = np.where(np.isin(col_index,100))

# Change color index to actual color values
colors_to_plot = list()
for i in range(len(col_index)):
    if col_index[i]  == 4:
        colors_to_plot.append("red")
    elif col_index[i]  == 16:
        colors_to_plot.append("blue")
    elif col_index[i] == 40:
        colors_to_plot.append("purple")
    elif col_index[i] == 100:
        colors_to_plot.append("black")



# Equal n10
equal_n10_d01_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:10" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "equal" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.1" + r'\b', i):
            val_check += 1
        if val_check == 3:
            equal_n10_d01_index.append(e)
            
equal_n10_d001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:10" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "equal" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.01" + r'\b', i):
            val_check += 1
        if val_check == 3:
            equal_n10_d001_index.append(e)
        
equal_n10_d0001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:10" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "equal" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.001" + r'\b', i):
            val_check += 1
        if val_check == 3:
            equal_n10_d0001_index.append(e)

        
fig, ax = plt.subplots(3)        
for i in equal_n10_d01_index:
    ax[0].plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color=colors_to_plot[i], linestyle='-', label = col_index[i])
    ax[0].fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    handles, labels = ax[0].get_legend_handles_labels()
    a, b = zip(*sorted(zip(handles, labels), key=lambda x: int(x[1])))
    ax[0].legend(loc = "upper left")
    ax[0].legend(a, b)
    
    ax[0].plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color=colors_to_plot[i], linestyle='-')
    ax[0].fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)
    ax[0].set_ylim([0,50])
    ax[0].grid(True, which='major', axis='both')

for i in equal_n10_d001_index:
    ax[1].plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color=colors_to_plot[i], linestyle='-', label = col_index[i])
    ax[1].fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    handles, labels = ax[1].get_legend_handles_labels()
    a, b = zip(*sorted(zip(handles, labels), key=lambda x: int(x[1])))
    ax[1].legend(loc = "upper left")
    ax[1].legend(a, b)
    
    ax[1].plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color=colors_to_plot[i], linestyle='-')
    ax[1].fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)
    ax[1].set_ylim([0,50])

for i in equal_n10_d0001_index:
    ax[2].plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color=colors_to_plot[i], linestyle='-', label = col_index[i])
    ax[2].fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    handles, labels = ax[2].get_legend_handles_labels()
    a, b = zip(*sorted(zip(handles, labels), key=lambda x: int(x[1])))
    ax[2].legend(loc = "upper left")
    ax[2].legend(a, b)
    
    ax[2].plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color=colors_to_plot[i], linestyle='-')
    ax[2].fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)
    ax[2].set_ylim([0,50])
ax[0].set_title("Fourier Smoothed Sigma 10 Equal") 

# Stabalizing n10
stabalizing_n10_d01_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:10" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "stabalizing" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.1" + r'\b', i):
            val_check += 1
        if val_check == 3:
            stabalizing_n10_d01_index.append(e)
            
stabalizing_n10_d001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:10" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "stabalizing" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.01" + r'\b', i):
            val_check += 1
        if val_check == 3:
            stabalizing_n10_d001_index.append(e)
        
stabalizing_n10_d0001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:10" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "stabalizing" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.001" + r'\b', i):
            val_check += 1
        if val_check == 3:
            stabalizing_n10_d0001_index.append(e)
        
fig, ax = plt.subplots(3)
for i in stabalizing_n10_d01_index:
    ax[0].plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax[0].fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax[0].plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax[0].fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)

for i in stabalizing_n10_d001_index:
    ax[1].plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax[1].fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax[1].plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax[1].fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)

for i in stabalizing_n10_d0001_index:
    ax[2].plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax[2].fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax[2].plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax[2].fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)
        
    
# Equal n1
equal_n1_d01_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:1" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "equal" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.1" + r'\b', i):
            val_check += 1
        if val_check == 3:
            equal_n1_d01_index.append(e)
            
equal_n1_d001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:1" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "equal" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.01" + r'\b', i):
            val_check += 1
        if val_check == 3:
            equal_n1_d001_index.append(e)
        
equal_n1_d0001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:1" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "equal" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.001" + r'\b', i):
            val_check += 1
        if val_check == 3:
            equal_n1_d0001_index.append(e)
        
fig, ax = plt.subplots()
for i in equal_n1_d01_index:
    ax.plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax.fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax.plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax.fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)

fig, ax = plt.subplots()
for i in equal_n1_d001_index:
    ax.plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax.fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax.plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax.fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)

fig, ax = plt.subplots()
for i in equal_n1_d0001_index:
    ax.plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax.fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax.plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax.fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)
        
# Stabalizing n1
stabalizing_n1_d01_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:1" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "stabalizing" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.1" + r'\b', i):
            val_check += 1
        if val_check == 3:
            stabalizing_n1_d01_index.append(e)
            
stabalizing_n1_d001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:1" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "stabalizing" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.01" + r'\b', i):
            val_check += 1
        if val_check == 3:
            stabalizing_n1_d001_index.append(e)
        
stabalizing_n1_d0001_index = list()
for e in range(len(unique_experiments)):
    spltName = unique_experiments[e].split("_")
    val_check = 0
    for i in spltName:
        if re.search(r'\b' + "n:1" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "stabalizing" + r'\b', i):
            val_check += 1
        if re.search(r'\b' + "disp:0.001" + r'\b', i):
            val_check += 1
        if val_check == 3:
            stabalizing_n1_d0001_index.append(e)
        
fig, ax = plt.subplots()
for i in stabalizing_n1_d01_index:
    ax.plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax.fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax.plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax.fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)

fig, ax = plt.subplots()
for i in stabalizing_n1_d001_index:
    ax.plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax.fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax.plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax.fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)

fig, ax = plt.subplots()
for i in stabalizing_n1_d0001_index:
    ax.plot(range(master_gamma_mean.shape[1]), master_gamma_mean[i,:], color='black', linestyle='-')
    ax.fill_between(range(master_gamma_mean.shape[1]), master_gamma_lc[i,:], master_gamma_uc[i,:], color='b', alpha=.1)
    
    ax.plot(range(master_alpha_mean.shape[1]), master_alpha_mean[i,:], color='green', linestyle='-')
    ax.fill_between(range(master_alpha_mean.shape[1]), master_alpha_lc[i,:], master_alpha_uc[i,:], color='y', alpha=.1)
        
        
        
        
    







plt.plot(range(len(rep_gamma[0,:])), rep_mean_gamma, color='black', linestyle='-')
plt.fill_between(range(len(rep_gamma[0,:])), rep_gamma_ci_lower, rep_gamma_ci_upper, color='b', alpha=.1)

plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[8,:], color='black', linestyle='-')
plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[13,:], color='blue', linestyle='-')
plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[17,:], color='green', linestyle='-')

plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[8,:], color='black', linestyle='-')
plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[13,:], color='blue', linestyle='-')
plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[17,:], color='green', linestyle='-')


plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[12,:], color='black', linestyle='-')
plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[15,:], color='blue', linestyle='-')

plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[12,:], color='black', linestyle='-')
plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[15,:], color='blue', linestyle='-')

plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[9,:], color='black', linestyle='-')
plt.plot(range(len(master_gamma_mean[0,:])), master_gamma_mean[10,:], color='blue', linestyle='-')

plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[9,:], color='black', linestyle='-')
plt.plot(range(len(master_alpha_mean[0,:])), master_alpha_mean[10,:], color='blue', linestyle='-')



# Filter experiments 
# By niche
index_n_10 = ["n_10" in i for i in files] # index 10
index_n_1 = [True if i == False else False for i in index_n_10] # index 1

n1_files = list(compress(files, index_n_1))
n10_files = list(compress(files, index_n_10))

# Filter by end member
#  n1
index_n_1_equal = ["equal" in i for i in n1_files] # n 1 equal
index_n_1_stabalizing = ["stabalizing" in i for i in n1_files] # n 1 stable
equal_n1_files = list(compress(n1_files, index_n_1_equal))
stabalizing_n1_files = list(compress(n1_files, index_n_1_stabalizing))

# n10
index_n_10_equal = ["equal" in i for i in n10_files] # n 10 equal
index_n_10_stabalizing = ["stabalizing" in i for i in n10_files] # n 10 stable
equal_n10_files = list(compress(n10_files, index_n_10_equal))
stabalizing_n10_files = list(compress(n10_files, index_n_10_stabalizing))
 
# Import for plotting
equal_n1 = [np.load(i) for i in equal_n1_files]
stabalizing_n1 = [np.load(i) for i in stabalizing_n1_files]

# Convert npz to list of arrays
equal_n1_list = list()
for i in range(len(equal_n1)):
    equal_n1_list.append([equal_n1[i][k] for k in equal_n1[i]])

stabalizing_n1_list = list()
for i in range(len(stabalizing_n1_list)):
    stabalizing_n1_list.append([stabalizing_n1[i][k] for k in stabalizing_n1[i]])


# Calculate alpha and gamma

# Sigma 10 
sigma10_gamma_equal = list()
sigma10_alpha_equal = list()

end_slice = len(MCME_niche10_equal[0]) - 1

for i, output in enumerate(MCME_niche10_equal):
    t_gamma = MC_Properies_Functions.gamma_diversity(output)[-1]
    t_alpha = MC_Properies_Functions.alpha_richness(output)[-1]
    sigma10_gamma.append(t_gamma)
    sigma10_alpha.append(t_alpha)
    

# Sigma 1
n1_gamma_equal = list()
n1_alpha_equal = list()

end_slice = len(equal_n1_list[0]) - 1
for i, output in enumerate(equal_n1_list):
    t_gamma = MC_Properies_Functions.gamma_diversity(output)[-1]
    t_alpha = MC_Properies_Functions.alpha_richness(output)[-1]
    n1_gamma_equal.append(t_gamma)
    n1_alpha_equal.append(t_alpha)
    
# Plot of gamma and alpha diversity
x = np.arange(0,len(n1_gamma_equal))
plt.plot(x, n1_gamma_equal, label = "gamma")
plt.plot(x, n1_alpha_equal, c = "green", label = "alpha")
plt.legend(loc="upper left")

    
    
    
    