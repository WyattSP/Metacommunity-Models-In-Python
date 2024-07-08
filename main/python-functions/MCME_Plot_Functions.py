#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 07:39:16 2024

@author: wyattpetryshen
"""

import matplotlib.pyplot as plt
import numpy as np
#from treelib import Node, Tree
import seaborn as sns
import pandas as pd
import math

def mcme_plot(population_list, patch):
    # This is technically the incorrect number of generations
    if len(patch) == 1:
        # Function of selected save intervals
        x = np.arange(0,len(population_list))
        # Stack save points into array
        pop_arr = np.dstack(population_list)
        # Plot dynamics for desired patch
        for i in np.arange(0, pop_arr.shape[1]):
            plt.plot(x, pop_arr[patch,i,:])
        plt.show()
    else:
        pop_arr = np.dstack(population_list)
        plot_rows = pop_arr.shape[0]
        fig, axarr = plt.subplots(plot_rows, 1, sharex=True, figsize=(12,6))

        x = np.arange(0,len(population_list))

        count = 0
        for ax in axarr.ravel():
            for j in np.arange(0, pop_arr.shape[1]):
                ax.plot(x, pop_arr[count,j,:])
            count += 1
        fig.tight_layout()

def species_richness_plot(N_save, time_step_length):
    # Get x-axis length
    x_axis_len_values =  np.arange(0, len(N_save)) * time_step_length
    # Number of species per time step
    total_living_species = [len(np.where(i.any(axis=0))[0]) for i in N_save]
    # Get orgination
    # origination = [total_living_species[i+1] - total_living_species[i] for i in range(len(total_living_species)-1)]

    # Plot
    fig, ax = plt.subplots()
    ax.plot(x_axis_len_values, total_living_species, c = "blue")
    ax.set_xlabel('Simulation Time')
    ax.set_ylabel('Change in Total Species Per Time Step')

    plt.show()

def plot_origination_extinction(N_save, plot = False):
    generations = len(N_save)

    # Origination
    origination = np.empty((generations, 2))
    for i in np.arange(1, generations):
        next_o = N_save[i].shape[1] - N_save[i-1].shape[1]
        origination[i, 0] = i
        origination[i, 1] = next_o

    # Extinction
    extinction = np.empty((generations, 2))
    for i in np.arange(1, generations):
        next_e = (len((N_save[i] != 0).any(axis=0)) - sum((N_save[i] != 0).any(axis=0))) - (len((N_save[i-1] != 0).any(axis=0)) - sum((N_save[i-1] != 0).any(axis=0)))
        extinction[i, 0] = i
        extinction[i, 1] = next_e

    if plot == True:
        fig, ax = plt.subplots()
        ax.plot(origination[:,0], origination[:,1], c = "blue")
        ax.plot(extinction[:,0], extinction[:,1], c = "red")
        ax.set_title("Origination & Extinction")
        ax.set_xlabel('Simulation Time')
        ax.set_ylabel('Species')
        plt.show()

    return(origination, extinction)


def plot_phylogeny(phylogeny, divergence_time):
    # phylogeny[index] list object
    # pyhlogeny[index][0,:] ancestor
    # phylogeny[index][1,:] descendent

    #tree = Tree()


    #str(phylogeny[0][1])

    # first node
    #tree.create_node(str(phylogeny[0][0][0]),str(phylogeny[0][0][0]))
    # next node
    #tree.create_node(str(phylogeny[0][1][0]),str(phylogeny[0][1][0]), parent = str(phylogeny[0][0][0]))

    #
    #tree.show()
    return()

def get_richness(N_values):
    shape = N_values.shape
    out_richness= np.zeros(shape[0])
    for n_rows in range(shape[0]):
        val = sum(N_values[n_rows,:] != 0)
        out_richness[n_rows] = val
    return(out_richness)

def plot_richness_map(N_values, coords, coord_index, return_pivot = False):

    lon_in = [coords.Lon[i] for i in coord_index]
    lat_in = [coords.Lat[i] for i in coord_index]
    r_in = get_richness(N_values)

    in_data = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "R": r_in})

    in_data = in_data.pivot_table(index='Lat', columns='Lon', values='R')

    sns.heatmap(in_data, annot=False, fmt="g", cmap='viridis')

    if return_pivot == True:
        return(in_data)

# Comparative Hill Number Methodology
# Chao, A., Ricotta, C. 2019. Quantifying evenness and linking it to diversity, beta diversity, and similarity
# Ecology, 100(12): e02852.


# Next step is calculate angles for all hill number values from end of simulation step and plot on map using below function
# Bascially just need to create a function the creats a patch x q order - 1 matrix
# Iterate over all q orders for each path using patch_angle diff and put into output matrix


def plot_hill_number_maps(Hill_Numbers, coords, coord_index, return_pivot = False):
    # Coordinates for plotting
    lon_in = [coords.Lon[i] for i in coord_index]
    lat_in = [coords.Lat[i] for i in coord_index]
    # Loop through columns for plotting
    fig, axs = plt.subplots(1, Hill_Numbers.shape[1], figsize=(12, 8))
    #min_vals = np.zeros(Hill_Numbers.shape[1])
    #max_vals = np.zeros(Hill_Numbers.shape[1])
    for q_values in np.arange(0,Hill_Numbers.shape[1],1):
        plot_vals = Hill_Numbers[:,q_values]
        in_data = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "H": plot_vals})
        in_data = in_data.pivot_table(index='Lat', columns='Lon', values='H')
        axs[q_values].imshow(in_data, cmap='viridis', interpolation='nearest')
        axs[q_values].set_title("Order: %d" % q_values)
        #min_vals[q_values] = min(plot_vals)
        #max_vals[q_values] = max(plot_vals)
    # Find the min and max of all colors for use in setting the color scale.
    #norm = colors.Normalize(vmin=vmin, vmax=vmax)
    #fig.colorbar(axs[0], ax=axs, orientation='horizontal', fraction=.1)
    if return_pivot == True:
        return(in_data)

def plot_hill_difference_map(delta_hill, coords, coord_index, return_pivot = False):
    # Coordinates for plotting
    lon_in = [coords.Lon[i] for i in coord_index]
    lat_in = [coords.Lat[i] for i in coord_index]
    # Loop through columns for plotting
    fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    in_data = pd.DataFrame({"Lon": lon_in, "Lat": lat_in, "dH": delta_hill})
    in_data = in_data.pivot_table(index='Lat', columns='Lon', values='dH')
    axs.imshow(in_data, cmap='viridis', interpolation='nearest')
    axs.set_title("Q Order Difference Map")
    if return_pivot == True:
        return(in_data)
    
    
    
    
    
    


