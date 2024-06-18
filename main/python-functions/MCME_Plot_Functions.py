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





