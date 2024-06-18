# Import libraries
import os
os.chdir('/home/wp288/project/MCME/python-functions')

import importlib.machinery
MCME = importlib.machinery.SourceFileLoader('MCME_Allopatric','/home/wp288/project/MCME/python-functions/MCME_Allopatric.py').load_module()
MCME_Functions = importlib.machinery.SourceFileLoader('MCME_Functions','/home/wp288/project/MCME/python-functions/MCME_Functions.py').load_module()
Climate_Import = importlib.machinery.SourceFileLoader('Climate_Import','/home/wp288/project/MCME/python-functions/Climate_Import.py').load_module()

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

##############################
#      Import base climate   #
##############################
Tmax_4kyr = Climate_Import.import_climate_data("/home/wp288/project/MCME/climate-input/Tmax.mon.4kyr.tif", "Geotif", "rasterio")

###################
# Set parameters  #
###################
S = 25  # set sarting species number
max_growth_r = 5 # set max growth rate
M = 165 # set patch number

###########################################
#  Import Habitat Continental Coordinates #
###########################################

# Metacommunity on North America
Cont_Coords = pd.read_csv("/home/wp288/project/MCME/climate-input/World_Continent_Coords.csv")

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
#       Fourier Smooth Climate    #
###################################
# Set here to 501 steps for M patches
climate_4kya = Climate_Import.return_climate_array(coords, 400, 901, Tmax_4kyr) - 273.5
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
niche_optimum_ini_4 = MCME_Functions.initial_species_niche_optimum(S, M, climate_4kya[0,:])
niche_optimum_ini_16 = MCME_Functions.initial_species_niche_optimum(S, M, fClimate_16kya[0,:])
niche_optimum_ini_40 = MCME_Functions.initial_species_niche_optimum(S, M, fClimate_40kya[0,:])
niche_optimum_ini_100 = MCME_Functions.initial_species_niche_optimum(S, M, fClimate_100kya[0,:])
nichi_ini_matrix_list = list([niche_optimum_ini_4, niche_optimum_ini_16, niche_optimum_ini_40, niche_optimum_ini_100])

# Dispersal ability
dispersal_rate_ini = [0.01]

# Niche Breadth
niche_breadth = [2]

# Specify the speciation rate
speciation_threshold = [0.5]

# Species interaction end-members
end_member = ["equal", "stabalizing"]

# Equal species interactions
# Stabalizing species interactions
alpha_matrix_equal = MCME_Functions.initalize_aij("equal", 0.0, 1.0, S, 0.3)
alpha_matrix_stabalizing = MCME_Functions.initalize_aij("stabalizing", 0.0, 1.0, S, 0.3)
alpha_matrix_list = list([alpha_matrix_equal, alpha_matrix_stabalizing])

###################################
#        Simulation Loops         #
###################################

# Save Path
save_path = "/home/wp288/project/MCME/save-directory/time-experiment/allopatric_fs/"

# Save name vector
run_num = 1
num_sims = len(climate_input_names) * len(niche_breadth) * len(speciation_threshold) * len(end_member) * len(dispersal_rate_ini) * 15
for reps in range(15):
    for t_i, t in enumerate(climate_input_names):
        for n in niche_breadth:
            for sp in speciation_threshold:
                for end_i, end_m in enumerate(end_member):
                    for disp in dispersal_rate_ini:
                        # Set save name
                        name = "Rep:" + str(reps) + "_Clim:" + str(t) + "_n:" + str(n) + "_sp:" + str(sp) + "_" + str(end_m) + "_disp:" + str(round(disp,3))

                        # Run simulation
                        results, niche_opt, alpha_matrix, phylogeny, divergence_time, patch_origin = MCME.MCME(time_in,
                                                                                           S,
                                                                                           max_growth_r,
                                                                                           sp,
                                                                                           alpha_matrix_list[end_i],
                                                                                           distance_matrix,
                                                                                           climate_input_list[t_i],
                                                                                           nichi_ini_matrix_list[t_i],
                                                                                           disp,
                                                                                           n,
                                                                                           end_m)
                        # Set full save path
                        full_save_path = save_path + "Richness_" + name + ".npz"
                        np.savez(full_save_path, *results)
                        phylo_full_save_path = save_path + "Phylo_" + name + ".npz"
                        np.savez(phylo_full_save_path, *phylogeny)
                        patch_full_save_path = save_path + "PatchOrigin_" + name + ".npz"
                        np.savez(patch_full_save_path, *patch_origin)
                        # Print stage
                        print("%d" % run_num, " of %d" % num_sims)
                        run_num += 1
