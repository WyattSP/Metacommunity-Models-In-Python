#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 18:30:46 2024

@author: wyattpetryshen
"""

import numpy as np
import scipy as scipy
import matplotlib.pyplot as plt
from scipy import integrate         

# Let's first model a simple predator-prey system
# time
t = np.arange(0,1000,1)
# populations
N0 = np.array([100, 100])
# grow rates
r_m = np.random.uniform(0, 1, 2)

# Interaction coefficients
alpha_m = np.random.uniform(0, 1, (2, 2))
np.fill_diagonal(alpha_m, 1.0)
 
 
def N2_lokta(N, t = 0):
    x_out = np.array([ r_m[0] * (1 / (1 + alpha_m[0,1] * N[1])),
                      r_m[1] * (1 / (1 + alpha_m[1,0] * N[0])) ])
    return(x_out)
    
# integrate
N, infodict = integrate.odeint(N2_lokta, N0, t, full_output=True)


N1, N2 = N.T
plt.plot(t, N1, 'r-', label='N1')
plt.plot(t, N2  , 'b-', label='N2')
plt.grid()
plt.xlabel('time')
plt.ylabel('population')


            density = N @ alpha_matrix

            # Species density change
            lambda_vector = r_ix * N * (1.0 / (1.0 + density))
            lambda_vector[lambda_vector < 0.0] = 0.0
            
results[0].shape
alpha_matrix[0].shape

results[0][0,:]


# rename the alpha_matrix[0] to alpha_matrix; may through error if definition alpha_matrix does not update with main loop
def density_dependent_interactions(t, Y, alpha_in, r_growth):
    # Can only input sinlge habitat patch 
    density = Y @ alpha_in
    lambda_vector = r_growth  * (1.0 / (1.0 + density))
    return(lambda_vector)
    

r_g = [1,1,1,1,1,0.5,0.5,0.5,0.5,0.5]

N = np.zeros(results[0].shape)
for i in range(results[0].shape[0]):
    # Run RK45
    X = integrate.solve_ivp(density_dependent_interactions, [0,100], results[0][i,:], method="RK45", dense_output = True, args = (alpha_matrix[0], r_g, ))
    # Save last time slice to matrix
    N[i,:] = X.y[:,-1]

#X = integrate.solve_ivp(density_dependent_interactions, [0,1000], results[0][0,:], method="RK45",dense_output = True)


for s in range(X.y.shape[0]):
    plt.plot(X.t, X.y[s])

plt.plot(N)


len(N)






