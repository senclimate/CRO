import numpy as np
import sys

### This module generates time-evolving RO parameters and external forcing for use in the RO solver. ###
########################################################################################################
def par_processing(par_in, dt, NT): 
    num_params = len(par_in)
    par_out = np.ones((NT, num_params))

    ro_time = np.linspace(0, (NT - 1) * dt, NT)
   
    for i, par in enumerate(par_in):
        params = par_in[par]
        if np.isscalar(params): # if 'n_T'/'n_h'/'n_g' is given as a scalar instead of list, convert it to list
            params = [params]   
        
        if (par == 'n_T') or (par == 'n_h') or (par == 'n_g'):
            par_out[:, i] = params[0]
        else:
            w = 2 * np.pi / 12
            par_out[:, i] = (params[0] + params[1] * np.sin(w * ro_time + params[2]) +
                             params[3] * np.sin(2 * w * ro_time + params[4]))
    return par_out # array of shape (NT, 16) or (NT, 2) containing time-evolving parameters
########################################################################################################
########################################################################################################
