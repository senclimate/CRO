from RO_fitting_MLE_white import *
from RO_fitting_MLE_red import *

def RO_fitting_MLE(T, h, par_option_T, par_option_h, par_option_noise, dt):
    n_T = par_option_noise[0]
    n_h = par_option_noise[1]

    if n_T == 1 and n_h == 1:
        par = RO_fitting_MLE_white(T, h, par_option_T, par_option_h, par_option_noise, dt)
    elif n_T == 0 and n_h == 0:
        par = RO_fitting_MLE_red(T, h, par_option_T, par_option_h, par_option_noise, dt)
    else:
        raise ValueError("Error: Mixed noise types not supported by MLE")
    
    return par
