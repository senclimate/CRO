import numpy as np
from RO_fitting_LR import *
from RO_fitting_MLE import *
from regress_std import *
from func_default_fitting_method import *
from RO_fitting_MLE_white import *
from RO_fitting_MLE_red import *
import sys

def RO_fitting(T, h, par_option_T, par_option_h, par_option_noise, method_fitting=None, dt=None):
    fitting_option_red = "ARn"   # "LR" or "AR1" or "ARn"
                                 # this option is defined internally in RO_fitting 
                                 # and then passed to RO_fitting_LR
                                 # meaning it is not exposed as a top-level fitting option to the user.

    ### Checking fitting setups ###
    ###############################
    print("---------------------------------------------------------------------------------")
    print("Welcome to CRO Fitting! Your fitting setups:")
    print("---------------------------------------------------------------------------------")

    if dt == None:
         dt = 1.0
         print(f" - Data time step is not given, defaulting to: dt = 1.0 months.")
    print(f" - Time series length: N = len(T)*dt = {len(T)*dt} months.")

    print(f" - Prescribed terms: {par_option_T}. \n"
      f"                     {par_option_h}. \n"
      "   0 - Do not prescribe. \n"
      "   1 - Prescribe only the annual mean. \n"
      "   3 - Prescribe the annual mean and annual seasonality. \n"
      "   5 - Prescribe the annual mean, annual seasonality, and semi-annual seasonality.")

    if par_option_noise['T'] != par_option_noise['h']:
        raise ValueError(f"par_option_noise['T'] = {par_option_noise['T']}, "
                         f"par_option_noise['h'] = {par_option_noise['h']}\n"
                          "Fitting methods for T and h equations should be the same.")

    print(f" - Noise options: {par_option_noise}.")
    if (par_option_noise['T'] == 'red') and (par_option_noise['h'] == 'red'):
        print(f" - Fitting method for T and h red noises: {fitting_option_red}.\n"
              f"   This option is defined internally within RO_fitting.py.\n"
              f"   Options available are: LR or AR1 or ARn.")
    ##################################
    ##################################



    # Convert parameter options
    par_option_T_vals = list(par_option_T.values())
    par_option_h_vals = list(par_option_h.values())
    noise_keys = ['T', 'h', 'T_type']
    
    noise_map = {
        'white': 1,
        'red': 0,
        'additive': 2,
        'multiplicative': 0,
        'multiplicative-H': 1
    }
    par_option_noise_array = [noise_map[str(par_option_noise[k])] for k in noise_keys]


    print(f" - Fitting method for T and h main equations: {method_fitting}.")
    if method_fitting == None:
        method_fitting = func_default_fitting_method(par_option_T_vals, par_option_h_vals, par_option_noise_array)


    # Perform fitting
    if method_fitting == "LR-F":
        par = RO_fitting_LR(T, h, par_option_T_vals, par_option_h_vals, par_option_noise_array, dt, "F", "LR", fitting_option_red)
    elif method_fitting == "LR-C":
        par = RO_fitting_LR(T, h, par_option_T_vals, par_option_h_vals, par_option_noise_array, dt, "C", "LR", fitting_option_red)
    elif method_fitting == "LR-F-MAC":
        par = RO_fitting_LR(T, h, par_option_T_vals, par_option_h_vals, par_option_noise_array, dt, "F", "MAC", fitting_option_red)
    elif method_fitting == "MLE":
        par = RO_fitting_MLE(T, h, par_option_T_vals, par_option_h_vals, par_option_noise_array, dt)
    else:
        raise ValueError(f"Unknown fitting method: {method_fitting}")

    if method_fitting in ("LR-F", "LR-F-MAC", "MLE"):
        if par[13][0] == 1: # Do Ito to Strantonovich Conversion only when white noise is used
            if par[15][0] == 0:   # Ito to Strantonovich Conversion for multiplicative noise
                par[0][0] = par[0][0] - 0.5 * (par[8][0] * par[10][0])**2
            elif par[15][0] == 1: # Ito to Strantonovich Conversion for Heaviside multiplicative noise
                par[0][0] = par[0][0] - 0.25 * (par[8][0] * par[10][0])**2
            elif par[15][0] == 2: # No correction needed for additive noise
                par[0][0] = par[0][0]


    ### Round negligible seasonal amplitudes of linear/nonlinear  ###
    ### parameters to zero.                                       ###
    for i in range(0,8):
        if par[i][1] < 10**(-4):
            par[i][1] = 0.0
            par[i][2] = 0.0
        if par[i][3] < 10**(-4):
            par[i][3] = 0.0
            par[i][4] = 0.0
    #################################################################
    #################################################################
   

    Arr = []
    for i in range(0, 13):
        arr = par[i]
        arr = arr[arr != 0]
        Arr.append(arr)


    param_keys = [
    "R", "F1", "F2", "epsilon", "b_T", "c_T", "d_T", "b_h",
    "sigma_T", "sigma_h", "B", "m_T", "m_h", "n_T", "n_h", "n_g"]
    param_values = [Arr[0], Arr[1], Arr[5], Arr[6], Arr[2], Arr[3], Arr[4], Arr[7], 
                    Arr[8], Arr[9], Arr[10], Arr[11], Arr[12], 
                    [int(par[13, 0])], [int(par[14, 0])], [int(par[15, 0])]]

    par = dict(zip(param_keys, param_values))
    par = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in par.items()}


    ### Displaying fitting results ###
    ##################################
    print("---------------------------------------------------------------------------------")
    print("All steps are successfully completed!")
    print("---------------------------------------------------------------------------------")

    return par

