import numpy as np
import sys
from par_processing import *
from RO_integral import *

### Shapes of RO parameters and external forcing ###
####################################################
par_shapes = ({'R': 5, 'F1': 5, 'F2': 5, 'epsilon': 5, 
               'b_T': 5, 'c_T': 5, 'd_T': 5, 'b_h': 5, 
               'sigma_T': 5, 'sigma_h': 5, 'B': 5, 
               'm_T': 5, 'm_h': 5, 
               'n_T': 1, 'n_h': 1, 'n_g': 1})
EF_shapes = {'E_T': 5, 'E_h': 5}
####################################################
####################################################


### Checking input par and EF have expected shapes ###
######################################################
def has_same_shape(other_dict, par_shapes):
    for key in par_shapes:
        if not isinstance(other_dict[key], list):
            other_dict[key] = [other_dict[key]]
        key_len = len(other_dict[key])
        #print(key, key_len)
        if not key in ['n_T', 'n_h', 'n_g']:
            if key_len not in [0, 1, 3, 5]:
                raise ValueError(f"Shape mismatch for key '{key}': expected length 0, 1, 3, or 5, but got {key_len}.")
            elif key_len in [0, 1, 3, 5]:
                other_dict[key] = other_dict[key] + [0.0] * (5 - len(other_dict[key]))
        else:
            if key_len != 1:
                raise ValueError(f"Shape mismatch for key '{key}': expected length 1, but got {key_len}.")
    return True
######################################################
######################################################


### RO solver ###
#################
def RO_solver(par, IC, N, NE, NM="EH", dt=0.1, saveat=1.0, savemethod="sampling", EF=None, noise_custom=None):
    NT = int(round(N / dt))
    step = int(round(saveat / dt))
    ratio = saveat / dt

    ### Checking simulation setups ###
    print("---------------------------------------------------------------------------------")
    print("Welcome to the CRO Solver! Your simulation setup is as follows:")
    print("---------------------------------------------------------------------------------")

    print(f" - Total simulation length: N = {N} months")
    print(f" - Number of ensemble members: NE = {NE}")
    print(f" - Numerical integration time step: dt = {dt} months (default: 0.1)")
    print(f" - Data output interval: saveat = {saveat} months (default: 1.0)")


    ### Check timestep settings ###
    if dt > saveat:
        raise ValueError(f"dt = {dt}, saveat = {saveat}: The numerical time step (dt) must be less than or equal to the data saving interval (saveat).")
    if not ratio.is_integer():
        raise ValueError(f"dt = {dt}, saveat = {saveat}: saveat must be divisible by dt so that saveat/dt is an integer.")
    ###############################

    ### Check Initial conditions ###
    if np.array(IC).shape != (2,):
        raise ValueError(f"Invalid 'IC' input: expected an array with shape (2,), but got {IC}.")
    else:
        print(f" - Initial conditions: IC = [T0, h0] = {IC}")
    ################################

    ### Check input parameter shapes ###
    if not has_same_shape(par, par_shapes):
        raise ValueError("Input parameters do not have the expected shape.")
    else:
        print(f" - Input parameters have the expected shapes.")
    ####################################

    ### Check noise parameters ###
    if np.array(par['n_T']).item() == 1:
        print(" - 'n_T' = 1: White noise forcing in the T equation; 'm_T' is ignored.")
    elif np.array(par['n_T']).item() == 0:
        print(f" - 'n_T' = 0: Red noise forcing in the T equation with m_T = {par['m_T']}.")
        if (np.all(np.array(par['m_T']) == 0)):
            raise ValueError(f"'m_T' needs to have at least one non-zero argument for red noise forcing.")
    else:
        raise ValueError(f"n_T = {np.array(par['n_T']).item()} is an invalid option: it has to be either '1' (red noise) or '0' (white noise).")

    if np.array(par['n_h']).item() == 1:
        print(" - 'n_h' = 1: White noise forcing in the h equation; 'm_h' is ignored.")
    elif np.array(par['n_h']).item() == 0:
        print(f" - 'n_h' = 0: Red noise forcing in the h equation with m_h = {par['m_h']}.")
        if (np.all(np.array(par['m_h']) == 0)):
            raise ValueError(f"'m_h' needs to have at least one non-zero argument for red noise forcing.")
    else:
        raise ValueError(f"n_h = {np.array(par['n_h']).item()} is an invalid option: it has to be either '1' (red noise) or '0' (white noise).")   

    if np.array(par['n_g']).item() == 2:
        print(" - 'n_g' = 2: Additive noise is used in the T equation; 'B' is ignored.")
    elif np.array(par['n_g']).item() == 0:
        print(" - 'n_g' = 0: Multiplicative noise (1 + B*T) is used in the T equation.")
    elif np.array(par['n_g']).item() == 1:
        print(" - 'n_g' = 1: Multiplicative noise with the heaviside function (1 + B*H*T) is used in the T equation.")
    else:
        raise ValueError(f"n_g = {np.array(par['n_g']).item()} is an invalid option. "
                          "It must be one of the following:\n"
                          "  0 — multiplicative\n"
                          "  1 — multiplicative + Heaviside\n"
                          "  2 — additive")

    ##############################

    ### Check numerical integration method ###
    if NM == "EM":
        print(f" - Numerical integration method: NM = '{NM}' (Euler–Maruyama method)")
        print("   The default method is NM = 'EH' (Euler–Heun method)")
    elif NM == "EH":
        print(f" - Numerical integration method: NM = '{NM}' (Euler–Heun method; default)")
    else:
        raise ValueError(
            f"{NM} is an invalid NM input. Expected 'EM' (Euler–Maruyama method) or "
            f"'EH' (Euler–Heun method)"
    )
    ##########################################

    ### Check data saving method ###
    if savemethod == "sampling":
        print(f" - Data saving method: savemethod = {savemethod} (default)")
    elif savemethod == "mean":
        print(f" - Data saving method: savemethod = {savemethod} (default is 'sampling')")
    else:
        raise ValueError(f"{savemethod} is an invalid savemethod: it has to be either 'sampling' or 'mean'.")
    ################################

    ### Check external forcing ###
    if EF is None:
        print(" - External forcing is not given, therefore using\n   EF = {'E_T': [0.0, 0.0, 0.0, 0.0, 0.0], 'E_h': [0.0, 0.0, 0.0, 0.0, 0.0]}.")
        EF = {'E_T': [0.0, 0.0, 0.0, 0.0, 0.0], 'E_h': [0.0, 0.0, 0.0, 0.0, 0.0]}
    elif has_same_shape(EF, EF_shapes):
        print(f" - External forcing is given as\n  EF = {EF}.")
    else:
        raise ValueError(f"- External forcing does not match the shape of {EF_shapes}.")
    ##############################

    ### Check noise ###
    if noise_custom is None:
        print(f" - noise_custom = {noise_custom}: System-generated noise is used and changes at every run.")
    elif np.isscalar(noise_custom):
        print(f" - noise_custom = {noise_custom}: Seed number {noise_custom} is used to generate\n"
              f"   the noises shared for all ensemble members.")

    elif noise_custom.shape == (NT-1,4):
        print(f" - You provided the customized noises with shapes (N/dt-1, 4) = ({NT-1}, 4).")
    else:
        raise ValueError(f"Invalid 'noise_custom' input: expected 'None', a scalar (int or float), "
                         f"or an array with shape (N/dt-1, 4) = ({NT-1}, 4).")
    ###################
    print("---------------------------------------------------------------------------------")
    ##################################
    ##################################



    ### preprocessing the input parameters ###
    ##########################################
    T_out = np.zeros((NT, NE))                   # reformatting into shapes (NT, NE)
    h_out = np.zeros((NT, NE))                   # reformatting into shapes (NT, NE)
    noise_out = np.full((NT - 1, 4, NE), np.nan) # creating dummy noise arrays with shapes (NT, 4, NE)
                                                 # (NT, 0, NE) -> for main T equation
                                                 # (NT, 1, NE) -> for main h equation
                                                 # (NT, 2, NE) -> for noise T equation
                                                 # (NT, 3, NE) -> for noise h equation
    if (NM == "EM") and (par['n_T'][0] == 1): # Do Ito to Strantonovich Conversion
                                              # when Euler-Maruya scheme is used with white noise forcing
        if par['n_g'][0] == 0:   # Ito to Strantonovich Conversion for multiplicative noise
            par['R'][0] = par['R'][0] + 0.5 * (par['sigma_T'][0] * par['B'][0])**2
        elif par['n_g'][0] == 1: # Ito to Strantonovich Conversion for Heaviside multiplicative noise
            par['R'][0] = par['R'][0] + 0.25 * (par['sigma_T'][0] * par['B'][0])**2
        elif par['n_g'][0] == 2: # No correction needed for additive noise
            par['R'][0] = par['R'][0]
    par_in = par_processing(par, dt, NT)         # returns (NT, 16)
    EF_in = par_processing(EF, dt, NT)           # returns (NT, 2)
    ##########################################
    ##########################################


    ### Numerical integration ###
    ##############################
    print(f"Numerical integration is being performed:")
    print("---------------------------------------------------------------------------------")
    for i in range(NE):
        print(f" - Ensemble member: {i+1}/{NE}")
        if noise_custom is None:
            noise_array = np.random.randn(NT - 1, 4) # It generates noises with shapes (NT, 4, NE)
        elif np.isscalar(noise_custom):
            np.random.seed(noise_custom)             # It generates noises with seed number = noise_custom
            noise_array = np.random.randn(NT - 1, 4) # Note that noises are same for all the ensemble members 
        elif noise_custom.shape == noise_out[:, :, 0].shape:
            noise_array = noise_custom               # It uses user-provided noises
             
        T, h, noise = RO_integral(par_in, EF_in, NM, NT, dt, IC, noise_array) # numerical integration
        T_out[:, i] = T
        h_out[:, i] = h
        noise_out[:, :, i] = noise
    print("---------------------------------------------------------------------------------")
    ##############################
    ##############################


    ### save data ###
    #################
    if (savemethod == "sampling") or (savemethod is None):
        T_out = T_out[::step, :]
        h_out = h_out[::step, :]
    elif savemethod == "mean":
        #T_out = T_out[:NT-NT%step,:]
        #h_out = h_out[:NT-NT%step,:]
        #T_out = T_out.reshape(step, -1, NE).mean(axis=0)
        #h_out = h_out.reshape(step, -1, NE).mean(axis=0)
        T_out = T_out[int(np.ceil(step/2)):NT-int(step/2),:]
        h_out = h_out[int(np.ceil(step/2)):NT-int(step/2),:]
        T_out = T_out.reshape(-1, step, NE).mean(axis=1)
        h_out = h_out.reshape(-1, step, NE).mean(axis=1)
        if NE > 1:
            T_out = np.vstack([T_out[0:1, :], T_out])
            h_out = np.vstack([h_out[0:1, :], h_out])
        elif NE == 1:
            T_out = np.concatenate([[T_out[0]], T_out])
            h_out = np.concatenate([[h_out[0]], h_out])
    #################
    #################
    print("All steps are successfully completed!")
    print("---------------------------------------------------------------------------------")
 
    return T_out, h_out, noise_out
#################
#################

