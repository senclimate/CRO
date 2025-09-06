import numpy as np
import sys
import matplotlib.pyplot as plt

def RO_solver_analytic(par, IC, N, NE, dt=0.1, saveat=1.0, savemethod="sampling", noise_custom=None):
    NT = int(round(N / dt))
    step = int(round(saveat / dt))
    ratio = saveat / dt

    ### Checking simulation setups ###
    ##################################
    print("---------------------------------------------------------------------------------")
    print("Welcome to the CRO Analytical Solver!")
    print("Ensure that the same arguments are provided as for RO_solver,\n"
          "except that 'NM' and 'EF' are not required.")
    print("---------------------------------------------------------------------------------")

    ### Check timestep settings ###
    if dt > saveat:
        raise ValueError(f"dt={dt} and saveat={saveat}: The numerical time step (dt) must be smaller than or equal to the data saving interval (saveat).")
    if not ratio.is_integer():
        raise ValueError(f"dt = {dt}, saveat = {saveat}: saveat must be divisible by dt so that saveat/dt results in an integer.")
    ###############################

    ### extracting linear and noise amplitude parameters ###
    ########################################################
    R_value = par['R'][0]
    F1_value = par['F1'][0]
    epsilon_value = par['epsilon'][0]
    F2_value = par['F2'][0]
    sigma_T_value = par['sigma_T'][0]
    sigma_h_value = par['sigma_h'][0]
    ########################################################
    ########################################################


    ### initial setups ###
    ######################
    To, ho = IC
    t = np.arange(NT) * dt
    gr = (R_value - epsilon_value) / 2
    w = np.sqrt(4 * F1_value * F2_value - (R_value + epsilon_value)**2) / 2

    T_anal_out = np.zeros((NT, NE))
    h_anal_out = np.zeros((NT, NE))

    T_det = np.exp(gr * t) * (To * np.cos(w * t) +
        To * (R_value + epsilon_value) / (2 * w) * np.sin(w * t) +
        ho * (F1_value / w) * np.sin(w * t))

    h_det = np.exp(gr * t) * (ho * np.cos(w * t) -
        ho * (R_value + epsilon_value) / (2 * w) * np.sin(w * t) -
        To * (F2_value / w) * np.sin(w * t))
    ######################
    ######################


    ### analytical calculations for each ensemble member ###
    ########################################################
    noise_out = []
    for j in range(NE):
        if noise_custom is None:
            noise_array = np.random.randn(NT - 1, 4) # It generates noises with shapes (NT, 4, NE)
        elif np.isscalar(noise_custom):
            np.random.seed(noise_custom)             # It generates noises with seed number = noise_custom
            noise_array = np.random.randn(NT - 1, 4) # Note that noises are same for all the ensemble members
        elif noise_custom.shape == np.random.randn(NT - 1, 4).shape:
            noise_array = noise_custom               # It uses user-provided noises
        else:
            raise ValueError(f"Invalid 'noise_custom' input: expected 'None', a scalar (int or float), "
                             f"or an array with shape (N/dt-1, 4) = ({NT-1}, 4).")
        noise_out.append(noise_array)
        w_T = sigma_T_value * noise_array[:, 0] / np.sqrt(dt)
        w_h = sigma_h_value * noise_array[:, 1] / np.sqrt(dt)

        T_sto = np.zeros(NT) 
        h_sto = np.zeros(NT) 

        for i in range(1, NT): 
            tau = t[:i] 
            delta = t[i] - tau
            sin = np.sin(w * delta)
            cos = np.cos(w * delta)

            T_sto[i] = np.sum(dt * np.exp(gr * delta) * (
                w_T[:i] * cos +                                        
                w_T[:i] * (R_value + epsilon_value) / (2 * w) * sin +
                w_h[:i] * (F1_value / w) * sin
            ))

            h_sto[i] = np.sum(dt * np.exp(gr * delta) * (
                w_h[:i] * cos -                                        
                w_h[:i] * (R_value + epsilon_value) / (2 * w) * sin -
                w_T[:i] * (F2_value / w) * sin
            ))

        T_anal = T_det + T_sto 
        h_anal = h_det + h_sto 

        T_anal_out[:, j] = T_anal 
        h_anal_out[:, j] = h_anal 
    noise_out = np.array(noise_out)
    noise_out = np.transpose(noise_out, (1, 2, 0))
    ########################################################
    ########################################################


    ### save data ###
    #################
    if savemethod == "sampling":
        T_anal_out = T_anal_out[::step, :]
        h_anal_out = h_anal_out[::step, :]
    elif savemethod == "mean":
        T_anal_out = T_anal_out[int(np.ceil(step/2)):NT-int(step/2),:]
        h_anal_out = h_anal_out[int(np.ceil(step/2)):NT-int(step/2),:]
        T_anal_out = T_anal_out.reshape(-1, step, NE).mean(axis=1)
        h_anal_out = h_anal_out.reshape(-1, step, NE).mean(axis=1)
        if NE > 1:
            T_anal_out = np.vstack([T_anal_out[0:1, :], T_anal_out])
            h_anal_out = np.vstack([h_anal_out[0:1, :], h_anal_out])
        elif NE == 1:
            T_anal_out = np.concatenate([[T_anal_out[0]], T_anal_out])
            h_anal_out = np.concatenate([[h_anal_out[0]], h_anal_out])
    print("All steps are successfully completed!")
    print("---------------------------------------------------------------------------------")
    #################
    #################
    return T_anal_out, h_anal_out, noise_out


