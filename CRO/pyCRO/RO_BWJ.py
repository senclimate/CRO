import numpy as np

def RO_BWJ(par):
    # Extract annual mean values
    R_value = par['R'][0]
    F1_value = par['F1'][0]
    epsilon_value = par['epsilon'][0]
    F2_value = par['F2'][0]

    # Calculate growth rate and frequency
    gr = (R_value - epsilon_value) / 2
    w = np.sqrt(4 * F1_value * F2_value - (R_value + epsilon_value)**2) / 2

    # Return complex index
    return gr + 1j * w
