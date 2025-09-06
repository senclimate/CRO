import numpy as np
import sys

def RO_std_analytic(par):
    R_value = par['R'][0]
    F1_value = par['F1'][0]
    epsilon_value = par['epsilon'][0]
    F2_value = par['F2'][0]
    sigma_T_value = par['sigma_T'][0]
    sigma_h_value = par['sigma_h'][0]


    # Precompute useful terms
    numerator_T = ((F1_value * F2_value - epsilon_value * R_value + epsilon_value**2) * sigma_T_value**2 +
                   (F1_value**2) * sigma_h_value**2)
    denominator = 2 * (-R_value + epsilon_value) * (F1_value * F2_value - R_value * epsilon_value)
    T_std = np.sqrt(numerator_T / denominator)

    numerator_h = ((F2_value**2) * sigma_T_value**2 +
                   (F1_value * F2_value - epsilon_value * R_value + R_value**2) * sigma_h_value**2)
    h_std = np.sqrt(numerator_h / denominator)

    return T_std, h_std
