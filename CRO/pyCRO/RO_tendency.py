import numpy as np
import sys

def RO_tendency(par, T, h, xi_T, xi_h, EF):
    # Unpack parameters
    R, F1, F2, epsilon = par[0:4]
    b_T, c_T, d_T, b_h = par[4:8]
    sigma_T, sigma_h, B = par[8:11]
    m_T, m_h = par[11:13]
    n_T, n_h, n_g = map(int, par[13:16])  # Ensure integer indexing

    # Unpack external forcing
    E_T, E_h = EF


    # Noise factor
    if n_g == 0:
        g_T_factor = B * T
    elif n_g == 1:
        h_T_factor = max(T, 0)  # Heaviside(T)
        g_T_factor = B * h_T_factor
    elif n_g == 2:
        g_T_factor = 0.0
    else:
        raise ValueError("Invalid value for n_g")


    # dT/dt
    if n_T == 0:  # red noise
        f_T = R*T + F1*h + b_T*(T**2) - c_T*(T**3) + d_T*T*h + sigma_T*(1 + g_T_factor)*xi_T + E_T
        g_T = 0.0
    elif n_T == 1:  # white noise
        f_T = R*T + F1*h + b_T*(T**2) - c_T*(T**3) + d_T*T*h + E_T
        g_T = sigma_T * (1 + g_T_factor)
    else:
        raise ValueError("Invalid value for n_T")


    # dh/dt
    if n_h == 0:  # red noise
        f_h = -F2*T - epsilon*h - b_h*(T**2) + sigma_h*xi_h + E_h
        g_h = 0.0
    elif n_h == 1:  # white noise
        f_h = -F2*T - epsilon*h - b_h*(T**2) + E_h
        g_h = sigma_h
    else:
        raise ValueError("Invalid value for n_h")

    # d(xi_T)/dt
    if n_T == 0:
        f_xi_T = -m_T * xi_T
        g_xi_T = np.sqrt(2 * m_T)
    else:
        f_xi_T = np.nan
        g_xi_T = np.nan

    # d(xi_h)/dt
    if n_h == 0:
        f_xi_h = -m_h * xi_h
        g_xi_h = np.sqrt(2 * m_h)
    else:
        f_xi_h = np.nan
        g_xi_h = np.nan

    return f_T, g_T, f_h, g_h, f_xi_T, g_xi_T, f_xi_h, g_xi_h
