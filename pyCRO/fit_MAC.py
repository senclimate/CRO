import sys
import numpy as np

def fit_MAC(T, h, par):
    # Extract annual mean values (first element from each parameter)
    R = par[0,0]
    F1 = par[1,0]
    epsilon = par[5,0]
    F2 = par[6,0]
    b_T = par[2,0]
    c_T = par[3,0]
    d_T = par[4,0]
    B_ref = par[10,0]
    sigma_T_ref = par[8,0]
    n_g = par[15,0]

    T = T - np.mean(T)
    h = h - np.mean(h)

    if n_g == 0:  # linear
        T2 = np.mean(T**2)
        T3 = np.mean(T**3)
        T4 = np.mean(T**4)
        T5 = np.mean(T**5)
        T6 = np.mean(T**6)
        hT = np.mean(h * T)
        hT2 = np.mean(h * T**2)
        hT3 = np.mean(h * T**3)
        hT4 = np.mean(h * T**4)

        TA =  R*T2 + F1*hT  + b_T*T3 - c_T*T4 + d_T*hT2
        T2A = R*T3 + F1*hT2 + b_T*T4 - b_T*T2**2 - c_T*T5 + c_T*T3*T2 + d_T*hT3 - d_T*hT*T2
        T3A = R*T4 + F1*hT3 + b_T*T5 - b_T*T3*T2 - c_T*T6 + c_T*T3**2 + d_T*hT4 - d_T*hT*T3

        k = -(2/3)*(T3A - 1.5*T2A*T3/T2 - 3*TA*T2) / (T4 - T3**2/T2 - T2**2)
        if np.isnan(np.sqrt(-2 * TA - k * T2)):
            sigma_T = sigma_T_ref
            B = B_ref
        else:
            sigma_T = np.sqrt(-2 * TA - k * T2)
            B = -(T2A + k * T3) / (2 * sigma_T**2 * T2)


    elif n_g == 1:  # Heaviside-linear
        T2 = np.mean(T**2)
        T3 = np.mean(T**3)
        T4 = np.mean(T**4)
        T5 = np.mean(T**5)
        T6 = np.mean(T**6)
        hT = np.mean(h * T)
        hT2 = np.mean(h * T**2)
        hT3 = np.mean(h * T**3)
        hT4 = np.mean(h * T**4)

        TA =  R*T2 + F1*hT +  b_T*T3 - c_T*T4 + d_T*hT2
        T2A = R*T3 + F1*hT2 + b_T*T4 - b_T*T2**2 - c_T*T5 + c_T*T3*T2 + d_T*hT3 - d_T*hT*T2
        T3A = R*T4 + F1*hT3 + b_T*T5 - b_T*T3*T2 - c_T*T6 + c_T*T3**2 + d_T*hT4 - d_T*hT*T3

        T1p = T
        T2p = T**2
        T3p = T**3
        T4p = T**4

        T1p = T1p[T1p >= 0]
        T2p = T2p[T2p >= 0]
        T3p = T3p[T3p >= 0]
        T4p = T4p[T4p >= 0]

        T1p = np.mean(T1p)
        T2p = np.mean(T2p)
        T3p = np.mean(T3p)
        T4p = np.mean(T4p)

        k1 = (4*TA*T2 + T2A*(2*T3p/T2p - T1p*T2/T2p) - (4/3)*T3A) / (
            2*T4p - 2*T2p*T2 - 2*(T3p**2)/T2p + T2*T3p*T1p/T2p
        )
        k2 = -(T2A + T3p*k1) / (2 * T2p)

        B = k1 / k2
        sigma_T = np.sqrt(k2 / B)

        if np.isnan(B) or np.isnan(sigma_T):
            sigma_T = sigma_T_ref
            B = B_ref

    par_out = par[:]
    par_out[8, 0] = sigma_T  # index 9 in MATLAB (0-based here)
    par_out[10, 0] = B       # index 11 in MATLAB

    return par_out
