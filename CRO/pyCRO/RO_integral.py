import numpy as np
import sys
from RO_tendency import *
from EM_scheme import *

### RO integrater (Euler-Maruyama and Euler-Huen) ###
#####################################################
def RO_integral(par, EF, NM, NT, dt, IC, noise_in):
    T = np.ones(NT)
    h = np.ones(NT)
    xi_T = np.ones(NT)
    xi_h = np.ones(NT)

    T[0] = IC[0]
    h[0] = IC[1]
    xi_T[0] = 0.0
    xi_h[0] = 0.0

    noise_out = np.full((NT - 1, 4), np.nan)



    for i in range(NT - 1):
        f_T, g_T, f_h, g_h, f_xi_T, g_xi_T, f_xi_h, g_xi_h = RO_tendency(par[i], T[i], h[i], xi_T[i], xi_h[i], EF[i])
 

        if NM == "EM":
            T[i + 1], noise_out[i, 0] = EM_scheme(T[i], f_T, g_T, dt, noise_in[i, 0])
            h[i + 1], noise_out[i, 1] = EM_scheme(h[i], f_h, g_h, dt, noise_in[i, 1])
            xi_T[i + 1], noise_out[i, 2] = EM_scheme(xi_T[i], f_xi_T, g_xi_T, dt, noise_in[i, 2])
            xi_h[i + 1], noise_out[i, 3] = EM_scheme(xi_h[i], f_xi_h, g_xi_h, dt, noise_in[i, 3])

        elif NM == "EH":
            T_s, noise_out[i, 0] = EM_scheme(T[i], f_T, g_T, dt, noise_in[i, 0])
            h_s, noise_out[i, 1] = EM_scheme(h[i], f_h, g_h, dt, noise_in[i, 1])
            xi_T_s, noise_out[i, 2] = EM_scheme(xi_T[i], f_xi_T, g_xi_T, dt, noise_in[i, 2])
            xi_h_s, noise_out[i, 3] = EM_scheme(xi_h[i], f_xi_h, g_xi_h, dt, noise_in[i, 3])

            f_T_s, g_T_s, f_h_s, g_h_s, f_xi_T_s, g_xi_T_s, f_xi_h_s, g_xi_h_s = RO_tendency(
                par[i + 1], T_s, h_s, xi_T_s, xi_h_s, EF[i + 1]
            )

            T[i + 1], _ = EM_scheme(T[i], 0.5 * (f_T + f_T_s), 0.5 * (g_T + g_T_s), dt, noise_out[i, 0])
            h[i + 1], _ = EM_scheme(h[i], 0.5 * (f_h + f_h_s), 0.5 * (g_h + g_h_s), dt, noise_out[i, 1])
            xi_T[i + 1], _ = EM_scheme(xi_T[i], 0.5 * (f_xi_T + f_xi_T_s), 0.5 * (g_xi_T + g_xi_T_s), dt, noise_out[i, 2])
            xi_h[i + 1], _ = EM_scheme(xi_h[i], 0.5 * (f_xi_h + f_xi_h_s), 0.5 * (g_xi_h + g_xi_h_s), dt, noise_out[i, 3])

    return T, h, noise_out
#####################################################
#####################################################
