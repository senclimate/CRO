import numpy as np
import sys
from scipy.optimize import fsolve
from scipy.optimize import least_squares
import warnings

def wrap_to_pi(angle):
    return (angle + np.pi) % (2 * np.pi) - np.pi

def RO_fitting_MLE_white(T, h, par_option_T, par_option_h, par_option_noise, dt, raw_output=False):
    T = np.array(T).flatten()
    h = np.array(h).flatten()
    n = len(T)

    par_option_T = np.array(par_option_T)
    par_option_h = np.array(par_option_h)
    par_option_noise = np.array(par_option_noise)
    n_T, n_h, n_g = par_option_noise

    par_option_T_onoff = (par_option_T >= 1).astype(int)
    par_option_h_onoff = (par_option_h >= 1).astype(int)
    par_option_T_h_onoff = np.concatenate([par_option_T_onoff, par_option_h_onoff])

    par_option_T_season = np.where(par_option_T > 1)[0]
    par_option_h_season = np.where(par_option_h > 1)[0]
    par_option_T_h_season = np.where(np.concatenate([par_option_T, par_option_h]) > 1)[0]
    par_option_T_h = np.concatenate([par_option_T, par_option_h])

    # Construct MT and Mh
    MT = np.column_stack([T, h, T**2, -T**3, T * h, np.zeros((n, 3))]) * dt * np.concatenate([par_option_T_onoff, np.zeros(3)])
    Mh = np.column_stack([np.zeros((n, 5)), -T, -h, -T**2]) * dt * np.concatenate([np.zeros(5), par_option_h_onoff])


    w = 2 * np.pi / 12
    t = np.arange(n) * dt
    t_shift = t + 0.5*dt


    for idx in par_option_T_season:
        base = MT[:, idx]
        if par_option_T[idx] == 3:
            MT_add = np.column_stack([base * np.sin(w * t_shift), base * np.cos(w * t_shift)])
        elif par_option_T[idx] == 5:
            MT_add = np.column_stack([base * np.sin(w * t_shift), base * np.cos(w * t_shift),
                                      base * np.sin(2 * w * t_shift), base * np.cos(2 * w * t_shift)])
        else:
            continue
        MT = np.hstack([MT, MT_add])
        Mh = np.hstack([Mh, np.zeros_like(MT_add)])

    for idx in par_option_h_season:
        base = Mh[:, 5 + idx]
        if par_option_h[idx] == 3:
            Mh_add = np.column_stack([base * np.sin(w * t_shift), base * np.cos(w * t_shift)])
        elif par_option_h[idx] == 5:
            Mh_add = np.column_stack([base * np.sin(w * t_shift), base * np.cos(w * t_shift),
                                      base * np.sin(2 * w * t_shift), base * np.cos(2 * w * t_shift)])
        else:
            continue
        Mh = np.hstack([Mh, Mh_add])
        MT = np.hstack([MT, np.zeros_like(Mh_add)])



    # Assume MT and Mh are shape (1200, 8)
    M = np.stack([MT.T, Mh.T], axis=0)  
    M = M[:, :, :-1]                    
    Mtr = np.transpose(M, (1, 0, 2))    


    on_T = np.concatenate([par_option_T_onoff, np.zeros(3)])
    on_h = np.concatenate([np.zeros(5), par_option_h_onoff])


    ind = np.where((on_T == 0) & (on_h == 0))[0]

    M = np.delete(M, ind, axis=1)
    Mtr = np.delete(Mtr, ind, axis=0)



    x = np.stack([np.diff(T), np.diff(h)], axis=0) 
    x = x[:, np.newaxis, :]                        

    sigma_T = 1.0
    sigma_h = 1.0
    sigma_B = 0.0

    if n_g == 0:
        noise_g = np.ones_like(T) 
    elif n_g == 1:
        noise_g = np.heaviside(T, 0) 
    elif n_g == 2:
        noise_g = np.zeros_like(T) 


    if n_g in [0, 1]:  # multiplicative
        n_iter = 100
    elif n_g == 2:     # additive
        n_iter = 1


    for my_iter in range(n_iter):
        # Setting R Matrix 
        RT = np.column_stack([((sigma_T + sigma_B * noise_g * T)**2), np.zeros(n)]) * dt  
        Rh = np.column_stack([np.zeros(n), np.ones(n) * sigma_h**2]) * dt 
        R = np.stack([RT, Rh], axis=2)         
        R = np.transpose(R, (2, 1, 0))        

        with np.errstate(divide='ignore', invalid='ignore'):
            Rinv = 1 / R
        Rinv[np.isinf(Rinv)] = 0

        R = R[:, :, :-1]
        Rinv = Rinv[:, :, :-1]

        thetao_A = np.einsum('ikt,klt,ljt->ij', Mtr, Rinv, M)
        thetao_B = np.einsum('ikt,klt,ljt->ij', Mtr, Rinv, x)
        thetao = np.linalg.solve(thetao_A, thetao_B)


        # Calculate R Matrix
        R_A = x - np.einsum('ijk,jl->ilk', M, thetao)
        R_B = np.transpose(R_A, (1, 0, 2))
        # Calculate Noise for Additive Noise RO
        R = np.sum(R_A * R_B, axis=2) / (n - 1)

        # Calculate Noise for Multiplicative Noise RO
        if my_iter == 0:
            sigma_T = np.sqrt(R[0, 0] / dt)  # % Initial Guess for sigma_T using Additive Noise
            sigma_B = 0.0                    # % Initial Guess for sigma_B

        

        def func(sigma_x):
            #print(sigma_x[0], sigma_x[1])
            r0 = (R_A[0, 0, :] * R_B[0, 0, :])
            r1 = -dt * (sigma_x[0] + sigma_x[1] * noise_g[:-1] * T[:-1])**2
            r2 = (sigma_x[0] + sigma_x[1] * noise_g[:-1] * T[:-1])**(-3) * T[:-1]
            r3 = (sigma_x[0] + sigma_x[1] * noise_g[:-1] * T[:-1])**(-3)
            r_val1 = np.sum((r0 + r1) * r2) / len(T)
            r_val2 = np.sum((r0 + r1) * r3) / len(T)
            val1 = r_val1
            val2 = r_val2

            return np.array([val1, val2])

        res = least_squares(func, x0=[sigma_T, sigma_B], bounds=([1e-6, -1], [np.inf, np.inf]))
        sigma_T, sigma_B = res.x
        sigma_h = np.sqrt(R[1, 1] / dt)




    # Post-Processing Output
    on_indices = np.where(par_option_T_h_onoff == 1)[0]
    thetao_ann = thetao[:len(on_indices)]
    thetao_season = thetao[len(on_indices):]

    par_T_h = [[0.0] * 5 for _ in range(8)]
    for idx, param_idx in enumerate(on_indices):
        par_T_h[param_idx][0] = float(thetao_ann[idx])


    if len(par_option_T_h_season) > 0:
        ind_base = 0
        for i in range(len(par_option_T_h_season)):
            idx = par_option_T_h_season[i]
            if par_option_T_h[idx] == 3:
                ind_sin = ind_base
                ind_cos = ind_base + 1

                A = np.sqrt(thetao_season[ind_sin]**2 + thetao_season[ind_cos]**2).item()
                phi = np.mod(np.arctan2(thetao_season[ind_cos], thetao_season[ind_sin]), 2 * np.pi).item()

                if np.pi < phi:
                    phi = phi - 2*np.pi

                par_T_h[idx][1:3] = [A, phi]

                ind_base = ind_cos + 1

            elif par_option_T_h[idx] == 5:
                ind_sin_a = ind_base
                ind_cos_a = ind_base + 1
                ind_sin_sa = ind_base + 2
                ind_cos_sa = ind_base + 3

                A_a = np.sqrt(thetao_season[ind_sin_a]**2 + thetao_season[ind_cos_a]**2).item()
                phi_a = np.mod(np.arctan2(thetao_season[ind_cos_a], thetao_season[ind_sin_a]), 2 * np.pi).item()
                A_sa = np.sqrt(thetao_season[ind_sin_sa]**2 + thetao_season[ind_cos_sa]**2).item()
                phi_sa = np.mod(np.arctan2(thetao_season[ind_cos_sa], thetao_season[ind_sin_sa]), 2 * np.pi).item()

                if np.pi < phi_a:
                    phi_a = phi_a - 2*np.pi
                if np.pi < phi_sa:
                    phi_sa = phi_sa - 2*np.pi

                par_T_h[idx][1:5] = [A_a, phi_a, A_sa, phi_sa]

                ind_base = ind_cos_sa + 1


    par_T_h = np.array(par_T_h)
    par_noise = np.array([[sigma_T, 0, 0, 0, 0],
                 [sigma_h, 0, 0, 0, 0],
                 [sigma_B/sigma_T, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0],
                 [n_T, 0, 0, 0, 0],
                 [n_h, 0, 0, 0, 0],
                 [n_g, 0, 0, 0, 0]])
    par_combined = np.vstack([par_T_h, par_noise])
    return par_combined


