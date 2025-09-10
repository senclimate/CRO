import sys
import numpy as np

from scipy.optimize import fsolve
from scipy.optimize import least_squares
from scipy.interpolate import interp1d
from scipy.linalg import solve, inv
import warnings

from .utils import wrap_to_pi

########################################################################################################
def fit_MLE_white(T, h, par_option_T, par_option_h, par_option_noise, dt, raw_output=False):
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



def fit_MLE_red(T, h, par_option_T, par_option_h, par_option_noise, dt):

    par_option_T = np.array(par_option_T)
    par_option_h = np.array(par_option_h)
    par_option_noise = np.array(par_option_noise)


    n_T, n_h, n_g = par_option_noise.astype(int)


    # Setting Option
    par_option_T_onoff = (par_option_T >= 1).astype(int)
    par_option_h_onoff = (par_option_h >= 1).astype(int)
    par_option_T_h_onoff = np.concatenate([par_option_T_onoff, par_option_h_onoff])

    par_option_T_season = np.where(par_option_T > 1)[0]
    par_option_h_season = np.where(par_option_h > 1)[0]
    par_option_T_h_season = np.where(np.concatenate([par_option_T, par_option_h]) > 1)[0]


    par_option_T_h = np.concatenate([par_option_T, par_option_h])


    # Initial Guess for White Noise Fitting Value
    par_white_in = fit_MLE_white(T, h, par_option_T, par_option_h, [1, 1, n_g], dt)


    # only select annual mean value
    R, F1, b_T, c_T, d_T, F2, epsilon, b_h = par_white_in[0:8, 0]
    sigma_T, sigma_h, sigma_B =  par_white_in[8:11, 0]
    sigma_B = sigma_B*sigma_T
    m_T = 1.0
    m_h = 1.0
    
    # Setting Seasonal Parameters
    par_season_T = par_white_in[0:5, 1:5]   
    par_season_h = par_white_in[5:8, 1:5]



    sigma_T2 = 0.01
    sigma_h2 = 0.01

    # Interpolation of T and h
    dt_interp = 0.01 # empirically dt_interp should be 0.01 or smaller 
    t = np.arange(len(T)) * dt
    if dt > dt_interp:
        time_interp = np.arange(0, t[-1] + dt_interp, dt_interp)
        T = interp1d(t, T, kind='linear')(time_interp)
        h = interp1d(t, h, kind='linear')(time_interp)
        dt = dt_interp
        t = time_interp


    w = 2 * np.pi / 12
    t_shift = t + 0.5*dt
    n = len(T)

    # Setting x diff 
    x = np.stack([np.diff(T), np.diff(h)]) 
    x = np.expand_dims(x, axis=1)          


    if n_g == 0:   # Linear (B*T)
        noise_g = np.ones_like(T) 
    elif n_g == 1: # Heavisdie (B*H(T)*T)
        noise_g = np.heaviside(T, 0) 
    elif n_g == 2: # Additive (B=0)
        noise_g = np.zeros_like(T) 

    # --- Iterative estimation ---
    print("-------------------------------------------------------------------")
    print("Hang tight — red noise MLE fitting can take a bit!")
    print("-------------------------------------------------------------------")
    niter = 10 #100
    for _ in range(niter): 
        print(f"{_ + 1} out of {niter} iterations")
        # Setting CGNS Matrix
        A0_T = R*T + F1*h + b_T*T**2 - c_T*T**3 + d_T*T*h
        A0_h = -F2*T - epsilon*h - b_h*T**2


        Tv = np.array([T, h, T**2, -T**3, T*h]).T
        hv = np.array([-T, -h, -T**2]).T



        A0_T_add = ( par_season_T[0,0] * Tv[:,0] * np.sin(w*t_shift) + par_season_T[0,1] * Tv[:,0] * np.cos(w*t_shift)
                   + par_season_T[0,2] * Tv[:,0]*np.sin(2*w*t_shift) + par_season_T[0,3] * Tv[:,0]*np.cos(2*w*t_shift)
                   + par_season_T[1,0] * Tv[:,1] * np.sin(w*t_shift) + par_season_T[1,1] * Tv[:,1] * np.cos(w*t_shift)
                   + par_season_T[1,2] * Tv[:,1]*np.sin(2*w*t_shift) + par_season_T[1,3] * Tv[:,1]*np.cos(2*w*t_shift)
                   + par_season_T[2,0] * Tv[:,2] * np.sin(w*t_shift) + par_season_T[2,1] * Tv[:,2] * np.cos(w*t_shift)
                   + par_season_T[2,2] * Tv[:,2]*np.sin(2*w*t_shift) + par_season_T[2,3] * Tv[:,2]*np.cos(2*w*t_shift)
                   + par_season_T[3,0] * Tv[:,3] * np.sin(w*t_shift) + par_season_T[3,1] * Tv[:,3] * np.cos(w*t_shift)
                   + par_season_T[3,2] * Tv[:,3]*np.sin(2*w*t_shift) + par_season_T[3,3] * Tv[:,3]*np.cos(2*w*t_shift)
                   + par_season_T[4,0] * Tv[:,4] * np.sin(w*t_shift) + par_season_T[4,1] * Tv[:,4] * np.cos(w*t_shift)
                   + par_season_T[4,2] * Tv[:,4]*np.sin(2*w*t_shift) + par_season_T[4,3] * Tv[:,4]*np.cos(2*w*t_shift))
        A0_T = A0_T + A0_T_add

        A0_h_add = ( par_season_h[0,0] * hv[:,0] * np.sin(w*t_shift) + par_season_h[0,1] * hv[:,0] * np.cos(w*t_shift)
                   + par_season_h[0,2] * hv[:,0]*np.sin(2*w*t_shift) + par_season_h[0,3] * hv[:,0]*np.cos(2*w*t_shift)
                   + par_season_h[1,0] * hv[:,1] * np.sin(w*t_shift) + par_season_h[1,1] * hv[:,1] * np.cos(w*t_shift)
                   + par_season_h[1,2] * hv[:,1]*np.sin(2*w*t_shift) + par_season_h[1,3] * hv[:,1]*np.cos(2*w*t_shift)
                   + par_season_h[2,0] * hv[:,2] * np.sin(w*t_shift) + par_season_h[2,1] * hv[:,2] * np.cos(w*t_shift)
                   + par_season_h[2,2] * hv[:,2]*np.sin(2*w*t_shift) + par_season_h[2,3] * hv[:,2]*np.cos(2*w*t_shift))
        A0_h = A0_h + A0_h_add


        A0 = np.stack([A0_T, A0_h])

        A1_T = np.column_stack([
               sigma_T + sigma_B * noise_g * T,
               np.zeros_like(T)])

        A1_h = np.column_stack([
               np.zeros_like(T),
               sigma_h * np.ones_like(T)])

        A1 = np.stack([A1_T, A1_h], axis=2)
        A1 = np.transpose(A1, (2, 1, 0)) 

        # Define noise model
        B1 = np.diag([sigma_T2, sigma_h2])
        a0 = np.zeros((2, 1))
        a1 = np.diag([-m_T, -m_h])
        b2 = np.identity(2)



        mu_f = np.zeros((2, n))
        R_f = np.zeros((2, 2, n))
        R_f[:, :, 0] = 0.01 * np.eye(2)


        # Forward filter
        for i in range(n - 1):
            A0_sel = A0[:, i]
            A0_sel = A0_sel[:, np.newaxis]             

            A1_sel = A1[:, :, i]          

            R_f_sel = R_f[:, :, i] 

            mu_f_sel = mu_f[:, i]          
            mu_f_sel = mu_f_sel[:, np.newaxis]

            x_sel = x[:, :, i] 

            dummy1 = a0 + a1 @ mu_f[:, [i]]
            dummy2 = R_f_sel @ A1_sel.T
            dummy3 = np.linalg.inv(B1 @ B1.T)
            dummy4 = x_sel / dt - (A0_sel + A1_sel @ mu_f_sel)
            mu_f_tendency = dummy1 + dummy2 * dummy3 @ dummy4

            R_f_tendency = (a1 @ R_f_sel
                          + R_f_sel @ a1.T
                          + b2 @ b2.T
                          - (R_f_sel @ A1_sel.T) @ np.linalg.inv(B1 @ B1.T) @ (A1_sel @ R_f_sel))

            mu_f[:, i+1] = mu_f[:, i] + (mu_f_tendency[:, 0] * dt)
            R_f[:, :, i+1] = R_f[:, :, i] + (R_f_tendency * dt)

        # Backward smoother
        mu_s = np.zeros((2, n))
        R_s = np.zeros((2, 2, n))
  

        mu_s[:, -1] = mu_f[:, -1]
        R_s[:, :, -1] = R_f[:, :, -1]

 
        for i in range(n - 2, -1, -1):
            #print(i)
            R_f_sel = R_f[:, :, i+1]
            mu_f_sel = mu_f[:, i+1]
            mu_f_sel = mu_f_sel[:, np.newaxis]
            R_s_sel = R_s[:, :, i+1]
            mu_s_sel = mu_s[:, i+1]
            mu_s_sel = mu_s_sel[:, np.newaxis]


            mu_s_tendency = -a0 - a1 @ mu_s_sel + (b2 @ b2.T) @ np.linalg.inv(R_f_sel) @ (mu_f_sel - mu_s_sel)
            R_s_tendency = -(a1 + (b2 @ b2.T) @ np.linalg.inv(R_f_sel)) @ R_s_sel \
                           - R_s_sel @ (a1.T + (b2 @ b2.T) @ R_f_sel) \
                           + (b2 @ b2.T)


            mu_s[:, i] = mu_s[:, i+1] + (mu_s_tendency * dt).flatten()
            R_s[:, :, i] = R_s[:, :, i+1] + R_s_tendency * dt
 

        # Calculate C matrix and y moments
        C = np.zeros_like(R_f)
        I = np.eye(2)
        A1_dt = I + a1 * dt
        B2 = b2 @ b2.T

        for i in range(R_f.shape[2]):
            R_f_sel = R_f[:, :, i]
            denom = B2 * dt + A1_dt @ R_f_sel @ A1_dt.T
            C[:, :, i] = R_f_sel @ A1_dt.T @ np.linalg.inv(denom)


        yi = mu_s.copy()
        yiyi = np.zeros_like(R_s)
        for i in range(R_s.shape[2]):
            R_s_sel = R_s[:, :, i]
            mu_s_sel = mu_s[:, i].reshape(-1, 1)  
            yiyi[:, :, i] = R_s_sel + mu_s_sel @ mu_s_sel.T


        yi1yi = np.zeros((2, 2, R_s.shape[2] - 1))
        for i in range(R_s.shape[2] - 1):
            C_i = C[:, :, i]
            mu_next = mu_s[:, i + 1].reshape(-1, 1)  
            mu_now = mu_s[:, i].reshape(-1, 1)       
            yi1yi[:, :, i] = R_s[:, :, i + 1] @ C_i.T + mu_next @ mu_now.T


        # Setting M and Rinv 
        Rinv_single = np.diag([1 / sigma_T2**2, 1 / sigma_h2**2, 1, 1]) * (1 / dt)
        Rinv = np.repeat(Rinv_single[:, :, np.newaxis], n, axis=2)  

        xi_T = yi[0, :].reshape(-1, 1)  
        xi_h = yi[1, :].reshape(-1, 1)  

        T_col = T.reshape(-1, 1)
        h_col = h.reshape(-1, 1)
        xi_T_col = xi_T.reshape(-1, 1)
        xi_h_col = xi_h.reshape(-1, 1)
        noise_g_col = noise_g.reshape(-1, 1)

        MT_raw = np.hstack([T_col, h_col, T_col**2, -T_col**3, T_col * h_col, xi_T_col, xi_T_col * noise_g_col * T_col, np.zeros((T.shape[0], 6))])
        MT_weights = np.concatenate([par_option_T_onoff, [1, 1], np.zeros(6)])
        MT = MT_raw * dt * MT_weights

        Mh_raw = np.hstack([np.zeros((T.shape[0], 7)), -T_col, -h_col, -T_col**2, xi_h_col, np.zeros((T.shape[0], 2))])
        Mh_weights = np.concatenate([np.zeros(7), par_option_h_onoff, [1], np.zeros(2)])
        Mh = Mh_raw * dt * Mh_weights

        MxiT_raw = np.hstack([np.zeros((T.shape[0], 11)), -xi_T_col, np.zeros((T.shape[0], 1))])
        MxiT_weights = np.concatenate([np.zeros(11), [1], [0]])
        MxiT = MxiT_raw * dt * MxiT_weights

        Mxih_raw = np.hstack([np.zeros((T.shape[0], 12)),-xi_h_col])
        Mxih_weights = np.concatenate([np.zeros(12), [1]])
        Mxih = Mxih_raw * dt * Mxih_weights


        if len(par_option_T_season) > 0:
            for i in range(len(par_option_T_season)):
                idx = par_option_T_season[i]  # column index to modulate

                MT_col = MT[:, idx]  # column to modulate, shape (n,)
        
                if par_option_T[idx] == 3:  # annual only
                    MT_add = np.column_stack([MT_col * np.sin(w * t_shift), MT_col * np.cos(w * t_shift)])
                elif par_option_T[idx] == 5:  # annual + semiannual
                    MT_add = np.column_stack([
                        MT_col * np.sin(w * t_shift),
                        MT_col * np.cos(w * t_shift),
                        MT_col * np.sin(2 * w * t_shift),
                        MT_col * np.cos(2 * w * t_shift)])
                else:
                    continue  # skip unsupported codes
        
                # Stack MT_add and corresponding zero-padding to other matrices
                MT = np.hstack([MT, MT_add])
                zeros_like_MT_add = np.zeros_like(MT_add)
                Mh = np.hstack([Mh, zeros_like_MT_add])
                MxiT = np.hstack([MxiT, zeros_like_MT_add])
                Mxih = np.hstack([Mxih, zeros_like_MT_add])

        if len(par_option_h_season) > 0:
            for i in range(len(par_option_h_season)):
                idx = par_option_h_season[i]  # Index within par_option_h (0-based)
                mh_col_idx = 7 + idx          # Adjusted index in Mh to find the correct feature column

                Mh_col = Mh[:, mh_col_idx]    # shape (n,)

                if par_option_h[idx] == 3:  # annual only
                    Mh_add = np.column_stack([Mh_col * np.sin(w * t_shift), Mh_col * np.cos(w * t_shift)])
                elif par_option_h[idx] == 5:  # annual + semiannual
                    Mh_add = np.column_stack([Mh_col * np.sin(w * t_shift),
                                          Mh_col * np.cos(w * t_shift),
                                          Mh_col * np.sin(2 * w * t_shift),
                                          Mh_col * np.cos(2 * w * t_shift)])
                else:
                    continue  # skip unsupported cases

                # Stack Mh_add and corresponding zero-padding to other matrices
                Mh = np.hstack([Mh, Mh_add])
                zeros_like_Mh_add = np.zeros_like(Mh_add)
                MT = np.hstack([MT, zeros_like_Mh_add])
                MxiT = np.hstack([MxiT, zeros_like_Mh_add])
                Mxih = np.hstack([Mxih, zeros_like_Mh_add])

        M = np.stack([MT, Mh, MxiT, Mxih], axis=2)  
        M = np.transpose(M, (2, 1, 0))  
        Mtr = np.transpose(M, (1, 0, 2))  

        if n_g == 0 or n_g == 1:  # multiplicative
            par_option_matrix = np.array([
                          np.concatenate([par_option_T_onoff, [1], [1], np.zeros(6)]),           # MT
                          np.concatenate([np.zeros(7), par_option_h_onoff, [1], np.zeros(2)]),   # Mh
                          np.concatenate([np.zeros(11), [1], [0]]),                               # MxiT
                          np.concatenate([np.zeros(11), [0], [1]])                                # Mxih
                          ])
        elif n_g == 2:  # additive
            par_option_matrix = np.array([
                          np.concatenate([par_option_T_onoff, [1], [0], np.zeros(6)]),
                          np.concatenate([np.zeros(7), par_option_h_onoff, [1], np.zeros(2)]),
                          np.concatenate([np.zeros(11), [1], [0]]),
                          np.concatenate([np.zeros(11), [0], [1]])
                          ])



        total_T = np.sum(np.where(par_option_T > 1, par_option_T - 1, 0))
        total_h = np.sum(np.where(par_option_h > 1, par_option_h - 1, 0))

        par_option_matrix = np.array([
        np.concatenate([par_option_matrix[0], np.ones(total_T), np.zeros(total_h)]),
        np.concatenate([par_option_matrix[1], np.zeros(total_T), np.ones(total_h)]),
        np.concatenate([par_option_matrix[2], np.zeros(total_T + total_h)]),
        np.concatenate([par_option_matrix[3], np.zeros(total_T + total_h)])
        ])

        zero_col_inds = np.where(np.sum(par_option_matrix, axis=0) == 0)[0]
        M = np.delete(M, zero_col_inds, axis=1)     
        Mtr = np.delete(Mtr, zero_col_inds, axis=0)

        # Calculate thetao_A
        temp = np.einsum('ijk,jlk->ilk', Mtr, Rinv)
        thetao_A = np.einsum('ijk,jlk->ilk', temp, M)

        nM = thetao_A.shape[0] - total_T - total_h
        nT = np.sum(par_option_T_onoff == 1)


        thetao_A[nT, nT, :] = dt * yiyi[0, 0, :] / (sigma_T2**2)
        thetao_A[nM - 3, nM - 3, :] = dt * yiyi[1, 1, :] / (sigma_h2**2)
        thetao_A[nM - 2, nM - 2, :] = dt * yiyi[0, 0, :] / 1
        thetao_A[nM - 1, nM - 1, :] = dt * yiyi[1, 1, :] / 1

        if n_g == 0 or n_g == 1:
            thetao_A[nT + 1, nT + 1, :] = dt * (noise_g * T**2) * yiyi[0, 0, :] / (sigma_T2**2)

        thetao_A = thetao_A[:, :, :-1]

        # Calculate thetao_B
        x_NaN = np.zeros_like(x)  
        x_all = np.concatenate([x, x_NaN], axis=0)  

        Mtr_sub = Mtr[:, :, :-1]          
        Rinv_sub = Rinv[:, :, :-1]        
        x_sub = x_all[:, :, :]            

        temp = np.einsum('ijk,jlk->ilk', Mtr_sub, Rinv_sub)  

        thetao_B = np.einsum('ijk,jlk->ilk', temp, x_sub)    
        thetao_B[nM - 2, 0, :] = -yi1yi[0, 0, :] + yiyi[0, 0, :-1]
        thetao_B[nM - 1, 0, :] = -yi1yi[1, 1, :] + yiyi[1, 1, :-1]


        # Sum across the time dimension (axis=2), from ind_res to end-ind_res+1
        thetao_A_sum = np.sum(thetao_A[:, :, :], axis=2)
        thetao_B_sum = np.sum(thetao_B[:, :, :], axis=2)
        thetao_B_sum = thetao_B_sum.squeeze()

        thetao = np.linalg.solve(thetao_A_sum, thetao_B_sum)

        M_sliced = M[:, :, :-1]             
        thetao_squeezed = np.squeeze(thetao)  
        M_thetao = np.matmul(M_sliced.transpose(2, 0, 1), thetao_squeezed).transpose(1, 0)  
        x_sliced = x_all[:, :, :] if x_all.ndim == 3 else x_all[:, :]
        R_A = x_sliced - M_thetao[:, np.newaxis, :]  
        R_B = np.transpose(R_A, (1, 0, 2))  
        R_AB = np.einsum('ijk,jlk->ilk', R_A, R_B) 


        adjust1 = ((sigma_T * xi_T[:-1, 0]) + sigma_B * noise_g[:-1] * xi_T[:-1, 0] * T[:-1]) * dt
        adjust1 = adjust1 ** 2
        R_AB[0, 0, :] = R_AB[0, 0, :] - np.squeeze(adjust1)

        adjust2 = ((sigma_h * xi_h[:-1, 0]) * dt) 
        adjust2 = adjust2 ** 2
        R_AB[1, 1, :] = R_AB[1, 1, :] - np.squeeze(adjust2)

        R_AB[0, 0, :] += (
        np.squeeze(yiyi[0, 0, :-1]) *
        (dt * (sigma_T + sigma_B * noise_g[:-1] * T[:-1])) ** 2)

        R_AB[1, 1, :] += (
        np.squeeze(yiyi[1, 1, :-1]) *
        (dt * sigma_h) ** 2)


        # Compute the thetao term (same index used in both the square and linear terms)
        theta_idx_1 = len(thetao) - 2 - total_T - total_h
        theta_term_1 = thetao[theta_idx_1] * dt - 1
        theta_idx_2 = len(thetao) - 1 - total_T - total_h
        theta_term_2 = thetao[theta_idx_2] * dt - 1

        R_AB[2, 2, :] = (
        np.squeeze(yiyi[0, 0, 1:]) +
        np.squeeze(yiyi[0, 0, :-1]) * (theta_term_1 ** 2) +
        2 * np.squeeze(yi1yi[0, 0, :]) * theta_term_1)

        R_AB[3, 3, :] = (
        np.squeeze(yiyi[1, 1, 1:]) +
        np.squeeze(yiyi[1, 1, :-1]) * (theta_term_2 ** 2) +
        2 * np.squeeze(yi1yi[1, 1, :]) * theta_term_2)

        # Sum over time (axis=2), then divide by the total number of time steps in R_AB
        R_M = np.sum(R_AB[:, :, : ], axis=2) / R_AB.shape[2]


        # Save to Parameter
        thetao_out = np.zeros((13 + total_T + total_h, 1))

        selected_indices = np.where(np.sum(par_option_matrix, axis=0) == 1)[0]

        thetao_out[selected_indices, 0] = thetao.flatten()

        R =       thetao_out[0, 0]  
        F1 =      thetao_out[1, 0]  
        b_T =     thetao_out[2, 0]  
        c_T =     thetao_out[3, 0]  
        d_T =     thetao_out[4, 0]  
        sigma_T = thetao_out[5, 0]  
        sigma_B = thetao_out[6, 0]
        B =       sigma_B/sigma_T   
        F2 =      thetao_out[7, 0]  
        epsilon = thetao_out[8, 0]  
        b_h =     thetao_out[9, 0]  
        sigma_h = thetao_out[10, 0] 
        m_T =     thetao_out[11, 0] 
        m_h =     thetao_out[12, 0] 
        thetao_out = list(thetao_out.flatten())

        par_T_h = []
        par_order_array=[0, 1, 2, 3, 4, 7, 8, 9]
        j = 0
        for i in range(0, len(par_option_T_h)):
            seasonal_list = []
            if (par_option_T_h[i] == 0) or (par_option_T_h[i] == 1):
                seasonal_list.extend([0, 0, 0, 0])
            elif(par_option_T_h[i] == 3):
                seasonal_list = thetao_out[8+j:8+j+2]
                seasonal_list.extend([0, 0])
                j += 2
            elif(par_option_T_h[i] == 5):
                seasonal_list = thetao_out[8+j:8+j+4]
                j += 4
            sin_a, cos_a, sin_sa, cos_sa = seasonal_list
            A_a = np.sqrt(sin_a**2 + cos_a**2)
            phi_a = np.mod(np.arctan2(cos_a, sin_a), 2 * np.pi)
            A_sa = np.sqrt(sin_sa**2 + cos_sa**2)
            phi_sa = np.mod(np.arctan2(cos_sa, sin_sa), 2 * np.pi)

            par_T_h.append([thetao_out[par_order_array[i]], A_a, phi_a, A_sa, phi_sa])
        par_T_h = np.array(par_T_h)
        sigma_T2 = 0.01
        sigma_h2 = 0.01

    par_noise = np.array([[sigma_T / np.sqrt(2 * m_T), 0, 0, 0, 0],
                 [sigma_h / np.sqrt(2 * m_h), 0, 0, 0, 0],
                 [sigma_B/sigma_T, 0, 0, 0, 0],
                 [m_T, 0, 0, 0, 0],
                 [m_h, 0, 0, 0, 0],
                 [n_T, 0, 0, 0, 0],
                 [n_h, 0, 0, 0, 0],
                 [n_g, 0, 0, 0, 0]])
    par_combined = np.vstack([par_T_h, par_noise])

    return par_combined



def fit_MLE(T, h, par_option_T, par_option_h, par_option_noise, dt):
    """
    Estimate recharge oscillator (RO) parameters using **Maximum Likelihood Estimation (MLE)**

    Parameters
    ----------
    T : array_like
        Time series of SST anomalies (1D).
    h : array_like
        Time series of thermocline depth anomalies (1D).
    par_option_T : array_like
        Switches for which SST-related parameters to include  
        (0=off, 1=on, 3=annual, 5=annual+semiannual).
    par_option_h : array_like
        Switches for which h-related parameters to include (same coding).
    par_option_noise : array_like of length 3
        Noise structure options `[n_T, n_h, n_g]`:  
        - n_T : noise option for T-equation  
        - n_h : noise option for h-equation  
        - n_g : type of multiplicative noise term  
            * 0 = linear multiplicative (`B*T`)  
            * 1 = Heaviside multiplicative (`B*H(T)*T`)  
            * 2 = additive (`B=0`)
    dt : float
        Time step (months).

    Returns
    -------
    par : ndarray
        Final estimated parameter set including deterministic RO parameters,
        seasonal amplitudes/phases, and red noise variances.

    Notes
    -----
    - Implements iterative forward filtering and backward smoothing 
      to estimate hidden noise states.  
    - This function may take a long time to converge (default 10 iterations, 
      could be increased for accuracy).  
    """
    
    n_T = par_option_noise[0]
    n_h = par_option_noise[1]

    if n_T == 1 and n_h == 1:
        par = fit_MLE_white(T, h, par_option_T, par_option_h, par_option_noise, dt)
    elif n_T == 0 and n_h == 0:
        par = fit_MLE_red(T, h, par_option_T, par_option_h, par_option_noise, dt)
    else:
        raise ValueError("Error: Mixed noise types not supported by MLE")
    
    return par