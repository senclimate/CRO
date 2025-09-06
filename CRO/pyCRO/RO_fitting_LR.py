import numpy as np
from RO_fitting_LR import *
from RO_fitting_MAC import *
from regress_std import *
import sys

def wrap_to_pi(angle):
    return (angle + np.pi) % (2 * np.pi) - np.pi

def heaviside(x):
    return np.where(x > 0, 1, 0)

def RO_fitting_LR(T, h, par_option_T, par_option_h, par_option_noise, dt, tend_option, fitting_option_B, fitting_option_red):
    T = np.ravel(T)
    h = np.ravel(h)

    par_option_T = np.array(par_option_T)
    par_option_h = np.array(par_option_h)
    par_option_noise = np.array(par_option_noise)

    n_T, n_h, n_g = par_option_noise.astype(int)

    par_option_T_onoff = (par_option_T >= 1).astype(int)
    par_option_h_onoff = (par_option_h >= 1).astype(int)


    # Derivative
    if tend_option == "F":
        Tdiff = np.diff(T) / dt
        hdiff = np.diff(h) / dt
        Ts = T[:-1]
        hs = h[:-1]
    elif tend_option == "C":
        Tdiff = (T[2:] - T[:-2]) / (2 * dt)
        hdiff = (h[2:] - h[:-2]) / (2 * dt)
        Ts = T[1:-1]
        hs = h[1:-1]

    # Design matrices
    Tv = np.column_stack([Ts, hs, Ts**2, -Ts**3, Ts * hs]) * par_option_T_onoff
    hv = np.column_stack([-Ts, -hs, -Ts**2]) * par_option_h_onoff


    # Seasonal features
    w = 2 * np.pi / 12
    t = np.arange(0, len(T)) * dt

    if tend_option == "F":
        t = t[:-1]
        t_shift = t + 0.5 * dt
    elif tend_option == "C": 
        t = t[1:-1]
        t_shift = t 



    def add_seasonal(X=Tv, par_option=par_option_T, option="T"):
        X_seasonal_extended = []
        for i in range(0, len(par_option)):
            if par_option[i] == 0:
                on_and_off = np.array([0, 0, 0, 0, 0])
            elif par_option[i] == 1:
                on_and_off = np.array([1, 0, 0, 0, 0])
            elif par_option[i] == 3:
                on_and_off = np.array([1, 1, 1, 0, 0])
            elif par_option[i] == 5:
                on_and_off = np.array([1, 1, 1, 1, 1])
            else:
                raise ValueError(f"Invalid {option} fitting option {par_option}: "
                                  "Value must be one of the following — 0 (constant zero), "
                                  "1 (constant non-zero), 3 (annual variation), or 5 (annual + semi-annual variation).")


       
            X_seasonal = np.column_stack([
                        X[:, i] * np.sin(w * t_shift),
                        X[:, i] * np.cos(w * t_shift),
                        X[:, i] * np.sin(2 * w * t_shift),
                        X[:, i] * np.cos(2 * w * t_shift)])
            X_seasonal = np.hstack([X[:, i].reshape(-1, 1), X_seasonal])
            X_seasonal = X_seasonal * on_and_off
            X_seasonal_extended.append(X_seasonal)

        X_seasonal_extended = np.array(X_seasonal_extended)
        X_seasonal_extended = np.transpose(X_seasonal_extended, (1, 0, 2))  # shape becomes (1199, 5, 5)
        if option == "T":
        	X_seasonal_extended = X_seasonal_extended.reshape(X.shape[0], 25) 
        if option == "h":
            X_seasonal_extended = X_seasonal_extended.reshape(X.shape[0], 15)
        return X_seasonal_extended

    Tv = add_seasonal(Tv, par_option_T, option="T")
    hv = add_seasonal(hv, par_option_h, option="h")

    par_T, _, res_T = regress_std(Tdiff, Tv)
    par_h, _, res_h = regress_std(hdiff, hv)


    def extract_seasonal(params, option):
        my_params = []
        if option == "T":
            index = 5
        elif option == "h":
            index = 3
        
        params = params.reshape(index, 5)
        for i in range(0, index):
            A_a = np.sqrt(params[i,1]**2 + params[i,2]**2)
            phi_a = np.mod(np.arctan2(params[i,2], params[i,1]), 2*np.pi)
 
            A_sa = np.sqrt(params[i,3]**2 + params[i,4]**2)
            phi_sa = np.mod(np.arctan2(params[i,4], params[i,3]), 2*np.pi)

            if np.pi < phi_a:
                phi_a = phi_a - 2*np.pi
            if np.pi < phi_sa:
                phi_sa = phi_sa - 2*np.pi
            my_param = [params[i,0], A_a, phi_a, A_sa, phi_sa]
            my_params.append(my_param)
        my_params = np.array(my_params) 

        return my_params


    par_T = extract_seasonal(par_T, option="T")
    par_h = extract_seasonal(par_h, option="h")


    # Noise fitting
    res_T_norm = res_T * np.sqrt(dt)
    res_h_norm = res_h * np.sqrt(dt)

    if n_T == 1 and n_h == 1 and n_g == 2:    # (white, white, additive)
        sigma_T = np.std(res_T_norm)
        sigma_h = np.std(res_h_norm)
        B = m_T = m_h = 0.0 # np.nan 

    elif n_T == 0 and n_h == 0 and n_g == 2:  # (red, red, additive)
        if fitting_option_red == "LR":
            res_Tdiff = np.diff(res_T) / dt
            res_Ts = res_T[:-1]
            m_T, _, _ = regress_std(res_Tdiff, res_Ts)

            res_hdiff = np.diff(res_h) / dt
            res_hs = res_h[:-1]
            m_h, _, _ = regress_std(res_hdiff, res_hs)
            m_T = np.abs(m_T.item())
            m_h = np.abs(m_h.item())
        elif fitting_option_red == "AR1":
            m_T = -np.log(np.corrcoef(res_T[:-1], res_T[1:])[0, 1]) / dt
            m_h = -np.log(np.corrcoef(res_h[:-1], res_h[1:])[0, 1]) / dt
        elif fitting_option_red == "ARn":
            # m_T part
            r, lags = np.correlate(res_T - np.mean(res_T), res_T - np.mean(res_T), mode='full'), np.arange(-len(res_T)+1, len(res_T))
            r = r / np.max(r)
            lags_dt = lags * dt
            positive_lags = lags_dt >= 0
            r = r[positive_lags]
            lags = lags_dt[positive_lags]

            try:
                cutoff = np.where(r < 0)[0][0]
                r = r[:cutoff]
                lags = lags[:cutoff]
            except IndexError:
                pass

            y = -np.log(r)
            X = lags.reshape(-1, 1)
            m_T, _, _ = regress_std(y, X)

            # m_h part
            r, lags = np.correlate(res_h - np.mean(res_h), res_h - np.mean(res_h), mode='full'), np.arange(-len(res_h)+1, len(res_h))
            r = r / np.max(r)
            lags_dt = lags * dt
            positive_lags = lags_dt >= 0
            r = r[positive_lags]
            lags = lags_dt[positive_lags]

            try:
                cutoff = np.where(r < 0)[0][0]
                r = r[:cutoff]
                lags = lags[:cutoff]
            except IndexError:
                pass

            y = -np.log(r)
            X = lags.reshape(-1, 1)
            m_h, _, _ = regress_std(y, X)
            m_T = m_T.item() 
            m_h = m_h.item()
        else:
            print(f"Invalid fitting method for red noise: {fitting_option_red}")
            sys.exit()



        sigma_T = np.std(res_T)
        sigma_h = np.std(res_h)
        B = 0 #np.nan

    elif n_T == 1 and n_h == 1 and n_g in (0, 1): # (white, white, multiplicative (either full or Heaviside))
        n = 10 * 12
        res_T_var = np.zeros(10**5)
        T_var = np.zeros(10**5)

        for i in range(10**5):
            ind = np.random.randint(0, len(res_T_norm), size=n)
            res_T_var[i] = np.var(res_T_norm[ind], ddof=1)  # ddof=1 for unbiased (same as MATLAB)
    
            if n_g == 0.0:
                T_var[i] = np.var(Ts[ind], ddof=1)
            elif n_g == 1.0:
                T_var[i] = np.var(np.heaviside(Ts[ind], 0) * Ts[ind], ddof=1)



        X = np.column_stack([T_var, np.ones(len(T_var))])
        par_T_res, _, _ = regress_std(res_T_var, X)
        if par_T_res[0] < 0:
            B = 0.0
        else:
            B = np.sqrt(par_T_res[0] / par_T_res[1])
        sigma_T = np.sqrt(par_T_res[1])
        sigma_h = np.std(res_h_norm)
        m_T = m_h = 0.0 

    elif n_T == 0 and n_h == 0 and n_g in (0, 1): # (red, red, multiplicative (either full or Heaviside))
        n = 10 * 12
        res_T_var, T_var = [], []
        for _ in range(10000):
            ind = np.random.randint(0, len(res_T_norm), size=n)
            res_T_var.append(np.var(res_T_norm[ind]))
            if n_g == 0:
                T_var.append(np.var(Ts[ind]))
            else:
                T_var.append(np.var(heaviside(Ts[ind]) * Ts[ind]))
        X = np.column_stack([T_var, np.ones(len(T_var))])
        res_T_var = np.array(res_T_var) 
        par_T_res, _, _ = regress_std(res_T_var, X)
        sigma_T = np.sqrt(par_T_res[1])

        if par_T_res[0] < 0:
            B = 0.0
        else:
            B = np.sqrt(par_T_res[0] / par_T_res[1])
        xi = res_T / (sigma_T * (1 + B * Ts))
        m_T, _, _ = regress_std(np.diff(xi), xi[:-1].reshape(-1, 1))
        m_T = np.abs(m_T.item())
        res_hdiff = np.diff(res_h) / dt
        res_hs = res_h[:-1]
        m_h, _, _ = regress_std(res_hdiff, res_hs.reshape(-1, 1))
        m_h = np.abs(m_h.item())
        sigma_h = np.std(res_h)


    #if (not np.isfinite(m_T) or m_T < 0.) or (not np.isfinite(m_h) or m_h < 0.):
    #    raise ValueError(f"Estimated m_T = {m_T} or m_h = {m_h} is not a valid number.\n"
    #                     f"This error typically occurs when trying to fit red noise expressions into the time series "
    #                     f"generated with white noise.\nTry using the fitting option for white noise instead.")

    #if np.isnan(B):
    #        raise ValueError(f"B = np.sqrt({par_T_res[0]} / {par_T_res[1]}) is NaN.\n"
    #                         f"This error typically occurs when trying to fit too many annual and semi-annual components "
    #                         f"at the same time with multiplicative noise.\n"
    #                         f"It can also occur when attempting to fit red noise expressions to a time series generated with white noise.\n"
    #                         f"Try reducing the seasonality components and/or using the fitting option for white noise instead.")
        #sys.exit(


    par_noise = np.array([[sigma_T, 0, 0, 0, 0], 
                 [sigma_h, 0, 0, 0, 0], 
                 [B, 0, 0, 0, 0], 
                 [m_T, 0, 0, 0, 0], 
                 [m_h, 0, 0, 0, 0],
                 [n_T, 0, 0, 0, 0],
                 [n_h, 0, 0, 0, 0],
                 [n_g, 0, 0, 0, 0]])

    final_par = np.vstack([par_T, par_h, par_noise])


    if (n_g == 0 or n_g == 1) and fitting_option_B == "MAC":
        final_par = RO_fitting_MAC(T, h, final_par)  # Calculate sigma_T and B using MAC, overwrite par

    return final_par
