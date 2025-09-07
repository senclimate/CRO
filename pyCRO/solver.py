import sys
import numpy as np

#############################
### Euler-Maruyama scheme ###

def EM_scheme(x, f, g, dt, noise=None):
    """
    Euler-Maruyama scheme.

    Supports scalar or vectorized inputs:
    x, f, g can be scalars or numpy arrays of shape (NE,).

    If noise is None, generates standard normal random noise of same shape.
    If noise contains np.nan, replaces with fresh noise.

    Returns:
        xs : new x after one step
        noise : the actual noise used
    """
    if noise is None:
        noise = np.random.randn(*np.shape(x))
    else:
        # Replace any np.nan in noise with new random draws
        noise = np.where(np.isnan(noise), np.random.randn(*np.shape(x)), noise)

    dW = np.sqrt(dt) * noise
    xs = x + f * dt + g * dW
    return xs, noise


### This module generates time-evolving RO parameters and external forcing for use in the RO solver. ###
########################################################################################################
def par_processing(par_in, dt, NT): 
    num_params = len(par_in)
    par_out = np.ones((NT, num_params))

    ro_time = np.linspace(0, (NT - 1) * dt, NT)
   
    for i, par in enumerate(par_in):
        params = par_in[par]
        if np.isscalar(params): # if 'n_T'/'n_h'/'n_g' is given as a scalar instead of list, convert it to list
            params = [params]   
        
        if (par == 'n_T') or (par == 'n_h') or (par == 'n_g'):
            par_out[:, i] = params[0]
        else:
            w = 2 * np.pi / 12
            par_out[:, i] = (params[0] + params[1] * np.sin(w * ro_time + params[2]) +
                             params[3] * np.sin(2 * w * ro_time + params[4]))
    return par_out # array of shape (NT, 16) or (NT, 2) containing time-evolving parameters
########################################################################################################

### RO tendency equations ###
#############################
def RO_tendency(par, T, h, xi_T, xi_h, EF):
    """
    Vectorized RO tendency function.

    Supports scalar or array inputs (NE,).
    """
    # Unpack parameters
    R, F1, F2, epsilon = par[0:4]
    b_T, c_T, d_T, b_h = par[4:8]
    sigma_T, sigma_h, B = par[8:11]
    m_T, m_h = par[11:13]
    n_T, n_h, n_g = map(int, par[13:16])  # Ensure integer indexing

    # Unpack external forcing
    E_T, E_h = EF

    ### Noise factor ###
    ####################
    if n_g == 0:
        g_T_factor = B * T
    elif n_g == 1:
        g_T_factor = B * np.maximum(T, 0)  # Heaviside(T)
    elif n_g == 2:
        g_T_factor = 0.0
    else:
        raise ValueError("Invalid value for n_g")
    ####################

    ### dT/dt ###
    #############
    if n_T == 0:  # red noise
        f_T = R*T + F1*h + b_T*(T**2) - c_T*(T**3) + d_T*T*h + sigma_T*(1 + g_T_factor)*xi_T + E_T
        g_T = 0.0
    elif n_T == 1:  # white noise
        f_T = R*T + F1*h + b_T*(T**2) - c_T*(T**3) + d_T*T*h + E_T
        g_T = sigma_T * (1 + g_T_factor)
    else:
        raise ValueError("Invalid value for n_T")

    ### dh/dt ###
    ##############
    if n_h == 0:  # red noise
        f_h = -F2*T - epsilon*h - b_h*(T**2) + sigma_h*xi_h + E_h
        g_h = 0.0
    elif n_h == 1:  # white noise
        f_h = -F2*T - epsilon*h - b_h*(T**2) + E_h
        g_h = sigma_h
    else:
        raise ValueError("Invalid value for n_h")

    ### d(xi_T)/dt ###
    ##################
    if n_T == 0:
        f_xi_T = -m_T * xi_T
        g_xi_T = np.sqrt(2 * m_T)
    else:
        f_xi_T = np.nan
        g_xi_T = np.nan

    ### d(xi_h)/dt ###
    ##################
    if n_h == 0:
        f_xi_h = -m_h * xi_h
        g_xi_h = np.sqrt(2 * m_h)
    else:
        f_xi_h = np.nan
        g_xi_h = np.nan

    return f_T, g_T, f_h, g_h, f_xi_T, g_xi_T, f_xi_h, g_xi_h
#############################
#############################

### RO integrater (Euler-Maruyama and Euler-Huen) ###
#####################################################
def RO_integral(par, EF, NM, NT, dt, T0, h0, noise_all):
    """
    Vectorized RO integrator over ensemble members.

    Inputs:
        par      : (NT, 16) array
        EF       : (NT, 2) array
        T0, h0   : initial conditions (NE,)
        noise_all: (NT-1, 4, NE)
    """
    NE = T0.size

    T    = np.zeros((NT, NE))
    h    = np.zeros((NT, NE))
    xi_T = np.zeros((NT, NE))
    xi_h = np.zeros((NT, NE))

    T[0,:] = T0
    h[0,:] = h0

    noise_out = np.full((NT-1, 4, NE), np.nan)

    for i in range(NT - 1):
        f_T, g_T, f_h, g_h, f_xi_T, g_xi_T, f_xi_h, g_xi_h = RO_tendency(
            par[i], T[i], h[i], xi_T[i], xi_h[i], EF[i]
        )

        if NM == "EM":
            T[i+1], noise_out[i,0] = EM_scheme(T[i], f_T, g_T, dt, noise_all[i,0])
            h[i+1], noise_out[i,1] = EM_scheme(h[i], f_h, g_h, dt, noise_all[i,1])
            xi_T[i+1], noise_out[i,2] = EM_scheme(xi_T[i], f_xi_T, g_xi_T, dt, noise_all[i,2])
            xi_h[i+1], noise_out[i,3] = EM_scheme(xi_h[i], f_xi_h, g_xi_h, dt, noise_all[i,3])

        elif NM == "EH":
            T_s, noise_out[i,0] = EM_scheme(T[i], f_T, g_T, dt, noise_all[i,0])
            h_s, noise_out[i,1] = EM_scheme(h[i], f_h, g_h, dt, noise_all[i,1])
            xi_T_s, noise_out[i,2] = EM_scheme(xi_T[i], f_xi_T, g_xi_T, dt, noise_all[i,2])
            xi_h_s, noise_out[i,3] = EM_scheme(xi_h[i], f_xi_h, g_xi_h, dt, noise_all[i,3])

            f_T_s, g_T_s, f_h_s, g_h_s, f_xi_T_s, g_xi_T_s, f_xi_h_s, g_xi_h_s = RO_tendency(
                par[i+1], T_s, h_s, xi_T_s, xi_h_s, EF[i+1]
            )

            T[i+1], _ = EM_scheme(T[i], 0.5*(f_T+f_T_s), 0.5*(g_T+g_T_s), dt, noise_out[i,0])
            h[i+1], _ = EM_scheme(h[i], 0.5*(f_h+f_h_s), 0.5*(g_h+g_h_s), dt, noise_out[i,1])
            xi_T[i+1], _ = EM_scheme(xi_T[i], 0.5*(f_xi_T+f_xi_T_s), 0.5*(g_xi_T+g_xi_T_s), dt, noise_out[i,2])
            xi_h[i+1], _ = EM_scheme(xi_h[i], 0.5*(f_xi_h+f_xi_h_s), 0.5*(g_xi_h+g_xi_h_s), dt, noise_out[i,3])

    return T, h, noise_out
#####################################################
#####################################################


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
def RO_solver(par, IC, N, NE, NM="EH", dt=0.1, saveat=1.0, savemethod="sampling", EF=None, noise_custom=None, verbose=True):
    """
    Numerical solution of the Recharge Oscillator (RO) model with stochastic forcing.

    This function integrates the Recharge Oscillator (RO) system numerically using
    either the Euler–Maruyama (EM) or Euler–Heun (EH) scheme. It includes
    deterministic dynamics, stochastic noise (white or red), and optional
    external forcing. Ensemble simulations are supported, with results returned
    at user-specified output intervals.

    Parameters
    ----------
    par : dict
        Dictionary of model parameters with the following keys (as 1-element arrays):

        - ``'R'`` : float
            Damping parameter.
        - ``'F1'`` : float
            Feedback parameter relating thermocline depth to SST.
        - ``'epsilon'`` : float
            Thermocline damping parameter.
        - ``'F2'`` : float
            Feedback parameter relating SST to thermocline.
        - ``'sigma_T'`` : float
            Noise amplitude for SST.
        - ``'sigma_h'`` : float
            Noise amplitude for thermocline depth.
        - ``'B'`` : float
            Multiplicative noise coefficient (used when ``n_g`` = 0 or 1).
        - ``'n_T'`` : int
            Noise type in SST equation (``1`` = white noise, ``0`` = red noise).
        - ``'m_T'`` : array-like
            Memory kernel for red noise in SST equation (ignored if ``n_T=1``).
        - ``'n_h'`` : int
            Noise type in thermocline equation (``1`` = white noise, ``0`` = red noise).
        - ``'m_h'`` : array-like
            Memory kernel for red noise in thermocline equation (ignored if ``n_h=1``).
        - ``'n_g'`` : int
            Noise structure in SST equation:
            ``0`` = multiplicative, ``1`` = multiplicative + Heaviside, ``2`` = additive.
    IC : tuple of float
        Initial condition ``(T0, h0)``:
        
        - ``T0`` : initial SST anomaly.
        - ``h0`` : initial thermocline depth anomaly.
    N : float
        Total simulation length (time units).
    NE : int
        Number of ensemble members.
    NM : {'EM', 'EH'}, optional
        Numerical integration method. Default is ``'EH'`` (Euler–Heun).
    dt : float, optional
        Numerical integration time step. Default is ``0.1``.
    saveat : float, optional
        Output saving interval. Must be divisible by ``dt``. Default is ``1.0``.
    savemethod : {'sampling', 'mean'}, optional
        Method for saving results:
        
        - ``'sampling'`` : store values every ``saveat`` steps (default).
        - ``'mean'`` : average values within each ``saveat`` interval.
    EF : dict, optional
        External forcing with the following keys (as 1D arrays of length 5):
        
        - ``'E_T'`` : SST forcing coefficients.
        - ``'E_h'`` : Thermocline forcing coefficients.
        
        If ``None`` (default), no external forcing is applied.
    noise_custom : {None, int, ndarray}, optional
        Specification of stochastic noise:
        
        - ``None`` (default): generate new Gaussian noise for each ensemble.
        - ``int`` : random seed for reproducible noise (same across ensembles).
        - ``ndarray`` of shape (NT-1, 4, NE) : user-provided noise realizations.
    verbose : bool, optional
        If ``True`` (default), print detailed setup and progress messages.

    Returns
    -------
    T_out : ndarray of shape (N_out, NE)
        SST anomalies from the numerical integration, after applying the saving scheme.
    h_out : ndarray of shape (N_out, NE)
        Thermocline anomalies from the numerical integration, after applying the saving scheme.
    noise_out : ndarray of shape (NT-1, 4, NE)
        Realizations of Gaussian noise used in the simulation.

    Notes
    -----
    - White noise forcing corresponds to ``n_T = 1`` or ``n_h = 1``.
    - Red noise forcing (``n_T = 0`` or ``n_h = 0``) requires a nonzero memory kernel ``m_T`` or ``m_h``.
    - For the Euler–Maruyama method (``NM = 'EM'``) with multiplicative noise,
      an Ito-to-Stratonovich correction is applied automatically.

    Examples
    --------
    >>> par = {
    ...     'R':[0.5], 'F1':[1.2], 'epsilon':[0.3], 'F2':[0.8],
    ...     'sigma_T':[0.2], 'sigma_h':[0.1], 'B':[0.05],
    ...     'n_T':[1], 'm_T':[0], 'n_h':[1], 'm_h':[0], 'n_g':[0]
    ... }
    >>> IC = (0.1, -0.2)
    >>> T, h, noise = RO_solver(par, IC, N=50, NE=5, dt=0.1, saveat=1.0)
    >>> T.shape, h.shape
    ((51, 5), (51, 5))
    """
    
    NT = int(round(N / dt))
    step = int(round(saveat / dt))
    ratio = saveat / dt

    ### Checking simulation setups ###
    if verbose:
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
        raise ValueError(f"dt = {dt}, saveat = {saveat}: The numerical time step (dt) must be less than or equal to the data saving interval (saveat).")

    ### Check Initial conditions ###
    if np.array(IC).shape != (2,):
        raise ValueError(f"Invalid 'IC' input: expected an array with shape (2,), but got {IC}.")
    elif verbose:
        print(f" - Initial conditions: IC = [T0, h0] = {IC}")

    ### Check input parameter shapes ###
    if not has_same_shape(par, par_shapes):
        raise ValueError("Input parameters do not have the expected shape.")
    elif verbose:
        print(f" - Input parameters have the expected shapes.")

    ### Check noise parameters ###
    if verbose:
        if np.array(par['n_T']).item() == 1:
            print(" - 'n_T' = 1: White noise forcing in T; 'm_T' ignored.")
        elif np.array(par['n_T']).item() == 0:
            print(f" - 'n_T' = 0: Red noise forcing in T with m_T = {par['m_T']}.")
            if np.all(np.array(par['m_T']) == 0):
                raise ValueError(f"'m_T' needs to have at least one non-zero argument for red noise forcing.")
        else:
            raise ValueError(f"n_T = {np.array(par['n_T']).item()} is an invalid option: it has to be either '1' (red noise) or '0' (white noise).")

        if np.array(par['n_h']).item() == 1:
            print(" - 'n_h' = 1: White noise forcing in h; 'm_h' ignored.")
        elif np.array(par['n_h']).item() == 0:
            print(f" - 'n_h' = 0: Red noise forcing in h with m_h = {par['m_h']}.")
            if np.all(np.array(par['m_h']) == 0):
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

    ### Check numerical integration method ###
    if verbose:
        if NM == "EM":
            print(f" - Numerical integration method: NM = '{NM}' (Euler–Maruyama method)")
            print("   The default method is NM = 'EH' (Euler–Heun method)")
        elif NM == "EH":
            print(f" - Numerical integration method: NM = '{NM}' (Euler–Heun method; default)")
        else:
            raise ValueError(
            f"{NM} is an invalid NM input. Expected 'EM' (Euler–Maruyama method) or "
            f"'EH' (Euler–Heun method)")

    ### Check data saving method ###
    if verbose:
        if savemethod == "sampling":
            print(f" - Data saving method: savemethod = {savemethod} (default)")
        elif savemethod == "mean":
            print(f" - Data saving method: savemethod = {savemethod} (default is 'sampling')")
        else:
            raise ValueError(f"{savemethod} is an invalid savemethod: it has to be either 'sampling' or 'mean'.")

    ### Check external forcing ###
    if EF is None:
        if verbose:
            print(" - External forcing is not given, therefore using\n   EF = {'E_T': [0.0, 0.0, 0.0, 0.0, 0.0], 'E_h': [0.0, 0.0, 0.0, 0.0, 0.0]}.")
        EF = {'E_T': [0.0, 0.0, 0.0, 0.0, 0.0], 'E_h': [0.0, 0.0, 0.0, 0.0, 0.0]}
    elif has_same_shape(EF, EF_shapes):
        if verbose:
            print(f" - External forcing is given as\n  EF = {EF}.")
    else:
        raise ValueError(f"- External forcing does not match the shape of {EF_shapes}.")

    ### Check noise ###
    if verbose:
        if noise_custom is None:
            print(f" - noise_custom = {noise_custom}: System-generated noise is used and changes at every run.")
        elif np.isscalar(noise_custom):
            print(f" - noise_custom = {noise_custom}: seeded same noise.")
        elif noise_custom.shape == (NT-1,4,NE):
            print(f" - You provided the customized noises with shapes (N/dt-1, 4) = ({NT-1}, 4).")
        else:
            raise ValueError(f"Invalid 'noise_custom' input: expected 'None', a scalar (int or float), "
                         f"or an array with shape (N/dt-1, 4) = ({NT-1}, 4).")
        print("---------------------------------------------------------------------------------")

    ### preprocessing parameters ###
    if (NM == "EM") and (par['n_T'][0] == 1): # Do Ito to Strantonovich Conversion
                                              # when Euler-Maruya scheme is used with white noise forcing
        if par['n_g'][0] == 0:   # Ito to Strantonovich Conversion for multiplicative noise
            par['R'][0] = par['R'][0] + 0.5 * (par['sigma_T'][0] * par['B'][0])**2
        elif par['n_g'][0] == 1: # Ito to Strantonovich Conversion for Heaviside multiplicative noise
            par['R'][0] = par['R'][0] + 0.25 * (par['sigma_T'][0] * par['B'][0])**2
        elif par['n_g'][0] == 2: # No correction needed for additive noise
            par['R'][0] = par['R'][0]

    par_in = par_processing(par, dt, NT)
    EF_in = par_processing(EF, dt, NT)

    ### initial conditions & noise ###
    T0 = np.full(NE, IC[0])
    h0 = np.full(NE, IC[1])

    if noise_custom is None:
        noise_all = np.random.randn(NT-1, 4, NE)
    elif np.isscalar(noise_custom):
        np.random.seed(noise_custom)
        noise_all = np.random.randn(NT-1, 4, NE)
    else:
        noise_all = noise_custom

    ### integration ###
    if verbose:
        print(f"Numerical integration starting:")
        print("---------------------------------------------------------------------------------")
    T_out, h_out, noise_out = RO_integral(par_in, EF_in, NM, NT, dt, T0, h0, noise_all)
    if verbose:
        print("---------------------------------------------------------------------------------")

    ### save ###
    if (savemethod == "sampling") or (savemethod is None):
        T_out = T_out[::step, :]
        h_out = h_out[::step, :]
    elif savemethod == "mean":
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

    if verbose:
        print("All steps successfully completed!")
        print("---------------------------------------------------------------------------------")

    return T_out, h_out, noise_out

#################
#################
