import sys
import numpy as np

def RO_analytic_std(par):
    """
    Compute analytical standard deviation of T and h for the RO model.

    Parameters
    ----------
    par : dict
        Dictionary containing parameter arrays:
        'R', 'F1', 'epsilon', 'F2', 'sigma_T', 'sigma_h'.

    Returns
    -------
    T_std : float
        Standard deviation of T.
    h_std : float
        Standard deviation of h.
    """
    
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


def RO_analytic_solver(par, IC, N, NE, dt=0.1, saveat=1.0, savemethod="sampling", noise_custom=None):
    """
    Analytical solution of the Recharge Oscillator (RO) model with deterministic and stochastic forcing.

    This function computes the analytical solution of the Recharge Oscillator (RO)
    model for a given parameter set, initial condition, and noise realization(s).
    Both deterministic and stochastic contributions are included. 
    The stochastic integrals are evaluated using vectorized convolution like summations across all ensemble members.

    Parameters
    ----------
    par : dict
        Dictionary of model parameters, each as a one-element array:
        
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
    IC : tuple of float
        Initial condition ``(To, ho)``:
        
        - ``To`` : initial SST anomaly.
        - ``ho`` : initial thermocline depth anomaly.
    N : float
        Total simulation length (time units).
    NE : int
        Number of ensemble members.
    dt : float, optional
        Numerical integration time step. Default is ``0.1``.
    saveat : float, optional
        Output saving interval. Must be divisible by ``dt``. Default is ``1.0``.
    savemethod : {'sampling', 'mean'}, optional
        Method for saving results:
        
        - ``'sampling'`` : take samples every ``saveat``.
        - ``'mean'`` : average values over each ``saveat`` interval.
    noise_custom : {None, int, ndarray}, optional
        Specification of stochastic noise:
        
        - ``None`` (default): new Gaussian noise is generated for each ensemble.
        - ``int`` : seed for reproducible noise (same across ensembles).
        - ``ndarray`` of shape ``(NT-1, 4)`` : user-provided noise, repeated across ensembles.
        - ``ndarray`` of shape ``(NT-1, 4, NE)`` : user-provided noise for each ensemble.

    Returns
    -------
    T_anal_out : ndarray of shape (N_out, NE)
        SST anomalies including deterministic and stochastic forcing, 
        after applying the saving scheme. The final row is NaN-padded 
        for dimensional consistency.
    h_anal_out : ndarray of shape (N_out, NE)
        Thermocline anomalies including deterministic and stochastic forcing, 
        after applying the saving scheme. The final row is NaN-padded.
    noise_array : ndarray of shape (NT-1, 4, NE)
        Realizations of Gaussian noise used in the simulation.

    Notes
    -----
    The analytical solution consists of:
    
    - Deterministic part: exponential-sinusoidal functions of time.
    - Stochastic part: weighted sums (discrete convolutions) of noise 
      with exponential and sinusoidal kernels.
    
    Ensemble members are evaluated simultaneously using vectorized 
    NumPy broadcasting, avoiding explicit Python loops.

    Examples
    --------
    >>> par = {'R':[0.5], 'F1':[1.2], 'epsilon':[0.3], 'F2':[0.8],
    ...        'sigma_T':[0.2], 'sigma_h':[0.1]}
    >>> IC = (0.1, -0.2)
    >>> T, h, noise = RO_solver_analytic(par, IC, N=50, NE=10, dt=0.1, saveat=1.0)
    >>> T.shape, h.shape
    ((51, 10), (51, 10))
    """
    
    NT = int(round(N / dt))
    step = int(round(saveat / dt))
    ratio = saveat / dt

    ### Checking simulation setups ###
    ##################################
    print("---------------------------------------------------------------------------------")
    print("Welcome to the CRO Analytical Solver!")
    print("Ensure that the same arguments are provided as for RO_solver,\n except that 'NM' and 'EF' are not required.")
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

    ########################################################
    # Noise generation (vectorized for all NE)
    if noise_custom is None:
        noise_array = np.random.randn(NT - 1, 4, NE)
    elif np.isscalar(noise_custom):
        np.random.seed(int(noise_custom))
        noise_array = np.random.randn(NT - 1, 4, NE)
    elif noise_custom.shape == (NT - 1, 4):
        noise_array = np.repeat(noise_custom[..., None], NE, axis=2)
    elif noise_custom.shape == (NT - 1, 4, NE):
        noise_array = noise_custom
    else:
        raise ValueError(f"Invalid 'noise_custom' input: expected None, scalar, shape ({NT-1},4) or ({NT-1},4,{NE}).")

    w_T = sigma_T_value * noise_array[:, 0, :] / np.sqrt(dt)  # (NT-1, NE)
    w_h = sigma_h_value * noise_array[:, 1, :] / np.sqrt(dt)  # (NT-1, NE)

    # Vectorized convolution-like integral for stochastic terms
    T_sto = np.zeros((NT - 1, NE))
    h_sto = np.zeros((NT - 1, NE))

    ### analytical calculations
    ########################################################
    for i in range(1, NT - 1):
        tau = t[:i+1]                     # (i+1,)
        delta = t[i] - tau                 # (i+1,)
        sin = np.sin(w * delta)            # (i+1,)
        cos = np.cos(w * delta)            # (i+1,)
        expgr = np.exp(gr * delta)         # (i+1,)

        # (i+1, NE) via broadcasting
        T_kernel = (
            w_T[:i+1, :] * cos[:, None] +
            w_T[:i+1, :] * (R_value + epsilon_value) / (2 * w) * sin[:, None] +
            w_h[:i+1, :] * (F1_value / w) * sin[:, None]
        )
        h_kernel = (
            w_h[:i+1, :] * cos[:, None] -
            w_h[:i+1, :] * (R_value + epsilon_value) / (2 * w) * sin[:, None] -
            w_T[:i+1, :] * (F2_value / w) * sin[:, None]
        )

        T_sto[i, :] = np.sum(dt * expgr[:, None] * T_kernel, axis=0)
        h_sto[i, :] = np.sum(dt * expgr[:, None] * h_kernel, axis=0)

    # Add deterministic and stochastic contributions
    T_anal_out = np.vstack([T_det[:-1, None] + T_sto, np.full((1, NE), np.nan)])
    h_anal_out = np.vstack([h_det[:-1, None] + h_sto, np.full((1, NE), np.nan)])

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
    return T_anal_out, h_anal_out, noise_array
