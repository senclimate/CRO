import os
import sys
import numpy as np
import pandas as pd

from .fit_LR import fit_LR
from .fit_MLE import fit_MLE

def func_default_fitting_method(par_option_T, par_option_h, par_option_noise, table_path='table_default_fitting_method.txt'):
    """
    Determine the default fitting method for the CRO (Community Recharge Oscillator) 
    model based on prescribed parameter and noise options.

    The function reads a lookup table (`table_default_fitting_method.txt`) that maps 
    combinations of seasonality, linearity, noise color, and noise amplitude type 
    to a recommended fitting method.

    Parameters
    ----------
    par_option_T : list of int
        Options for prescribed terms in the SST equation.
    par_option_h : list of int
        Options for prescribed terms in the thermocline equation.
    par_option_noise : list of int
        Noise configuration array: [noise_color_T, noise_color_h, noise_amp_type].
    table_path : str, optional
        Path to the tab-delimited lookup table text file (default: 'table_default_fitting_method.txt').

    Returns
    -------
    str
        Recommended fitting method, e.g., 'LR-F', 'LR-C', 'LR-F-MAC', or 'MLE'. 
        Defaults to 'LR-F' if no matching entry is found.

    Notes
    -----
    - Determines whether the model is seasonal or constant.
    - Determines whether the dynamics are linear or nonlinear.
    - Determines the noise type (red/white) and amplitude (multiplicative/additive).
    - The function prints the selected fitting method.
    """
    
    # Get path relative to the script file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    table_path = os.path.join(script_dir, table_path)
    
    # Determine seasonality
    seasonal_type = "seasonal" if any(p > 1 for p in par_option_T + par_option_h) else "constant"

    # Determine linearity
    if par_option_T[2] == 0 and par_option_T[3] == 0 and par_option_T[4] == 0 and par_option_h[2] == 0:
        det_type = "linear"
    else:
        det_type = "nonlinear"

    # Noise settings
    if par_option_noise[0] == 0:
        noise_color_type = "red"
    elif par_option_noise[0] == 1:
        noise_color_type = "white"

    if par_option_noise[2] == 0:
        noise_amp_type = "multiplicative"
    elif par_option_noise[2] == 1:
        noise_amp_type = "multiplicative-H"
    elif par_option_noise[2] == 2:
        noise_amp_type = "additive"

    # Load the method lookup table from a tab-delimited text file
    df = pd.read_csv(table_path, sep='\s+')

    # Filter for matching row
    matched = df[
        (df["seasonal_type"] == seasonal_type) &
        (df["det_type"] == det_type) &
        (df["noise_color_type"].astype(str) == noise_color_type) &
        (df["noise_amp_type"].astype(str) == noise_amp_type)
    ]

    if not matched.empty:
        print("Referring to table_default_fitting_method.txt and using "+matched["fitting_method"].iloc[0])
        return matched["fitting_method"].iloc[0]
    else:
        print("Warning: No matching fitting method found. Defaulting to 'LR-F'.")
        return "LR-F"

########################################################################################################

def RO_fitting(T, h, par_option_T, par_option_h, par_option_noise, method_fitting=None, dt=None, verbose=True):
    """
    Fit the CRO (Community Recharge Oscillator) model parameters to T and h time series.

    This function performs parameter fitting for the RO system, using the specified
    prescribed parameter options and noise characteristics. Fitting can be done 
    via linear regression (LR) or maximum likelihood estimation (MLE), with optional
    handling for multiplicative or additive noise. It can automatically select a
    default fitting method based on parameter options.

    Parameters
    ----------
    T : ndarray
        Time series of SST anomalies (1D array).
    h : ndarray
        Time series of thermocline anomalies (1D array).
    par_option_T : dict
        Prescribed options for SST equation parameters.
    par_option_h : dict
        Prescribed options for thermocline equation parameters.
    par_option_noise : dict
        Noise settings, e.g., {'T': 'red', 'h': 'red', 'T_type': 'multiplicative'}.
    method_fitting : str, optional
        Fitting method to use: 'LR-F', 'LR-C', 'LR-F-MAC', or 'MLE'. 
        If None, the method is determined automatically.
    dt : float, optional
        Time step of the input time series. Default is 1.0 (months).
    verbose : bool, optional
        If True, prints fitting information and progress.

    Returns
    -------
    dict
        Dictionary of fitted CRO parameters:
        Keys include:
        'R', 'F1', 'F2', 'epsilon', 'b_T', 'c_T', 'd_T', 'b_h',
        'sigma_T', 'sigma_h', 'B', 'm_T', 'm_h', 'n_T', 'n_h', 'n_g'.

    Notes
    -----
    - Automatically removes mean from input time series.
    - Applies Ito-to-Stratonovich correction for multiplicative noise.
    - Negligible seasonal amplitudes are rounded to zero.
    - Ensures consistent noise settings for T and h equations.
    - Prints a summary of the fitting process if `verbose` is True.

    Examples
    --------
    >>> par_option_T = {'mean': 1, 'seasonal': 3, 'semi_annual': 0, ...}
    >>> par_option_h = {'mean': 1, 'seasonal': 3, 'semi_annual': 0, ...}
    >>> par_option_noise = {'T': 'red', 'h': 'red', 'T_type': 'multiplicative'}
    >>> fitted_params = RO_fitting(T, h, par_option_T, par_option_h, par_option_noise)
    >>> print(fitted_params['sigma_T'])
    """

    fitting_option_red = "ARn"   # "LR" or "AR1" or "ARn"

    ### Checking fitting setups ###
    if verbose:
        print("---------------------------------------------------------------------------------")
        print("Welcome to CRO Fitting! Your fitting setups:")
        print("---------------------------------------------------------------------------------")

    if dt is None:
        dt = 1.0
        if verbose:
            print(f" - Data time step is not given, defaulting to: dt = 1.0 months.")
    if verbose:
        print(f" - Time series length: N = len(T)*dt = {len(T)*dt} months.")

        print(f" - Prescribed terms: {par_option_T}. \n"
              f"                     {par_option_h}. \n"
              "   0 - Do not prescribe. \n"
              "   1 - Prescribe only the annual mean. \n"
              "   3 - Prescribe the annual mean and annual seasonality. \n"
              "   5 - Prescribe the annual mean, annual seasonality, and semi-annual seasonality.")

    if par_option_noise['T'] != par_option_noise['h']:
        raise ValueError(f"par_option_noise['T'] = {par_option_noise['T']}, "
                         f"par_option_noise['h'] = {par_option_noise['h']}\n"
                         "Fitting methods for T and h equations should be the same.")

    if verbose:
        print(f" - Noise options: {par_option_noise}.")
        if (par_option_noise['T'] == 'red') and (par_option_noise['h'] == 'red'):
            print(f" - Fitting method for T and h red noises: {fitting_option_red}.\n"
                  f"   This option is defined internally within fit.py.\n"
                  f"   Options available are: LR or AR1 or ARn.")

    # Convert parameter options
    par_option_T_vals = list(par_option_T.values())
    par_option_h_vals = list(par_option_h.values())
    noise_keys = ['T', 'h', 'T_type']

    noise_map = {
        'white': 1,
        'red': 0,
        'additive': 2,
        'multiplicative': 0,
        'multiplicative-H': 1
    }
    par_option_noise_array = [noise_map[str(par_option_noise[k])] for k in noise_keys]

    if verbose:
        print(f" - Fitting method for T and h main equations: {method_fitting}.")

    if method_fitting is None:
        method_fitting = func_default_fitting_method(par_option_T_vals, par_option_h_vals, par_option_noise_array)

    # Perform fitting
    if method_fitting == "LR-F":
        par = fit_LR(T, h, par_option_T_vals, par_option_h_vals,
                            par_option_noise_array, dt, "F", "LR", fitting_option_red)
    elif method_fitting == "LR-C":
        par = fit_LR(T, h, par_option_T_vals, par_option_h_vals,
                            par_option_noise_array, dt, "C", "LR", fitting_option_red)
    elif method_fitting == "LR-F-MAC":
        par = fit_LR(T, h, par_option_T_vals, par_option_h_vals,
                            par_option_noise_array, dt, "F", "MAC", fitting_option_red)
    elif method_fitting == "MLE":
        par = fit_MLE(T, h, par_option_T_vals, par_option_h_vals,
                             par_option_noise_array, dt)
    else:
        raise ValueError(f"Unknown fitting method: {method_fitting}")

    if method_fitting in ("LR-F", "LR-F-MAC", "MLE"):
        if par[15][0] == 0:   # Ito to Stratonovich Conversion for multiplicative noise
            par[0][0] = par[0][0] - 0.5 * (par[8][0] * par[10][0])**2
        elif par[15][0] == 1: # Ito to Stratonovich Conversion for Heaviside multiplicative noise
            par[0][0] = par[0][0] - 0.25 * (par[8][0] * par[10][0])**2
        elif par[15][0] == 2: # No correction needed for additive noise
            par[0][0] = par[0][0]

    ### Round negligible seasonal amplitudes of linear/nonlinear parameters to zero ###
    for i in range(0, 8):
        if par[i][1] < 1e-4:
            par[i][1] = 0.0
            par[i][2] = 0.0
        if par[i][3] < 1e-4:
            par[i][3] = 0.0
            par[i][4] = 0.0

    # Filter out zero values
    Arr = []
    for i in range(0, 13):
        arr = par[i]
        arr = arr[arr != 0]
        Arr.append(arr)

    # Organize into dictionary
    param_keys = [
        "R", "F1", "F2", "epsilon", "b_T", "c_T", "d_T", "b_h",
        "sigma_T", "sigma_h", "B", "m_T", "m_h", "n_T", "n_h", "n_g"
    ]
    param_values = [
        Arr[0], Arr[1], Arr[5], Arr[6], Arr[2], Arr[3], Arr[4], Arr[7],
        Arr[8], Arr[9], Arr[10], Arr[11], Arr[12],
        [int(par[13][0])], [int(par[14][0])], [int(par[15][0])]
    ]

    par = dict(zip(param_keys, param_values))
    par = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in par.items()}

    ### Final print ###
    if verbose:
        print("---------------------------------------------------------------------------------")
        print("All steps are successfully completed!")
        print("---------------------------------------------------------------------------------")

    return par
