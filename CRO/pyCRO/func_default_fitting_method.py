import pandas as pd
import sys

def func_default_fitting_method(par_option_T, par_option_h, par_option_noise, table_path='./Data/table_default_fitting_method.txt'):
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
    df = pd.read_csv(table_path, delim_whitespace=True)

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




