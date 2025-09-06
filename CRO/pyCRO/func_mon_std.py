import numpy as np

def func_mon_std(X, dt=1.0):
    X = np.array(X).flatten()
    months_per_year = int(12 / dt)
    num_years = int(len(X) / months_per_year)
    
    if len(X) % months_per_year != 0:
        X = X[:months_per_year * num_years]  # trim excess
    
    X_mon = X.reshape((num_years, months_per_year))
    X_mon_std = np.std(X_mon, axis=0)
    return X_mon_std
