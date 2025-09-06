import numpy as np
import statsmodels.api as sm

def wrap_to_pi(angle):
    return (angle + np.pi) % (2 * np.pi) - np.pi

def heaviside(x):
    return np.where(x > 0, 1, 0)

def func_mon_std(X, dt=1.0):
    X = np.array(X).flatten()
    months_per_year = int(12 / dt)
    num_years = int(len(X) / months_per_year)
    
    if len(X) % months_per_year != 0:
        X = X[:months_per_year * num_years]  # trim excess
    
    X_mon = X.reshape((num_years, months_per_year))
    X_mon_std = np.std(X_mon, axis=0)
    return X_mon_std


def RO_BWJ(par):
    # Extract annual mean values
    R_value = par['R'][0]
    F1_value = par['F1'][0]
    epsilon_value = par['epsilon'][0]
    F2_value = par['F2'][0]

    # Calculate growth rate and frequency
    gr = (R_value - epsilon_value) / 2
    w = np.sqrt(4 * F1_value * F2_value - (R_value + epsilon_value)**2) / 2

    # Return complex index
    return gr + 1j * w


def regress_std(y, X, alpha=0.05):
    X = np.asarray(X)
    y = np.asarray(y).flatten()


    # Detect intercept column
    intercept_idx = np.where(np.all(X == 1, axis=0))[0]
    has_intercept = len(intercept_idx) > 0


    # Remove intercept for standardization
    X_noint = np.delete(X, intercept_idx, axis=1) if has_intercept else X

    mu_X = X_noint.mean(axis=0)
    sigma_X = X_noint.std(axis=0, ddof=0)

    if np.isscalar(sigma_X):
        sigma_X = np.array([sigma_X])
    sigma_X[sigma_X == 0] = 1
    X_std = (X_noint - mu_X) / sigma_X

    mu_y = y.mean()
    sigma_y = y.std(ddof=0)
    if sigma_y == 0:
        sigma_y = 1
    y_std = (y - mu_y) / sigma_y

    # Add intercept back
    if has_intercept:
        X_reg = np.column_stack([np.ones(len(X)), X_std])
    else:
        X_reg = X_std

    # Regress
    model = sm.OLS(y_std, X_reg).fit()
    b_std = model.params
    r_std = model.resid

    b = np.zeros(b_std.shape)
    if has_intercept:
        non_idx = [i for i in range(X.shape[1]) if i not in intercept_idx]
        b[non_idx] = b_std[1:] * (sigma_y / sigma_X)
        b[intercept_idx] = mu_y - np.sum((mu_X / sigma_X) * b_std[1:]) * sigma_y
    else:
        b = b_std * (sigma_y / sigma_X)

    if len(b) == 1:
        X = X.reshape(-1, 1)  
        b = b.reshape(1, -1)  
        y_hat = X @ b  
    else:
        y_hat = X @ b
    r = y - y_hat
    return b, None, r