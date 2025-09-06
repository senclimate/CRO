import numpy as np

def EM_scheme(x, f, g, dt, noise=None):
    if noise is None or np.isnan(noise):
        noise = np.random.randn()
    dW = np.sqrt(dt) * noise
    xs = x + f * dt + g * dW
    return xs, noise