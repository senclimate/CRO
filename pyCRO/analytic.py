import sys
import numpy as np

def RO_BWJ(par):
    """
    Compute the Recharge Oscillator Bjerknes–Wyrtki–Jin (BWJ) indices
    (linear ENSO growth rate and oscillation frequency).

    This function evaluates the linear stability of the recharge oscillator
    system, returning a complex eigenvalue composed of growth rate and
    oscillation frequency.

    Mathematical formulation
    ------------------------

    The BWJ system is defined as:

    .. math::

        BJ = \\frac{R - \\epsilon}{2}

        WF = \\frac{1}{2} \\sqrt{4 F_1 F_2 - (R + \\epsilon)^2}

    The complex eigenvalue is:

    .. math::

        \\lambda = BJ + i\\,WF

    where:
    
        - BJ: growth/decay rate
        
        - WF: oscillation frequency

    Parameters
    ----------
    par : dict
        Dictionary of model parameters. Must include:

        R : ndarray
            Recharge/discharge feedback parameter.
        F1 : ndarray
            Coupling coefficient for zonal wind–SST feedback.
        F2 : ndarray
            Coupling coefficient for thermocline feedback.
        epsilon : ndarray
            Damping (linear dissipation) term.

        Only the first element (annual mean value) is used.

    Returns
    -------
    complex
        Complex eigenvalue of the BWJ system:

        real part:
            BJ (growth rate, 1/month)
        imaginary part:
            WF (oscillation frequency, 1/month)

    Examples
    --------
    >>> par = {
    ...     "R": [0.5],
    ...     "F1": [1.2],
    ...     "F2": [1.0],
    ...     "epsilon": [0.3]
    ... }
    >>> RO_BWJ(par)
    (0.1+0.95j)
    """
    
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


def RO_analytic_std(par):
    """
    Compute analytical standard deviation of T and h for the Recharge Oscillator (RO) model.
    ONLY for linear RO model with white noise (annual mean parameters only)

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

