
import numpy as np
import matplotlib.pyplot as plt

def _seasonal_cycle(param):
    """
    Compute the monthly seasonal cycle (12 points) from parameter vector.
    
    Parameters
    ----------
    param : array-like
        - length 1: [R0]
        - length 3: [R0, R_A, R_pha]
        - length 5: [R0, R_A, R_pha, R_A2, R_pha2]
        - length 0: return zeros
    
    Returns
    -------
    ndarray
        Seasonal cycle values for 12 months.
    """
    t = np.arange(0.5, 12, step=1)  # month centers
    if len(param) == 1:
        return np.full(12, param[0])
    elif len(param) == 3:
        R0, R_A, R_pha = param
        return R0 + R_A * np.sin(2*np.pi/12*t + R_pha)
    elif len(param) == 5:
        R0, R_A, R_pha, R_A2, R_pha2 = param
        return (R0
                + R_A  * np.sin(2*np.pi/12*t + R_pha)
                + R_A2 * np.sin(2*np.pi/6*t  + R_pha2))
    elif len(param) == 0:
        return np.zeros(12)
    else:
        raise ValueError("Parameter must have length 1, 3, or 5.")


# Mapping parameter -> LaTeX label
_PARAM_LABELS = {
    "BJ":       r"BJ (month$^{-1}$)",  
    "R":       r"$R$ (month$^{-1}$)",  
    "F1":      r"$F_1$ (K m$^{-1}$ month$^{-1}$)",  
    "F2":      r"$F_2$ (m K$^{-1}$ month$^{-1}$)",  
    "epsilon": r"$\varepsilon$ (month$^{-1}$)",  
    "b_T":     r"$b_T$ (K$^{-1}$ month$^{-1}$)",  
    "c_T":     r"$c_T$ (K$^{-2}$ month$^{-1}$)",  
    "d_T":     r"$d_T$ (m$^{-1}$ month$^{-1}$)",  
    "b_h":     r"$b_h$ (K$^{-2}$ m month$^{-1}$)",  
    "sigma_T": r"$\sigma_T$ (K month$^{-1/2}$ or K month$^{-1}$)",  
    "sigma_h": r"$\sigma_h$ (m month$^{-1/2}$ or m month$^{-1}$)",  
    "B":       r"$B$ (K$^{-1}$)",  
}


def plot_RO_par(par, ax=None, keys=None, ncol=4, label=None):
    """
    Plot deterministic RO parameter seasonal cycles.
    
    Parameters
    ----------
    par : dict
        Dictionary of fitted CRO parameters. Each key (e.g., 'R', 'F1', etc.)
        maps to an array of length 1, 3, or 5.
    ax : matplotlib axis or array-like, optional
        Axis or array of axes to plot on. If None, creates a new figure/axes.
    keys : list of str, optional
        Which parameters to plot. If None, uses a default set.
    ncol : int, default=4
        Number of columns for subplot layout.
    """
    if keys is None:
        keys = ['BJ', 'R', 'epsilon', 'F1', 'F2', 'b_T', 'c_T', 'd_T', 'b_h']
    n = len(keys)
    nrow = int(np.ceil(n / ncol))

    if ax is None:
        fig, axs = plt.subplots(nrow, ncol, figsize=(3.5*ncol, 2*nrow), layout='compressed')
        axs = np.array(axs).reshape(-1)
    else:
        axs = np.atleast_1d(ax)

    for i, k in enumerate(keys):
        if k == 'BJ':
            sea_R = _seasonal_cycle(np.atleast_1d(par.get('R', [])))
            sea_Eps = _seasonal_cycle(np.atleast_1d(par.get('Eps', [])))
            seasonal = (sea_R - sea_Eps)/2
        else:
            seasonal = _seasonal_cycle(np.atleast_1d(par.get(k, [])))

        if label is None:
            lb = k
        else:
            lb = label
            
        axs[i].plot(np.arange(1, 13), seasonal, marker=".", label=lb)
        axs[i].set_ylabel(_PARAM_LABELS[k])
        axs[i].set_xticks(range(1, 13))
        # axs[i].set_xticklabels(["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"])
        axs[i].set_xticklabels(["J","F","M","A","M","J","J","A","S","O","N","D"])
        axs[i].grid(True, linestyle="--", alpha=0.5)
        # axs[i].legend()

    # Hide unused subplots if n < nrow*ncol
    for j in range(n, len(axs)):
        axs[j].axis("off")

    return axs


def plot_ens_RO_par(
    par,
    ax=None,
    keys=None,
    ncol=4,
    label=None,
    spread=None,
    center="mean",
    show_members=False,
    alpha=0.25,
    member_alpha=0.08,
):
    """
    Plot ensemble RO parameters.

    Parameters
    ----------
    par : dict or list of dict
        Either

        - fitted parameter dictionary (single fit)

        or

        - list of fitted parameter dictionaries returned by multiple
          ``RO_fitting`` calls.

    spread : {None, "std", "q10_q90", "iqr"}, optional
        Ensemble spread to plot.

    show_members : bool
        Plot all ensemble members.
    """

    if keys is None:
        keys = ['BJ', 'R', 'epsilon', 'F1', 'F2',
                'b_T', 'c_T', 'd_T', 'b_h']

    n = len(keys)
    nrow = int(np.ceil(n / ncol))

    if ax is None:
        fig, axs = plt.subplots(
            nrow,
            ncol,
            figsize=(3.5*ncol, 2*nrow),
            layout="compressed",
        )
        axs = np.asarray(axs).reshape(-1)
    else:
        axs = np.atleast_1d(ax)

    months = np.arange(1,13)

    ##################################################################
    # Single fit
    ##################################################################

    if isinstance(par, dict):

        for i, k in enumerate(keys):

            if k == "BJ":
                seasonal = (
                    _seasonal_cycle(np.atleast_1d(par["R"]))
                    - _seasonal_cycle(np.atleast_1d(par["epsilon"]))
                ) / 2
            else:
                seasonal = _seasonal_cycle(
                    np.atleast_1d(par.get(k, []))
                )

            axs[i].plot(
                months,
                seasonal,
                marker=".",
                lw=2,
                label=label if label else k,
            )

            axs[i].set_ylabel(_PARAM_LABELS[k])
            axs[i].set_xticks(months)
            axs[i].set_xticklabels(
                ["J","F","M","A","M","J",
                 "J","A","S","O","N","D"]
            )
            axs[i].grid(True, ls="--", alpha=0.5)

    ##################################################################
    # Ensemble fit
    ##################################################################

    else:

        fit_results = par

        for i, k in enumerate(keys):

            ax = axs[i]

            cycles = []

            for fit in fit_results:

                if k == "BJ":
                    cyc = (
                        _seasonal_cycle(np.atleast_1d(fit["R"]))
                        - _seasonal_cycle(np.atleast_1d(fit["epsilon"]))
                    ) / 2
                else:
                    cyc = _seasonal_cycle(
                        np.atleast_1d(fit.get(k, []))
                    )

                cycles.append(cyc)

            cycles = np.asarray(cycles)

            if center == "mean":
                yc = np.nanmean(cycles, axis=0)
            else:
                yc = np.nanmedian(cycles, axis=0)

            if show_members:
                ax.plot(
                    months,
                    cycles.T,
                    color="0.75",
                    lw=0.8,
                    alpha=member_alpha,
                )

            if spread == "std":
                ys = np.nanstd(cycles, axis=0)

                ax.fill_between(
                    months,
                    yc-ys,
                    yc+ys,
                    alpha=alpha,
                )

            elif spread == "iqr":

                q25 = np.nanquantile(cycles,0.25,axis=0)
                q75 = np.nanquantile(cycles,0.75,axis=0)

                ax.fill_between(
                    months,
                    q25,
                    q75,
                    alpha=alpha,
                )

            elif spread == "q10_q90":

                q10 = np.nanquantile(cycles,0.10,axis=0)
                q90 = np.nanquantile(cycles,0.90,axis=0)

                ax.fill_between(
                    months,
                    q10,
                    q90,
                    alpha=alpha,
                )

            ax.plot(
                months,
                yc,
                marker=".",
                lw=2,
                label=label if label else k,
            )

            ax.set_ylabel(_PARAM_LABELS[k])
            ax.set_xticks(months)
            ax.set_xticklabels(
                ["J","F","M","A","M","J",
                 "J","A","S","O","N","D"]
            )
            ax.grid(True, ls="--", alpha=0.5)

    for j in range(n, len(axs)):
        axs[j].axis("off")

    return axs