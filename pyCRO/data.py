import re
import os
import sys

import numpy as np
import xarray as xr

from importlib.resources import files

_MAT_FILENAME = str(files("pyCRO").joinpath("data", "CRO_parlib_v0.0.mat"))
_CESM1_FILENAME = files("pyCRO").joinpath("data", "CESM1_LENS_ENSO_timeseries.nc")
_CMIP6_FILENAME = files("pyCRO").joinpath("data", "CROdata_timeseries_CMIP6.nc")
_ORAS5_FILENAME = files("pyCRO").joinpath("data", "CROdata_timeseries_oras5.nc")

def _try_load_mat(fname):
    """Load .mat (v7 via scipy; fallback to v7.3 via mat73)."""
    try:
        from scipy.io import loadmat
        mat = loadmat(fname, squeeze_me=True, struct_as_record=False)
        # strip MATLAB metadata keys
        return {k: v for k, v in mat.items() if not k.startswith("__")}
    except Exception:
        # v7.3 (HDF5) fallback if mat73 is available
        try:
            import mat73
            return mat73.loadmat(fname)
        except Exception as e:
            raise RuntimeError(f"Failed to load {fname} with scipy and mat73") from e

def _to_str_array(x):
    """Normalize MATLAB string/cellstr/char arrays -> numpy array of Python str (same shape)."""
    x = np.asarray(x, dtype=object)
    out = np.empty(x.shape, dtype=object)
    it = np.nditer(out, flags=['multi_index', 'refs_ok'], op_flags=['writeonly'])
    while not it.finished:
        v = x[it.multi_index]
        if isinstance(v, str):
            s = v
        elif isinstance(v, bytes):
            s = v.decode("utf-8", errors="ignore")
        elif isinstance(v, np.ndarray) and v.dtype.kind in ("U", "S"):
            # MATLAB char array -> join characters
            s = "".join(map(str, v.tolist()))
        else:
            s = str(v)
        it[0] = s
        it.iternext()
    return out.astype(str)

def _as_col_cell(obj):
    """
    Ensure a (16,1) numpy object array from a MATLAB 16x1 cell stored inside S['par'] element.
    The element could already be object array (16,), (16,1), list of 16, etc.
    """
    arr = obj
    # Convert lists/tuples to np.object array
    if isinstance(arr, (list, tuple)):
        arr = np.array(arr, dtype=object)
    if isinstance(arr, np.ndarray):
        # flatten then reshape to (16,1)
        arr = arr.astype(object)
        arr = arr.reshape(-1, order="F")  # MATLAB-friendly flatten
        if arr.size != 16:
            raise ValueError(f"Expected 16 elements, got {arr.size}")
        return arr.reshape(16, 1, order="F")
    # Anything else: treat as scalar and fail
    raise TypeError("Unexpected parameter cell content type")


def par_load(data_name: str, ro_name: str):
    """
    Load precomputed CRO parameters from `CRO_parlib_v0.0.mat`.

    This function extracts parameter sets for different RO model configurations
    and datasets (e.g., CMIP6 historical runs, ORAS5, etc.).

    See details at `precomputed Library <../notebooks/model_library.html>`_

    Parameters
    ----------    
    data_name : str
        Dataset identifier. Supported options include:
        
            - "ORAS5"
            
            - "CMIP6-historical-1" ... "CMIP6-historical-48"
            
            - "CMIP6-historical-all" (returns all realizations)

    ro_name : str
        Recharge Oscillator model configuration. Supported options:
        
            - "Linear-White-Additive"
            
            - "Seasonal-Linear-White-Additive"
            
            - "Nonlinear-White-Additive"
            
            - "Seasonal-Nonlinear-White-Additive"
            
            - "Linear-White-Multiplicative"
            
            - "Seasonal-Linear-White-Multiplicative"
            
            - "Nonlinear-White-Multiplicative"
            
            - "Seasonal-Nonlinear-White-Multiplicative"

    Returns
    -------
    dict or list of dict
        If `data_name != "CMIP6-historical-all"`:
            Returns a dictionary containing CRO parameters:
                - R, F1, F2, epsilon, b_T, c_T, d_T, b_h
                - sigma_T, sigma_h, B
                - m_T, m_h, n_T, n_h, n_g

            Each value is a list (to support scalar or ensemble formats).

        If `data_name == "CMIP6-historical-all"`:
            Returns a list of parameter dictionaries, one per CMIP6 realization.

    Raises
    ------
    ValueError
        If `ro_name` is not recognized or `data_name` is invalid.

    Notes
    -----
    - The parameter library is stored in MATLAB file format (`.mat`).
    - Internally, parameters are indexed as:
        (RO configuration index, dataset index).
    - This function automatically handles scalar, vector, and empty entries
      from MATLAB cell arrays.

    Examples
    --------
    >>> par = par_load("ORAS5", "Linear-White-Additive")
    >>> par["R"]

    >>> pars = par_load("CMIP6-historical-all", "Nonlinear-White-Additive")
    >>> len(pars)
    """
    S = _try_load_mat(_MAT_FILENAME)
    par = S['par']
    # print(par.shape) # (8, 49)

    data_name = str(data_name)
    ro_name   = str(ro_name)


    # ---- exact single match ----
    if ro_name == "Linear-White-Additive":
        ro_name_index = 0
    elif ro_name == "Seasonal-Linear-White-Additive":
        ro_name_index = 1
    elif ro_name == "Nonlinear-White-Additive":
        ro_name_index = 2
    elif ro_name == "Seasonal-Nonlinear-White-Additive":
        ro_name_index = 3
    elif ro_name == "Linear-White-Multiplicative":
        ro_name_index = 4
    elif ro_name == "Seasonal-Linear-White-Multiplicative":
        ro_name_index = 5
    elif ro_name == "Nonlinear-White-Multiplicative":
        ro_name_index = 6
    elif ro_name == "Seasonal-Nonlinear-White-Multiplicative":
        ro_name_index = 7
    else:
        raise ValueError(f"Wrong input for RO_type")


    if data_name == "CMIP6-historical-all":
        my_parr = []
        for data_name_index in range(1,49):
            row = par[ro_name_index, data_name_index]

            R, F1, F2, epsilon, b_T, c_T, d_T, b_h, sigma_T, sigma_h, B, m_T, m_h, n_T, n_h, n_g = row

            def to_list(x):
                if isinstance(x, np.ndarray):
                    if x.size == 0:          # empty placeholder
                        return []
                    if x.ndim == 0:          # scalar-like
                        return [x.item()]
                    return x.ravel().tolist()
                # numpy or python scalar
                try:
                    return [x.item()]
                except AttributeError:
                    return [x] if np.isscalar(x) else [x]

            my_par = {
                'R': to_list(R),
                'F1': to_list(F1),
                'F2': to_list(F2),
                'epsilon': to_list(epsilon),
                'b_T': to_list(b_T),
                'c_T': to_list(c_T),
                'd_T': to_list(d_T),
                'b_h': to_list(b_h),
                'sigma_T': to_list(sigma_T),
                'sigma_h': to_list(sigma_h),
                'B': to_list(B),
                'm_T': to_list(m_T),
                'm_h': to_list(m_h),
                'n_T': to_list(n_T),
                'n_h': to_list(n_h),
                'n_g': to_list(n_g),
            }
            my_parr.append(my_par)
        return my_parr
    else:
        data_name_mapping = {"ORAS5": 0}
        data_name_mapping.update({f"CMIP6-historical-{i}": i for i in range(1, 49)})

        try:
            data_name_index = data_name_mapping[data_name]
        except KeyError:
            raise ValueError("Invalid input for `data_name` or `ro_name`")

        row = par[ro_name_index, data_name_index]

        R, F1, F2, epsilon, b_T, c_T, d_T, b_h, sigma_T, sigma_h, B, m_T, m_h, n_T, n_h, n_g = row

        def to_list(x):
            if isinstance(x, np.ndarray):
                if x.size == 0:          # empty placeholder
                    return []
                if x.ndim == 0:          # scalar-like
                    return [x.item()]
                return x.ravel().tolist()
            # numpy or python scalar
            try:
                return [x.item()]
            except AttributeError:
                return [x] if np.isscalar(x) else [x]

        my_par = {
            'R': to_list(R),
            'F1': to_list(F1),
            'F2': to_list(F2),
            'epsilon': to_list(epsilon),
            'b_T': to_list(b_T),
            'c_T': to_list(c_T),
            'd_T': to_list(d_T),
            'b_h': to_list(b_h),
            'sigma_T': to_list(sigma_T),
            'sigma_h': to_list(sigma_h),
            'B': to_list(B),
            'm_T': to_list(m_T),
            'm_h': to_list(m_h),
            'n_T': to_list(n_T),
            'n_h': to_list(n_h),
            'n_g': to_list(n_g),
        }

        return my_par


def ROdata_load(name: str):
    """
    Load precomputed CRO dataset time series.

    This function provides access to built-in ENSO-related datasets used in
    the CRO (Coupled Recharge Oscillator) framework, including ORAS5 reanalysis,
    CMIP6 historical simulations, and CESM1 LENS. 

    Parameters
    ----------
    name : str
        Name of the dataset to load.

        Supported options:

        - "CESM1_LENS"
            CESM1 Large Ensemble Niño3.4 and thermocline time series.

        - "CMIP6"
            Preprocessed CMIP6 historical ENSO time series (multi-model mean or ensemble form).

        - "ORAS5"
            ORAS5 ocean reanalysis ENSO time series.

    Returns
    -------
    str
        Absolute file path to the requested NetCDF dataset.

    Raises
    ------
    ValueError
        If `data_name` is not one of the supported options.

    Notes
    -----
    - Data are stored inside the `pyCRO.data` package directory.
    - Files are in NetCDF format and should be opened using `xarray.open_dataset`.
    - This function does not load the dataset into memory, only returns the path.

    Examples
    --------
    Load CESM1 LENS data:

    >>> ds = pyCRO.ROdata_load("CESM1_LENS")

    Load ORAS5 reanalysis:

    >>> ds = pyCRO.ROdata_load("ORAS5")
    """

    mapping = {
        "CESM1_LENS": _CESM1_FILENAME,
        "CMIP6": _CMIP6_FILENAME,
        "ORAS5": _ORAS5_FILENAME,
    }
    
    try:
        return xr.open_dataset(mapping[name])
    except KeyError:
        raise ValueError(
            f"Invalid dataname='{name}'. "
            f"Choose from {list(mapping.keys())}."
        )


##################################################################################

def ROdata_calc(
    sst_a: xr.DataArray,
    h_a: xr.DataArray,
    sst_regions=None,
    h_regions=None,
) -> xr.Dataset:
    """
    Compute Recharge Oscillator SST and thermocline/heat-content indices.

    Parameters
    ----------
    sst_a : xr.DataArray
        Sea surface temperature with dimensions ``(time, lat, lon)``.

    h_a : xr.DataArray
        Thermocline depth, heat content, SSH, or another recharge proxy with
        dimensions ``(time, lat, lon)``.

    sst_regions : list of str or dict, optional
        SST regions to compute.

        * ``None`` (default): compute all predefined SST regions.
        * list of str: predefined region names (e.g., ``["Nino34", "Nino3"]``).
        * dict: mapping from output names to custom region definitions.

    h_regions : list of str or dict, optional
        Thermocline/heat-content regions to compute.

        * ``None`` (default): compute all predefined thermocline regions.
        * list of str: predefined region names.
        * dict: mapping from output names to custom region definitions.

    Returns
    -------
    xr.Dataset
        Dataset containing the requested SST and thermocline indices.

    Examples
    --------
    Compute all predefined SST and thermocline indices::

        ds = ROdata_calc(sst, h)

    Compute only the Niño-3.4 SST index::

        ds = ROdata_calc(sst, h, sst_regions=["Nino34"])

    Compute selected thermocline indices::

        ds = ROdata_calc(sst, h, h_regions=["Hw", "He"])

    Define a custom SST region::

        ds = ROdata_calc(
            sst,
            h,
            sst_regions={
                "CP": {
                    "latS": -5,
                    "latN": 5,
                    "lonW": 170,
                    "lonE": 220,
                }
            },
        )
    """
    DEFAULT_SST_REGIONS = ["Nino34", "Nino3", "Nino4", "Nino12", "ColdTongue",]
    DEFAULT_H_REGIONS = ["Hm", "Hw", "He", "Hw1", "He1", "Hw2", "He2",]

    if sst_regions is None:
        sst_regions = DEFAULT_SST_REGIONS

    if h_regions is None:
        h_regions = DEFAULT_H_REGIONS

    # -----------------------------------------------------
    # Area-mean SST indices
    # -----------------------------------------------------
    ssti = area_average_regions(sst_a, sst_regions)

    # -----------------------------------------------------
    # Area-mean thermocline/heat content indices
    # -----------------------------------------------------
    hi = area_average_regions(h_a, h_regions)

    # -----------------------------------------------------
    # Align time
    # -----------------------------------------------------
    ssti, hi = xr.align(ssti, hi, join="inner")

    return xr.merge([ssti, hi])


##################################################################################

def area_average(x, region=None):
    '''
        cos-weighted area averaged fields
    '''
    if region is None:
        x_subset = x
    else:
        x_subset = _select_region(x, region)

    w = np.cos(np.deg2rad(x_subset.lat))
    w = w.broadcast_like(x_subset)
    aave = (x_subset * w).mean(dim=('lat', 'lon'))/w.mean(dim=('lat', 'lon'))

    if isinstance(region, dict):
        region_name = region.get('name') or region.get('long_name')
        if region_name:
            aave.attrs['long_name'] = region_name
    return aave


def area_average_regions(x, regions):
    """
    Compute cosine-weighted regional averages.

    Parameters
    ----------
    x : xr.DataArray
        Input field.

    regions : list of str or dict
        Regions to average.

        * list of str
            Predefined region names.

        * dict
            Mapping from output variable names to region definitions.

    Returns
    -------
    xr.Dataset
        Dataset containing one variable per region.
    """
    if isinstance(regions, (list, tuple)):
        region_dict = {name: name for name in regions}

    elif isinstance(regions, dict):
        region_dict = regions

    else:
        raise TypeError(
            "regions must be a list of region names or a dictionary."
        )

    out = []

    for name, region in region_dict.items():
        da = area_average(x, region)
        da.name = name
        out.append(da.to_dataset())

    return xr.merge(out)


def _select_region(x, region):
    '''
        select region boxes
    '''
    if isinstance(region, str):
        try:
            Rbox = _box_region_array(region)
            R_latS = Rbox['latS']
            R_latN = Rbox['latN']
            R_lonW = Rbox['lonW']
            R_lonE = Rbox['lonE']
            x_sel = x.sel(lat=slice(R_latS, R_latN), lon=slice(R_lonW, R_lonE))

        except ValueError:
            print("box_region_array error: undefined region string!")

    elif isinstance(region, dict) and len(region)>=4:
        R_latS = region['latS']
        R_latN = region['latN']
        R_lonW = region['lonW']
        R_lonE = region['lonE']
        x_sel = x.sel(lat=slice(R_latS, R_latN), lon=slice(R_lonW, R_lonE))

    elif isinstance(region, (list, tuple))  and len(region)>=4:
        R_lonW = region[0]
        R_lonE = region[1]
        R_latS = region[2]
        R_latN = region[3]
        x_sel = x.sel(lat=slice(R_latS, R_latN), lon=slice(R_lonW, R_lonE))
    else:

        raise ValueError('error in select_region: unsupported region type!')

    return x_sel

def _box_region_array(Rstr = 'Nino34'):
    '''
        return box region array of latS, latN, lonW, lonE
    '''
    region = {}

    ## ENSO SST indices
    region['Nino34']  = {'latS': -5, 'latN': 5, 'lonW': 190, 'lonE': 240, 'name': 'Niño3.4 (5°S-5°N, 170°W-120°W)'}
    region['Nino3']   = {'latS': -5, 'latN': 5, 'lonW': 210, 'lonE': 270, 'name': 'Niño3 (5°S-5°N, 150°W-90°W)'}
    region['Nino4']   = {'latS': -5, 'latN': 5, 'lonW': 160, 'lonE': 210, 'name': 'Niño4 (5°S-5°N, 160°E-150°W)'}
    region['Nino12']  = {'latS': -10, 'latN': 0, 'lonW': 270, 'lonE': 280, 'name': 'Niño1+2 (10°S-0°, 90°W-80°W)'}
    region['ColdTongue'] = {'latS': -6, 'latN': 6, 'lonW': 180, 'lonE': 270, 'name': 'ColdTongue (6°S-6°N, 180°-90°W)'}

    ## ENSO thermcline/heat content indices regions
    region['Hm'] = {'latS': -5, 'latN': 5, 'lonW': 120, 'lonE': 280, 'name': 'Hm (5°S-5°N, 120°E-80°W)'}

    ## origioanl defintion: https://www.pmel.noaa.gov/elnino/upper-ocean-heat-content-and-enso
    region['Hw'] = {'latS': -5, 'latN': 5, 'lonW': 120, 'lonE': 205, 'name': 'Hw (5°S-5°N, 120°E-155°W)'}
    region['He'] = {'latS': -5, 'latN': 5, 'lonW': 205, 'lonE': 280, 'name': 'He (5°S-5°N, 155°W-80°W)'}

    ## defintion in Zhao et al. (2021), west and east are the same size http://onlinelibrary.wiley.com/doi/abs/10.1029/2021GL094366
    region['Hw1'] = {'latS': -5, 'latN': 5, 'lonW': 120, 'lonE': 200, 'name': 'Hw1 (5°S-5°N, 120°E-160°W)'}
    region['He1'] = {'latS': -5, 'latN': 5, 'lonW': 200, 'lonE': 280, 'name': 'He1 (5°S-5°N, 160°W-80°W)'}

    ## this defintion may works for SSH data, see details in Zhao et al. (2021)
    region['Hw2'] = {'latS': -5, 'latN': 5, 'lonW': 120, 'lonE': 180, 'name': 'Hw2 (5°S-5°N, 120°E-180°)'}
    region['He2'] = {'latS': -5, 'latN': 5, 'lonW': 180, 'lonE': 280, 'name': 'He2 (5°S-5°N, 180°-80°W)'}

    return region.get(Rstr, "nothing")