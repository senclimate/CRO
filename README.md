# CRO: Community Recharge Oscillator Model

## What is CRO?

The recharge oscillator (RO) model is one of the leading theoretical frameworks for understanding the El Niño–Southern Oscillation (Jin, 1997, Jin et al., 2020, Vialard et al., 2025).

While many RO variants exist in the literature, **CRO** is an open-source Python/Matlab package for:

- solving the recharge oscillator equations  
- fitting model parameters to observational and climate model data  
- applying the model in research and teaching applications

Key methodological foundations are described in Kim et al. (2025).

---

## Key Features

- **Solver**: Numerical and analytical solvers for the RO equations.

- **Fitting**: Parameter estimation from observations, reanalysis, or climate model data.

- **Applications**: Simulations, sensitivity experiments, and analysis of ENSO dynamics.

---

## Installation

```bash
pip install pythonCRO
```

---

## Quick Start (Python)

```python
import pyCRO

print(pyCRO.__version__)
```

---

## Data Included

CRO includes curated datasets for benchmarking and teaching:

- Precomputed ORAS5 ENSO SST and WWV indices
- Precomputed CESM1 Large Ensemble ENSO time series  
- Precomputed CMIP6 ENSO time series
- Precomputed parameter libraries  

---

## Documentation

- Documentation: https://senclimate.github.io/CRO/python/  
- Source code: https://github.com/senclimate/CRO  

---

## Citation

If you use CRO in your research, please cite:

Kim, S.-K., Zhao, S., & et al. (2025). Community Recharge Oscillator (CRO) v1.0: an open-source Python and MATLAB package for solving, parameter fitting, and practical applications of the ENSO recharge oscillator. In preparation.
