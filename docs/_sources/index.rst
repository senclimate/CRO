************************************
CRO: Community Recharge Oscillator
************************************

:code:`CRO` is an open-source Python and MATLAB package for modeling the El Niño–Southern Oscillation (ENSO) using the recharge oscillator (RO) framework (Jin).  It provides tools for solving the RO equations, fitting parameters to data, and applying the model in research and forecasting.

**Key Features**

- **Solver**: Numerical and analytical solvers for the recharge oscillator equations.  
- **Fitting**: Parameter estimation from observations, reanalysis, or climate model data.  
- **Applications**: simulations, sensitivity experiments, and analysis of ENSO changes.  

.. warning ::
    This package is still in its early stage and under active development, and its API could be changed frequently.

|

.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card::  Installation
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/installation.png
        :link: ug-installation
        :link-type: doc

        Installation instructions.

    .. grid-item-card::  Model
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/setup.png
        :link: ug-model
        :link-type: doc

        :code:`CRO` model and parameters 

    .. grid-item-card::  Solver
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/setup.png
        :link: ug-solver
        :link-type: doc

        :code:`CRO` solver 

    .. grid-item-card::  Fitting
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/postprocessing.png
        :link: ug-fitting
        :link-type: doc

        :code:`CRO` fitting 

    .. grid-item-card::  Application
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/diags.png
        :link: ug-application
        :link-type: doc

        :code:`CRO` applications 

    .. grid-item-card::  API
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/api.png
        :link: ug-api
        :link-type: doc

        The essential API.

.. toctree::
   :maxdepth: 3
   :hidden:
   :caption: User Guide

   ug-installation
   ug-solver
   ug-fitting
   ug-application
   ug-api