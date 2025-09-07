************************************
CRO: Community Recharge Oscillator
************************************

:code:`CRO` is an open-source Python and MATLAB package for modeling the El Niño–Southern Oscillation (ENSO) using the recharge oscillator (RO) framework (see `Jin, 1997 <#ref-jin1997>`_; `Zhao et al., 2024 <#ref-zhao2024>`_; 
`Kim et al., in preparation <#ref-cro>`_).  It provides tools for solving the RO equations, fitting parameters to data, and applying the model in teaching, research and even forecasting.

**Key Features**

- **Solver**: Numerical and analytical solvers for the RO equations.  
- **Fitting**: Parameter estimation from observations, reanalysis, or climate model data.  
- **Applications**: Simulations, sensitivity experiments, and analysis of ENSO changes.  

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
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   ug-installation
   ug-model
   ug-solver
   ug-fitting
   ug-application
   ug-api


.. _references:

References
==========

.. _ref-jin1997:

Jin, F.-F. (1997). An equatorial ocean recharge paradigm for ENSO. 
*Part I: Conceptual model*. Journal of the Atmospheric Sciences, 54(7), 811–829. 
https://doi.org/10.1175/1520-0469(1997)054<0811:AEORPF>2.0.CO;2  

.. _ref-zhao2024:

Zhao, S., Stevenson, S., Kim, S.-K., Fedorov, A. V., & Vecchi, G. A. (2024).  
Pantropical interactions and predictability of ENSO in a nonlinear recharge oscillator framework. 
*Nature*, 626, 543–549. https://doi.org/10.1038/s41586-024-07009-1  

.. _ref-cro:

Kim, S.-K., Zhao, S., et al. (in preparation).  
*Community Recharge Oscillator (CRO) v1.0: an open-source Python and MATLAB package for solving, parameter fitting, and practical applications of the ENSO recharge oscillator.*  
