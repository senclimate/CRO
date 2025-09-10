************************************
CRO: Community Recharge Oscillator
************************************

What is CRO?
============

The recharge oscillator (RO) model is one of the leading theories for the El Niño–Southern Oscillation (ENSO) :cite:`jin1997, jin2020, vialard2025`. While the literature contains many RO variants and implementations, :code:`CRO` is an open-source Python and Matlab package for solving the RO equations, fitting parameters to observational or model data, and applying the model in teaching and research :cite:`kim2025`. 

Key Features
------------

- **Solver**: Numerical and analytical solvers for the RO equations.  
- **Fitting**: Parameter estimation from observations, reanalysis, or climate model data.  
- **Applications**: Simulations, sensitivity experiments, and analysis of ENSO dynamics.  


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
        :img-top: _static/model.png
        :link: ug-model
        :link-type: doc

        :code:`CRO` equations and parameters 

    .. grid-item-card::  Solver
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/solver.png
        :link: ug-solver
        :link-type: doc

        :code:`CRO` solver 

    .. grid-item-card::  Fitting
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/fitting.png
        :link: ug-fitting
        :link-type: doc

        :code:`CRO` fitting 

    .. grid-item-card::  Application
        :class-title: custom-title
        :class-body: custom-body
        :img-top: _static/application.png
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




Relationship with XRO
=====================

The eXtended nonlinear Recharge Oscillator :code:`XRO` (https://github.com/senclimate/XRO) was developed to study interactions between ENSO and other climate modes, providing advanced predictive capabilities :cite:`zhao2024`.

:code:`CRO`, in contrast, focuses on the core features and flexibility of the RO framework. It offers multiple fitting methods and supports varying degrees of seasonality in individual parameters.

Both models share the same theoretical RO foundation and belong to the same RO family.


References
==========

.. bibliography::
   :style: apa