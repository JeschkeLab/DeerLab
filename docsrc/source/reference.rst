.. _api_reference:

API Reference
============================================================


.. currentmodule:: deerlab

.. rubric:: Classes

.. autosummary::
    :toctree: _autosummary
    :template: custom_class_template.rst
    :nosignatures:
    
    Model
    Parameter
    Penalty
    FitResult
    UQResult


.. currentmodule:: deerlab

.. rubric:: Functions

.. autosummary::
    :toctree: _autosummary
    :template: custom_function_template.rst
    :nosignatures:

    fit
    merge
    link 
    lincombine
    relate 
    regoperator 
    selregparam 
    noiselevel 
    whitegaussnoise
    diststats
    correctphase
    bootstrap_analysis
    profile_analysis 
    snlls 
    fnnls
    cvxnnls
    
.. rubric:: Dipolar EPR functions

.. autosummary::
    :toctree: _autosummary
    :template: custom_function_template.rst
    :nosignatures:

    dipolarmodel 
    dipolarpenalty
    deerload
    dipolarkernel 
    dipolarbackground
    fftspec
    distancerange

.. rubric:: Utility functions 

All functions in ``deerlab.utils`` are private functions used internally by other functions.
Stable functionality is not guaranteed.

.. currentmodule:: deerlab.utils

.. autosummary::
    :toctree: _autosummary
    :template: custom_function_template.rst
    :nosignatures:

    store_pickle
    read_pickle
    hccm
    Jacobian
    nearest_psd
    movmean
    ovl
    der_snr
    
