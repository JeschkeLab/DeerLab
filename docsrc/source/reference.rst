.. _api_reference:

Reference Index
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
    ExperimentInfo


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
    goodness_of_fit
    
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

.. currentmodule:: deerlab

.. autosummary::
    :toctree: _autosummary
    :template: custom_function_template.rst
    :nosignatures:

    store_pickle
    read_pickle
    sophegrid
    choleskycovmat
    hccm
    Jacobian
    nearest_psd
    movmean
    ovl
    der_snr
    formatted_table
