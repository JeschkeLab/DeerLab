.. highlight:: matlab
.. _mixmodels:


***********************
:mod:`mixmodels`
***********************

Combine parametric models into one

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`newModel = mixmodels(Models)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **Models** - Input parametric models (cell array of function handles)

Returns
    *   **newModel** - Generated mixed parameteric model (function handle)

Usage
=========================================

.. code-block:: matlab

            newmodel = mixmodels({@model1,@model2, .., @modelN})

Combines the parametric model function handles ``@model1``, ``@model2, ...,  ``@modelN`` into a new parametric model function handle ``newmodel``. The models must be passed as a cell array of function handles and the parametric models must be of the type as the models distributed in DeerAnalysis2. The returned function handle can be used for parametric model fitting.

Example:

    .. code-block:: matlab

        myModel = mixmodels({@onegaussian,@tworice})