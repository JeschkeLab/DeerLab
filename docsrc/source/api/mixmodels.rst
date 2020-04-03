.. highlight:: matlab
.. _mixmodels:


***********************
:mod:`mixmodels`
***********************

Combine parametric models into one

-----------------------------



Syntax
=========================================

.. code-block:: matlab

    newmodel = mixmodels(models)



Parameters
    *   ``models`` - Input parametric models (*N*-element cell array of function handles)

Returns
    *   ``newModel`` - Generated mixed parametric model (function handle)

-----------------------------



Description
=========================================

.. code-block:: matlab

            newmodel = mixmodels({@model1,@model2, .., @modelN})

Combines the parametric model function handles ``model1``, ``model2``, ...,  ``modelN`` into a new parametric model function handle ``newmodel``. The models must be passed as a cell array of function handles and the parametric models must be of the type as the models distributed in DeerLab. The returned function handle can be used for parametric model fitting.

Information about the new model and its parameters can be accessed by

.. code-block:: matlab

            info = newmodel();

The mixed models will include an additional amplitude parameter for each additional model in the mix.

-----------------------------



Example
=========================================

If one mixes a single Gaussian model (2 parameters) with a WLC model (2 parameters) into a single model

.. code-block:: matlab

    newmodel = mixmodels({@dd_gauss,@dd_wormchain})

the resulting model ``newmodel`` will contain 5 parameters in the following order: 1 amplitude parameter, the 2 single-Gaussian parameters and the 2 WLC parameters. 
