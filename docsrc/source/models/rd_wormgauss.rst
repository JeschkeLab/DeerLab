.. highlight:: matlab
.. _rd_wormgauss:

***********************
:mod:`rd_wormgauss`
***********************

Worm-like chain model near the rigid limit with Gaussian convolution

Syntax
=========================================

.. code-block:: matlab

        info = rd_wormgauss()
        P = rd_wormgauss(r,param)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)

============== =============== ======== ======== ======== ===============================
 Variable       Symbol         Default   Lower   Upper       Description
============== =============== ======== ======== ======== ===============================
``param(1)``   :math:`L`       3.7      1.5       10        Length of the worm-like chain
``param(2)``   :math:`L_p`     10       2         100       Persistence length
``param(3)``   :math:`\sigma`  0.2      0.01      2         Gaussian standard deviation
============== =============== ======== ======== ======== ===============================

Example using default parameters:

.. image:: ../images/model_rd_wormgauss.png
   :width: 40%


Description
=========================================

.. code-block:: matlab

        info = rd_wormgauss()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = rd_wormgauss(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

