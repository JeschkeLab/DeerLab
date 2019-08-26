.. highlight:: matlab
.. _randcoil:

***********************
:mod:`randcoil`
***********************

Random-coil model for an unfolded peptide/protein

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = randcoil(t,param)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **r** - Distance axis (N-array)
    *   **param** - Model parameters
Returns
    *   **P** - Model distance distribution (N-array)

Model equation: :math:`P(r) = \frac{3}{(2\pi\nu_0)^{3/2}}4\pi r^2\exp(-\frac{3 r^2}{\nu_0})`

where :math:`\nu_0 = 3/(12\pi r_0 N \nu)^{3/2}` and :math:`r_0` = 0.198 nm.

========== =========== ======== ======== ======== ==================================
 Variable   Symbol     Default   Lower   Upper       Description
========== =========== ======== ======== ======== ==================================
param(1)   :math:`N`      50      2        1000    Number of residues between labels
param(2)   :math:`\nu`    0.602   0.33     1       Scaling exponent
========== =========== ======== ======== ======== ==================================

Usage
=========================================

.. code-block:: matlab

        info = randcoil()

Returns an ``info`` structure containing the specifics of the model:

* ``info.Model`` -  Full name of the parametric model.
* ``info.Equation`` -  Mathematical equation of the model.
* ``info.nParam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

.. code-block:: matlab

    P = randcoil(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

