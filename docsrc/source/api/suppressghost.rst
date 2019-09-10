.. highlight:: matlab
.. _suppressghost:

*********************
:mod:`suppressghost`
*********************

Ghost distance suppression in multi-spin systems

Syntax
=========================================

.. code-block:: matlab

   Vs = suppressghost(V,n)

Parameters
    *   ``V`` - Signal vector (N-array)
    *   ``n`` - number of radicals (scalar)
Returns
    *   ``Vs`` Power-scaled signal (N-array)

Description
=========================================

.. code-block:: matlab

   Vs = suppressghost(V,n)

Suppresses multi-spin contributions to the signal ``V`` by means of the power scaling approximation [1]_. The scaling is determined by the number of radicals ``n`` in the system. The function returns the power-scaled signal in ``Vs``.

References
=========================================

.. [1] von Hagens et al., Phys.Chem.Chem.Phys. 15, 2013, 5854-5866, https://doi.org/10.1039/C3CP44462G
