.. highlight:: matlab
.. _dd_randcoil:

***********************
:mod:`dd_randcoil`
***********************

Random-coil model for an unfolded peptide/protein


-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_randcoil()
        P = dd_randcoil(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_dd_randcoil.png
   :width: 35%

:math:`P(r) = \frac{3}{(2\pi\nu_0)^{3/2}}4\pi r^2\exp(-\frac{3 r^2}{\nu_0})`

where :math:`\nu_0 = 3/(12\pi r_0 N \nu)^{3/2}`

============== =========== ======== ======== ======== ==================================
 Variable       Symbol     Default   Lower   Upper       Description
============== =========== ======== ======== ======== ==================================
``param(1)``   :math:`N`      50      2        1000    Number of residues
``param(2)``   :math:`R_0`    0.20    0.10     0.40    Segment length
``param(3)``   :math:`\nu`    0.60    0.33     1.00    Scaling exponent
============== =========== ======== ======== ======== ==================================

Example using default parameters:

.. image:: ../images/model_dd_randcoil.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_randcoil()

Returns an ``info`` structure containing the information of the model parameters and boundaries.

* ``info(n).Index`` -  Index of the parameter in the ``param`` array.
* ``info(n).Parameter`` -  Description of the n-th parameter.
* ``info(n).Lower`` -  Lower bound of the n-th parameter.
* ``info(n).Upper`` -  Upper bound of the n-th parameter.
* ``info(n).Start`` -  Start value of the n-th parameter.

-----------------------------


.. code-block:: matlab

    P = dd_randcoil(r,param)

Computes the model distance distribution ``P`` of residue-to-residue distances ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.


References
=========================================

[1] N. C. Fitzkee, G. D. Rose, PNAS 2004, 101(34), 12497-12502
DOI: `10.1073/pnas.0404236101 <https://doi.org/10.1073/pnas.0404236101>`_