.. highlight:: matlab
.. _dd_shellsphere:


************************
:mod:`dd_shellsphere`
************************

Particles distributed on a sphere inside a spherical shell 

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_shellsphere()
        P = dd_shellsphere(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_dd_sphereshell.png
   :width: 25%

:math:`P(r) = \frac{3}{16R_1^3(R_2^3 - R_1^3)}\begin{cases} 12r^3R_1^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_1,R_2 - R_1) \\ 8r^2(R_2^3 - R_1^3) - 3r(R_2^2 - R_1^2)^2 - 6r^3(R_2 - R_1)(R_2 + R_1) \quad \text{for} \quad R_2-R_1 \leq r < 2R_1 \\ 16r^2R_1^3 \quad \text{for} \quad 2R_1\leq r < R_2 - R_1  \\  r^5 - 6r^3(R_2^2 + R_1^2) + 8r^2(R_2^3 + R_1^3) - 3r(R_2^2 - R1_2)^2 \quad \text{for} \quad \max(R_2-R_1,2R_1) \leq r < R_1+R_2 \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

with 

:math:`R_1 = R`

:math:`R_2 = R + w_1`


================ ============== ========= ======== ========= ===================================
 Variable         Symbol         Default   Lower    Upper       Description
================ ============== ========= ======== ========= ===================================
``param(1)``     :math:`R`       1.5       0.1        20         Sphere radius
``param(2)``     :math:`w`       0.5       0.1        20         Shell thickness
================ ============== ========= ======== ========= ===================================


Example using default parameters:

.. image:: ../images/model_dd_shellsphere.png
   :width: 650px


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_shellsphere()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_shellsphere(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

References
=========================================

[1] D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63.
`DOI:  10.1016/j.jmr.2013.01.007 <http://doi.org/10.1016/j.jmr.2013.01.007>`_