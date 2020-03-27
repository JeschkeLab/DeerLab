.. highlight:: matlab
.. _dd_shellshell:


************************
:mod:`dd_shellshell`
************************

Particles distributed on a spherical shell inside another spherical shell

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = dd_shellshell()
        P = dd_shellshell(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_dd_shellshell.png
   :width: 25%

:math:`P(r) = (R_1^3(R_2^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_2) - R_1^3(R_3^3 - R_1^3)P_\mathrm{BS}(r|R_1,R_3) - R_2^3(R_3^3 - R_2^3)P_\mathrm{BS}(r|R_2,R_3))/((R_3^3 - R_2^3)(R_2^3 - R_1^3))`

with 

:math:`P_\mathrm{BS}(r|R_i,R_j) = \frac{3}{16R_i^3(R_j^3 - R_i^3)}\begin{cases} 12r^3R_i^2 - r^5  \quad \text{for} \quad 0\leq r < \min(2R_i,R_j - R_i) \\ 8r^2(R_j^3 - R_i^3) - 3r(R_j^2 - R_i^2)^2 - 6r^3(R_j - R_i)(R_j + R_i) \quad \text{for} \quad R_j-R_i \leq r < 2R_i \\ 16r^2R_i^3 \quad \text{for} \quad 2R_i\leq r < R_j - R_i  \\  r^5 - 6r^3(R_j^2 + R_i^2) + 8r^2(R_j^3 + R_i^3) - 3r(R_j^2 - R1_2)^2 \quad \text{for} \quad \max(R_j-R_i,2R_i) \leq r < R_i+R_j \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`

and

:math:`R_1 = R`

:math:`R_2 = R + w_1`

:math:`R_3 = R + w_1 + w_2`



================ ============== ========= ======== ========= ===================================
 Variable         Symbol         Default   Lower    Upper       Description
================ ============== ========= ======== ========= ===================================
``param(1)``     :math:`R`       1.5       0.1        20         Inner shell radius
``param(2)``     :math:`w_1`     0.5       0.1        20         1st Shell thickness
``param(3)``     :math:`w_2`     0.5       0.1        20         2nd Shell thickness
================ ============== ========= ======== ========= ===================================


Example using default parameters:

.. image:: ../images/model_dd_shellshell.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = dd_shellshell()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = dd_shellshell(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

References
=========================================

[1] D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63.
`DOI:  10.1016/j.jmr.2013.01.007 <http://doi.org/10.1016/j.jmr.2013.01.007>`_