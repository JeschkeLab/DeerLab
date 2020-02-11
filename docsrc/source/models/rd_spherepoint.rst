.. highlight:: matlab
.. _rd_spherepoint:


************************
:mod:`rd_spherepoint`
************************

One particle distanced from particles distributed on a sphere


-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = rd_spherepoint()
        P = rd_spherepoint(r,param)

Parameters
    *   ``r`` - Distance axis (N-array)
    *   ``param`` - Model parameters
Returns
    *   ``P`` - Distance distribution (N-array)
    *   ``info`` - Model information (struct)

-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_rd_spherepoint.png
   :width: 25%

:math:`P(r) = \begin{cases} \frac{3r(R^2-(d-r)^2)}{4dR^3} \quad \text{for} \quad d-R \leq r < d+R \\ 0 \quad \text{for} \quad \text{otherwise}  \end{cases}`


================ ============== ========= ======== ========= ===================================
 Variable         Symbol         Default   Lower    Upper       Description
================ ============== ========= ======== ========= ===================================
``param(1)``     :math:`R`       1.5       0.1        20        Sphere radius
``param(2)``     :math:`d`       3.5       0.1        20        Distance to point
================ ============== ========= ======== ========= ===================================


Example using default parameters:

.. image:: ../images/model_rd_spherepoint.png
   :width: 40%


-----------------------------


Description
=========================================

.. code-block:: matlab

        info = rd_spherepoint()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    P = rd_spherepoint(r,param)

Computes the distance distribution model ``P`` from the axis ``r`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

References
=========================================

[1] D.R. Kattnig, D. Hinderberger, Journal of Magnetic Resonance, 230 (2013), 50-63.
`DOI:  10.1016/j.jmr.2013.01.007 <http://doi.org/10.1016/j.jmr.2013.01.007>`_