.. highlight:: matlab
.. _bg_exvol:

***********************
:mod:`bg_exvol`
***********************

Hard-shell excluded-volume background parametric model

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = bg_exvol()
        P = bg_exvol(r,param)

Parameters
    *   ``t`` - Time axis (N-array), in microseconds
    *   ``param`` - Model parameters
Returns
    *   ``B`` - Model background (N-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

This implements a hard-shell excluded-volume model, with distance of closest approach ``R`` (first parameter) and excited spin concentration ``c_x`` (second parameter, in micromole/L). To get the excited spin concentration, multiply the total spin concentration with the pump excitation efficiency.

The analytical expression for the decay function as a function of ``R`` and ``c_x`` is complicated, see Kattnig et al, J.Phys.Chem.B 2013, 117, 16542.

========== =================== ========= ============= ============= ================================================
 Variable   Symbol              Default   Lower bound   Upper bound      Description
========== =================== ========= ============= ============= ================================================
param(1)    :math:`R`              1          0.1            20           Exclusion radius
param(1)    :math:`c_x`            50         0.01          1000          Excited spin concentration (uM)
========== =================== ========= ============= ============= ================================================

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = bg_exvol()

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    B = bg_exvol(t,param)

Computes the background model ``B`` from the axis ``t`` according to the parameters array ``param``. The required parameters can also be found in the ``info`` structure.

