.. highlight:: matlab
.. _ex_4pdeer:


***********************
:mod:`ex_4pdeer`
***********************

4-pulse DEER experiment 

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = ex_4pdeer(t)
        pathways = ex_4pdeer(t,param)

Parameters
    *   ``t`` - Time axis (*M*-array)
    *   ``param`` - Model parameters (array)
Returns
    *   ``pathways`` - Dipolar pathways (matrix)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_ex_4pdeer.png
   :width: 550px


:math:`V(t) = [1-\lambda + \lambda D(t-T_0^{(1)})]B(t-T_0^{(1)}) = [1-\lambda + \lambda D(t)]B(t)`

:math:`K(t,r) = [1-\lambda + \lambda K(t-T_0^{(1)},r)]B(t-T_0^{(1)}) = [1-\lambda + \lambda K(t,r)]B(t)`

where :math:`T_0^{(1)}=0` is the refocusing time of the modulated dipolar pathway.


============== ================ ============ ============ ============ ================================================
 Variable        Symbol           Default       Lower        Upper                Description
============== ================ ============ ============ ============ ================================================
``param(1)``   :math:`\lambda`     0.3           0            1          Modulated pathway amplitude (modulation depth)
============== ================ ============ ============ ============ ================================================


Example of a simulated signal using default parameters:

.. image:: ../images/model_ex_4pdeer.png
   :width: 550px

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = ex_4pdeer(t)

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    pathways = ex_4pdeer(t,param)

Generates the dipolar pathways matrix ``pathways`` from the time-axis ``t`` and model parameters ``param``. 


