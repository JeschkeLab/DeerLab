.. highlight:: matlab
.. _exp_4pdeer:


***********************
:mod:`exp_4pdeer`
***********************

4-pulse DEER experiment 

-----------------------------


Syntax
=========================================

.. code-block:: matlab

        info = exp_4pdeer(t)
        K = exp_4pdeer(t,r,param,Bmodel)
        [K,B] = exp_4pdeer(t,r,param,Bmodel)

Parameters
    *   ``t`` - Time axis (*M*-array)
    *   ``r`` - Distance axis (*N*-array)
    *   ``param`` - Model parameters (array)
    *   ``Bmodel`` - Background model (function handle)
Returns
    *   ``K`` - Dipolar kernel (*MxN*-array)
    *   ``B`` - Experiment background (*M*-array)
    *   ``info`` - Model information (struct)


-----------------------------

Model
=========================================

.. image:: ../images/model_scheme_exp_4pdeer.png
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

.. image:: ../images/model_exp_4pdeer.png
   :width: 550px

-----------------------------


Description
=========================================

.. code-block:: matlab

        info = exp_4pdeer(t)

Returns an ``info`` structure containing the specifics of the model:

* ``info.model`` -  Full name of the parametric model.
* ``info.nparam`` -  Total number of adjustable parameters.
* ``info.parameters`` - Structure array with information on individual parameters.

-----------------------------


.. code-block:: matlab

    [K,B] = exp_4pdeer(t,r,param,Bmodel)

Computes the distance distribution model ``P`` from the time axis ``t`` and distance axis ``r`` according to the parameters array ``param``.  The required parameters can also be found in the ``info`` structure. 

The full background ``B`` is also computed from the basic background model ``Bmodel``, which has to be passed as a function of the time-axis ``t``. For example: 

.. code-block:: matlab

    Bmodel = @(t) bg_exp(t,k);
    [K,B] = exp_4pdeer(t,r,param,Bmodel)


