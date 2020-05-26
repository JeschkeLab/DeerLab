.. highlight:: matlab
.. _aptkernel:

*********************
:mod:`aptkernel`
*********************

Computes the dipolar interaction kernel and components required for the approximate Pake transformation (APT).

.. warning:: This is a legacy function. Its use is not recommended for routine or accurate data analysis.


-----------------------------


Syntax
=========================================

.. code-block:: matlab

    K = aptkernel(t)
    K = aptkernel(t,'Property',Value)

Parameters
    *   ``t`` - Time axis (*N*-element array)
Returns
    *   ``K`` - APT kernel elements (struct)

-----------------------------



Description
=========================================

.. code-block:: matlab

    K = aptkernel(t)

Generates a structure a APT kernel structure ``K`` containing the following fields:

    *   ``.Base`` - Time/Frequency dipolar kernel (*(N/2-2)xN*-element matrix) 
    *   ``.NormalizationFactor`` -  Normalization factors (*(N/2-2)*-element array)
    *   ``.FreqAxis`` - Frequency Axis (*(N/2-2)*-element array)
    *   ``.TimeAxis`` -  Time Axis (*N*-element array)
    *   ``.Crosstalk`` -  Crosstalk matrix *(N/2-2)x(N/2-2)*-element matrix)

This structure can be then be passed directly to the :ref:`apt` function for computing the APT. 

Since the APT is based on a time-domain to frequency-domain transformation, the APT kernel and the dipolar kernel obtained from the :ref:`dipolarkernel` function are not interchangeable.



-----------------------------


Additional Settings
=========================================
Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    K = aptkernel(t,'Property1',Value1,'Property2',Value2,___)


-``'ExcitationBandwidth'`` - Excitation bandwidth of the pulses in **MHz**. 
    If specified, its value is used in the compensation of limited excitation bandwidth of the experimental pulses. If not specified infinite excitation bandwidth is assumed.

    *Default:* empty, corresponding to infinite excitation bandwidth

    *Example:*

		.. code-block:: matlab

			K = aptkernel(t,'ExcitationBandwidth',100)   % 100 MHz excitation bandwidth

