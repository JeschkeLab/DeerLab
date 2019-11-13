.. highlight:: matlab
.. _specunits:


***********************
:mod:`specunits`
***********************

Specification of physical units for time and distance vectors

Syntax
=========================================

.. code-block:: matlab

    [ax,unitout] = specunits(ax,unit)


Parameters
    *   ``ax`` -  Time/Distance axis (array)
    *   ``unit`` -  Unit sufix (string)
Returns
    *  ``ax`` - Converted time/distance axis (array)
    *  ``unitout`` - Unit sufix of converted axis (string)


Description
=========================================

.. code-block:: matlab

            [t,unitout] = specunits(t,'ns')
            [t,unitout] = specunits(t,'us')
            [r,unitout] = specunits(r,'nm')
            [r,unitout] = specunits(r,'A')


Evaluates the input time-axis ``t`` or distance axis ``r`` and adapts the axis accordingly such that they are recognized by all functions as the unit prefix specified in the input. The new units of the output vector are specified in the second output ``unitout`` as a string.


.. important:: The automatic unit determination of DeerAnalysis described in the section :ref:`physicalunits` covers all conditions found in experimental studies. This function is designed to ensure proper functionality of DeerAnalysis if the user desires to conduct studies outside of those regimes. 

For example, the following distance-domain vector in Angstrom

.. code-block:: matlab

            r = linspace(0.5,12,100) %A

would be recognized as nanometers by the unit-heuristics of DeerAnalysis. By specifying the units via

.. code-block:: matlab

            [r,unitout] = specunits(r,'A')

the distance axis is converted to nanometers (as specifed in ``unitout = 'nm'``) such that it is properly processed by DeerAnaysis.