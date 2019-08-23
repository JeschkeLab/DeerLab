.. highlight:: matlab
.. _correctphase:


***********************
:mod:`correctphase`
***********************

Phase correction of complex-valued data

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`[S,p,io] = correctphase(C,p,oc)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **C** - Complex-valued signal (N-array)
    *   **p** - Correction phase (scalar)
    *   **oc** - Imaginary offset correction (boolean)
Returns
    *   **S** - Phase-corrected signal (N-array)
    *   **p** - Correction phase (scalar)
    *   **io**  - Imaginary offset (scalar)

Usage
=========================================

.. code-block:: matlab

     S = correctphase(C)

Performs a phase optimization on the complex-valued data ``C`` by minimization of the imaginary component of the data. The phase corrected data ``S`` is returned normalized.

.. code-block:: matlab

     S = correctphase(C,p)

The phase used for correction can be passed manually as a second argument ``p``.

.. code-block:: matlab

    S = correctphase(C,p,oc)

A third boolean argument ``oc`` can be passed to enable/diasable the fitting of a possible offset on the imaginary component of the data (defaults to ``false``).

.. image:: ./images/correctphase1.svg
