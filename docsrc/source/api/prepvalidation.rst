.. highlight:: matlab
.. _prepvalidation:

***********************
:mod:`prepvalidation`
***********************

Preparation of the validation parameter vectors

Syntax
=========================================

.. code-block:: matlab

    varparam = prepvalidation(param)
    varparam = prepvalidation(param,'Property',Value)

Parameters
    *   ``param`` - Validation parameter settings (struct array)

Returns
    *   ``varparam`` - Cell array of parameter combinations

Description
=========================================

.. code-block:: matlab

    varparam = prepvalidation(param)

Returns all the possible permutations of the parameters in the input structure array ``param``. All possible combinations are randomly permuted and returned as a cell array. The ``param`` structure must have to following fields for all of the N parameters:

*   ``param(n).name`` - Name of the parameter in the script (string)
*   ``param(n).values`` - Values to be adapted by the parameter (Cell or numerical array)

.. Important:: This function is automatically called by :ref:`validate`, where the ``param`` structure and not ``varparam`` must be passed.

For example, a case where two parameters ``varA`` and ``varB`` are to be validated. The first parameter should be validated times at for the strings ``'str1'`` and ``'str2'``  and the second parameter should be linearly validated 3 times in the range 5 to 15:

.. code-block:: matlab

    %Validation parameter varA
    param(1).name = 'varA';
    param(2).values = {'str1','str2'};

    %Validation parameter varB
    param(2).name = 'varB';
    param(2).values = linspace(5,15,3);

    varparam = prepvalidation(param);

which will compute all the possible combinations of the values of the first and the second parameter. The cell array ``varparam`` will contain the following permutations:

================== ============= ===========
    Combination #   ``varA``      ``varB``
================== ============= ===========
        1           ``'str1'``       5
        2           ``'str2'``       5
        3           ``'str1'``       10
        4           ``'str2'``       10
        5           ``'str1'``       15
        6           ``'str2'``       15
================== ============= ===========


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    varparam = prepvalidation(param,'Property1',Value1,'Property2',Value2)

RandPerm
    Specifies whether to randomly permute the validation parameters combinations.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

         varparam = prepvalidation(param,'RandPerm',false)

