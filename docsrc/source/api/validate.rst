.. highlight:: matlab
.. _validate:

***********************
:mod:`validate`
***********************

Statistical validation of script/function parameters.

Syntax
=========================================

.. code-block:: matlab

    [mean,std] = validate(op,param)
    [mean,std] = validate(op,param,file)
    [mean,std] = validate(op,param,file,'Property',Value)

Parameters
    *   ``op`` - Name of output parameter (string)
    *   ``param`` - Validation parameter settings (struct array)
    *   ``file`` - External file name (string)

Returns
    *   ``mean`` - Validated ``op`` mean value
    *   ``std`` - Validated ``op`` standard deviation

Description
=========================================

.. code-block:: matlab

    [mean,std] = validate(op,param)

Validates the objective parameter with name ``op`` in the script/function from which ``validate()`` is called. The validation settings are defined by the ``param`` structure array:

*   ``param(n).name`` - Name of the parameter in the script (string)
*   ``param(n).values`` - Values to be adapted by the parameter (Cell or numerical array)

The objective parameter ``op`` is then evaluated for all possible combinations of the parameter values. The mean and standard deviation of the different computed objective parameter values are returned as the output arguments ``mean`` and ``std``, respectively.

.. code-block:: matlab

    [mean,std] = validate(op,param,file)

If the code to be evaluated is contained in a different file than the `validate()` call, the name ``filename`` of said file can be specified as a third argument. The parameter names ``op`` and ``param(n).name`` must correspond to variable names in that file.

Optional Arguments
=========================================

Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    [mean,std] = validate(op,param,'Property1',Value1,'Property2',Value2)
    [mean,std] = validate(op,param,file,'Property1',Value1,'Property2',Value2)

RandPerm
    Specifies whether to randomly permute the validation parameters combinations.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

        [mean,std] = validate(op,param,'RandPerm',false)

AxisHandle
    Axis handle to plot the state of the validation results at each parameter combination.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        [mean,std] = validate(op,param,'AxisHandle',gca)

Example
=========================================

.. code-block:: matlab


    clc,clf
    %Parameters
    N = 200;
    regparam = 5;
    Lorder = 2;
    validationnoise = 0;

    %Preparation
    t = linspace(0,4,N);
    r = time2dist(t);
    P = rd_onegaussian(r,[4,0.3]);
    K = dipolarkernel(t,r);
    S = K*P;
    S = dipolarsignal(t,r,P,'noiselevel',0.05);
    L = regoperator(N,Lorder);

    %Add extra noise to validate its effects
    S = S + whitenoise(M,validationnoise);

    %Use Tikhonov regularization
    Pfit = fitregmodel(S,K,r,L,'tikh',regparam);

    %Define validation parameters
    ValParam(1).name = 'regparam'
    ValParam(1).value = logspace(-2,2,25)

    ValParam(2).name = 'validationnoise'
    ValParam(2).value = linspace(0.01,0.1,10)

    ValParam(3).name = 'Lorder'
    ValParam(3).value = [1 2];

    %Run the validation using the code above to calculate Pfit
    [Pmean,Pstd] = validate('Pfit',ValParam)
