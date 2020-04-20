.. _dipolarbackground:

*************************
:mod:`dipolarbackground`
*************************

Constructs the total multi-pathway background function

-------------------------------


Syntax
=========================================

.. code-block:: matlab

    B = dipolarbackground(t,pathinfo,Bmodel)
    B = dipolarbackground(___,'Property',Value)


Parameters
    *   ``t``        - Time axis vector (*N*-element array), in microseconds
    *   ``pathinfo`` - Modulation depths, refocusing times, and harmonics (*px2* or *px3* array) for multiple dipolar pathways
    *   ``Bmodel``        - Background model function handle
Returns
    *  ``B`` - Multi-pathway background (*N*-element array)

-------------------------------


Description
=========================================

.. code-block:: matlab

   B = dipolarbackground(t,pathinfo,Bmodel)

Computes the multipathway background ``B`` for the time axis ``t`` (in microseconds) corresponding to the dipolar pathways specified in ``pathinfo``. The total background is computed from the basis background function model specified in ``Bmodel``. This model must be a function handle of the type ``@(t) bg_model(t,params)``.

For a multi-pathway DEER signal (e.g, 4-pulse DEER with 2+1 contribution; 5-pulse DEER with 4-pulse DEER residual signal, and more complicated experiments), ``pathinfo`` is a 2D array that contains a list of modulation depths (amplitudes), refocusing times (in microseconds), and optional harmonics for all modulated pathway signals.

Each row of ``pathinfo`` contains two values: one modulation depth and one refocusing time. For a pathway with unmodulated signal, set the refocusing time to ``NaN``.

Optionally, the harmonic (1 = fundamental, 2 = first overtone, etc.) can be given as a third value in each row. This can be useful for modeling RIDME signals. If not given, the harmonic is 1 for all pathways.

Example:
	To specify the standard model for 4-pulse DEER with an unmodulated offset and a single dipolar pathway that refocuses at time 0, use

.. code-block:: matlab

    lambda = 0.4; % modulation depth main signal
    kappa = 0.3;
    
    pathinfo = [1-lambda NaN; lambda 0];
    
    % alternative input
    pathinfo(1,:) = [1-lambda NaN];    % unmodulated part, gives offset
    pathinfo(2,:) = [lambda 0];        % main modulation, refocusing at time zero
    
    Bmodel = @(t) bg_exp(t,kappa)
    
    B = dipolarbackground(t,pathinfo,Bmodel)



-------------------------------



Additional Settings
=========================================


Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    B = dipolarbackground(___,'Property1',Value1,'Property2',Value2,___)


- ``'OvertoneCoeffs'`` - RIDME overtone coefficients
    1D array containing the overtone coefficients for RIDME experiments. 
    
    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			B = dipolarbackground(___,'OvertoneCoeffs',[0.4 0.2 0.4])   % fundamental, 1st, and 2nd overtone


- ``'Renormalize'`` - Re-normalization of multi-pathway background
    The multi-pathway background does not necessarily satisfy ``V(0) == 1``. This option enables(``true``) or disables(``false``) a re-normalization to ensure that equality is satisfied.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

        B = dipolarbackground(___,'Renormalize',false)
