Tutorial
=======================================================

Generating datasets
------------------------------

It takes only a few lines to generate synthetic datasets:

.. code-block:: matlab
    

    t = linspace(0,3,200);          % Set up time axis
    r = time2dist(t);               % Set up distance axis
    P0 = rd_onegaussian(r,[3 0.2]); % Generate a one-Gaussian model distribution
    B = td_exp(t,0.2);              % Generate exponential background
    
    % Generate associated dipolar signal
    [V,S] = dipolarsignal(t,r,P0,'Background',B,'ModDepth',0.3,'NoiseLevel',0.02);
    
    subplot(2,1,1)
    plot(r,P0);
    subplot(2,1,2)
    plot(t,S,t,V)

Analyzing datasets
------------------------------

Also, it takes only a few lines to analyze an experimental dataset. Here is how to fit the above generated data using a streched-exponential background model followed by Tikhonov regularization:

.. code-block:: matlab
    
    % assuming data are in V and t
    [B,lam] = fitbackground(V,t,@td_exp);        % Fit background
    K = dipolarkernel(t,r,lam,B);                % Prepare dipolar kernel
    Pfit = fitregmodel(V,K,r,'tikhonov','aic');  % Run Tikhonov regularization
    Vfit = K*Pfit;                               % Calculate fitted time-domain signal
    
    subplot(2,1,1)
    plot(t,V,'.',t,Vfit,t,(1-lam)*B);
    subplot(2,1,2)
    plot(r,P0,r,Pfit);


.. _physicalunits:

Physical units
------------------------------

In DeerAnaysis are distance-domain variables are given in nanometers, whereas all time-domain variables are given in microseconds.