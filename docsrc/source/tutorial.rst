Tutorial
=======================================================

Generating datasets
------------------------------

It takes only a few lines to generate synthetic datasets:

.. code-block:: matlab
    
    t = linspace(0,3,200);          % Set up time axis
    r = time2dist(t);               % Set up distance axix
    P0 = rd_onegaussian(r,[3 0.2]); % Generate a one-Gaussian model distribution
    B = td_strexp(t,[0.2 3]);       % Generate stretched-exponential background

    % Generate associated dipolar signal
    [V,S] = dipolarsignal(t,r,P0,'Background',B,'ModDepth',0.3,'Noiselevel',0.02);

    subplot(2,1,1)
    plot(r,P0);
    subplot(2,1,2)
    plot(t,S,t,V)

Analyzing datasets
------------------------------

Also, it takes only a few lines to analyze an experimental dataset. Here is how to fit the above generated data using a streched-exponential background model followed by Tikhonov regularization:

.. code-block:: matlab
    
    % assuming data are in V and t
    [B,lam] = fitbackground(V,t,@td_strexp,1.5);   % Fit background
    K = dipolarkernel(t,r,B,lam);                  % Prepare dipolar kernel
    L = regoperator(numel(r),2);                   % Prepare 2nd order regularization operator
    Pfit = fitregmodel(S,K,r,L,'tikhonov',1);      % Run Tikhonov regularization
    
    subplot(2,1,1)
    plot(t,V,t,(1-lam)*B,t,F);
    subplot(2,1,2)
    plot(r,Pfit,r,P0);

