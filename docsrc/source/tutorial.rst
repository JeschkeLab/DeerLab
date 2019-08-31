Tutorial
=======================================================

Generating datasets
------------------------------

It takes only a few lines to generate synthetic datasets:

.. code-block:: matlab
    
    t = linspace(0,3,200);          % set up time axis
    r = time2dist(t);               % set up distance axix
    P0 = rd_onegaussian(r,[3 0.2]); % generate a one-Gaussian model distribution
    S = dipolarsignal(t,r,P0);      % generate associated dipolar signal
    B = td_strexp(t,[0.2 3]);       % generate stretched-exponential background
    lam = 0.3;
    V0 = ((1-lam)+lam*S).*B;        % generate complete dipolar signal
    V = V0 + randn(size(V0))*0.02;  % add noise
    
    subplot(2,1,1)
    plot(r,P0);
    subplot(2,1,2)
    plot(t,S,t,V)

Analyzing datasets
------------------------------

Also, it takes only a few lines to analyze an experimental dataset. Here is how to fit the above generated data using a streched-exponential background model followed by Tikhonov regularization:

.. code-block:: matlab
    
    % assuming data are in V and t
    [B,lam] = fitbackground(V,t,@td_strexp,1.5);   % fit background
    S = (V./B - (1-lam))/lam;                      % calculate form factor
    K = dipolarkernel(t,r);                        % prepare dipolar kernel
    L = regoperator(numel(r),2);                   % prepare regularization operator
    Pfit = fitregmodel(S,K,r,L,'tikhonov',1);      % run Tikhonov regularization
    
    subplot(2,1,1)
    plot(t,V,t,(1-lam)*B,t,F);
    subplot(2,1,2)
    plot(r,Pfit,r,P0);

