%============================================================================
% DeerLab Example:
% Global model fits with global, local and fixed parameters
%============================================================================

% This example shows how to fit multiple signals to a global model, which
% may depend on some parameters which need to be globally fitted, some
% locally and some might be fixed and not fitted. 
% This kind of script is the basis for the analysis of pseudo-titration
% curves or dose-response curves.

clc,clear,clf

% Generate two datasets
%-----------------------------------------------------------------------------

% For this example we will simulate a system containing two states A and B
% both havng a Gaussian distribution of known width but unknown mean
% distance. For this system we have two measurements V1 amd V2 measured
% under two different conditions leading to different fractions of states A
% and B. 

r = linspace(2,6,300); % distance axis [nm]
t1 = linspace(0,4,200); % time axis of first measurement [us]
t2 = linspace(0,6,150); % time axis of first measurement [us]

% Parameters
rmeanA = 3.45; % mean distance state A [nm]
rmeanB = 5.05; % mean distance state B [nm]
fwhmA = 0.5; % FWHM state A [nm]
fwhmB = 0.3; % FWHM state B [nm]

fracA1 = 0.8; % Molar fraction of state A under conditions 1
fracA2 = 0.2; % Molar fraction of state A under conditions 2
% The molar fraction of state B is not required as it follows fracB = 1 - fracA

% Generate the two distributions for conditions 1 & 2
P1 = dd_gauss2(r,[rmeanA fwhmA fracA1 rmeanB fwhmB]);
P2 = dd_gauss2(r,[rmeanA fwhmA fracA2 rmeanB fwhmB]);

% ...and the two corresponding signals
V1 = dipolarsignal(t1,r,P1,'noiselevel',0.01);
V2 = dipolarsignal(t2,r,P2,'noiselevel',0.01);
% (for the sake of simplicity no background and 100% modulation depth are assumed)


% Globa fit
%-----------------------------------------------------------------------------

% Now when considering such systems is always important to (1) identify the
% parameters which must be fitted and (2) identify which parameters are the
% same for all signals (global) and which are specific for a individual
% signal (local). 

% In this examples we have the following parameters:
%   - fixed: fwhmA, fwhmB (known paramters)
%   - global: rmeanA, rmeanB (same for both signals)
%   - local: fracA1, fracA2 (different for both signals/conditions)

% The next step is to construct the model function which describes our
% system (see function definition of myABmodel() below)

% Generate the corresponding dipolar kernels
K1 = dipolarkernel(t1,r);
K2 = dipolarkernel(t2,r);

%-----------------------------------------
%                Fit parameters 
%-----------------------------------------
%        [rmeanA rmeanB fracA1 fracA2]
%-----------------------------------------
par0 =   [2        2     0.5    0.5];
lower =  [1        1      0      0];
upper =  [20       20     1      1];
%-----------------------------------------

% ollect data for global fit into cell arrays
Ks = {K1,K2};
Vs = {V1,V2};
ts = {t1,t2};

% Prepare the model function, when using local fit parameters is very 
% important to request the |idx| variable (see below). 
model = @(t,par,idx) myABmodel(t,par,idx,r,Ks);

% Fit the global parametric model to both signals
parfit = fitparamodel(Vs,model,ts,par0,'lower',lower,'upper',upper,'multistart',50);

% The use of the option 'multistart' will help the solver to find the
% global minimum and not to get stuck at local minima.

% Get the fitted models 
[Vfit1,Pfit1] = myABmodel(t1,parfit,1,r,Ks);
[Vfit2,Pfit2] = myABmodel(t2,parfit,2,r,Ks);

% Plot results
%-----------------------------------------------------------------------------
subplot(221)
plot(r,P1,'k',r,Pfit1,'r')
axis tight,grid on
xlabel('r [nm]'),ylabel('P(r) [nm]^{-1}')
title('Conditions #1')
legend('truth','fit')

subplot(222)
plot(t1,V1,'k.',t1,Vfit1,'r')
xlabel('t [\mus]'),ylabel('V(t)')

subplot(223)
plot(r,P2,'k',r,Pfit2,'b')
axis tight,grid on
xlabel('r [nm]'),ylabel('P(r) [nm]^{-1}')
title('Conditions #2')
legend('truth','fit')

subplot(224)
plot(t2,V2,'k.',t2,Vfit2,'b')
xlabel('t [\mus]'),ylabel('V(t)')

% Model definition
%-----------------------------------------------------------------------------
function [Vfit,Pfit] = myABmodel(~,par,idx,r,Ks)

    % This function models the signals in our A-B system, and it is used to
    % simulate all signals passed to fitparamodel(). 
    % To be able to distinguish which signal is being analyzed at the
    % moment by fitparamodel(), the function will provide the variable
    % |idx| which identifys the signal, i.e.
    %   idx = 1 ---> Vs{1}
    %   idx = 2 ---> Vs{2}

    %Fixed parameters
    fwhmA = 0.5;
    fwhmB = 0.3;
    %Global parameters
    rmeanA = par(1);
    rmeanB = par(2);
    %Local parameters
    fracA1 = par(3);
    fracA2 = par(4);

    % Depending on which signal is being simulated at the moment, use the
    % corresponding fraction
    if idx == 1
        %Use fracA1 for signal Vs{1}
        fraction = fracA1;
    else
        %Use fracA2 for signal Vs{2}
        fraction = fracA2;
    end
    
    % Generate the signal-specific distribution
    Pfit = dd_gauss2(r,[rmeanA fwhmA fraction rmeanB fwhmB]);
    
    % Generate the corresponding signal (either V{1} or V{2}) by using the
    % appropiate dipolar kernel
    Vfit = Ks{idx}*Pfit;

end