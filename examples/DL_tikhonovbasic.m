%============================================================================
% DeerLab Example:
% Basic fitting of a 4-pulse DEER signal, parameter-free distribution
%============================================================================

% This example shows how to fit a simple 4-pulse DEER signal with a parameter-
% free distribution, a background, and a modulation amplitude.

clear, clc, clf

% Generate data
%-----------------------------------------------------------------------------
rng(1)
t = -0.1:0.02:4;                 % time axis, us
r = linspace(1.5,6,numel(t));    % distance axis, ns
param0 = [3 0.3 0.2 3.5 0.3 0.65 3.8 0.2]; % parameters for three-Gaussian model
P = dd_gauss3(r,param0);         % model distance distribution
lam = 0.5;                       % modulation depth
B = bg_hom3d(t,300,lam);         % background decay
Vexp = dipolarsignal(t,r,P,lam,B,'NoiseLevel',0.01);

% Run fit
%-----------------------------------------------------------------------------
[Vfit,Pfit,Bfit,parfit,parci] = ...
    fitsignal(Vexp,t,r,'P',@bg_hom3dex,@ex_4pdeer);


% Plot results
%-----------------------------------------------------------------------------
subplot(211)
plot(t,Vexp,'.',t,Vfit)
legend('data','fit')

subplot(212)
plot(r,P,r,Pfit);
legend('model','fit')
