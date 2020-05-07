%============================================================================
% DeerLab Example:
% Global fit of dipolar evolution functions using fitregmodel
%============================================================================

% This short example will introduce the use of global fitting (or global analysis) 
% using Tikhonov regularization via the fitregmodel() function.

clear,clc,clf

% Generating multiple datasets
%-----------------------------------------------------------------------------

% First, let's generate the three dipolar evolution functions D1,D2 and 
% D3 (i.e. without background and full modulation depth), but different length, 
% and noise levels. All of them arising from the same bimodal distance distribution. 
% We assume the modulation depths to be approximately equal between measurements.

%Model input parameters
k = 0.3;
lam = 0.35;
rmean1 = 3;
rmean2 = 3.5;
w1 = 0.3;
w2 = 0.3;
A1 = 0.3;

%Construct the axes
t1 = linspace(0,5,200);
t2 = linspace(0,2,100);
t3 = linspace(0,3,150);

%Generate distance distribution
r = linspace(2,5,100);
P = dd_gauss2(r,[rmean1 w1 A1 rmean2 w2]);

%Generate dipolar signals
D1 = dipolarsignal(t1,r,P,'noiselevel',0.06);
D2 = dipolarsignal(t2,r,P,'noiselevel',0.01);
D3 = dipolarsignal(t3,r,P,'noiselevel',0.04);

% Global fit
%-----------------------------------------------------------------------------

% Doing this global fit does not require any additional steps in comparison 
% to the local fits. In order to perform the fit in a global way, we only need 
% to pass all signals as a cell array |{D1,D2,D3}| as well as the dipolar kernels 
% for each signal.

%Construct the dipolar kernels
K1 = dipolarkernel(t1,r);
K2 = dipolarkernel(t2,r);
K3 = dipolarkernel(t3,r);

%Run global fit of time-domain parameteric model
Pfit = fitregmodel({D1,D2,D3},{K1,K2,K3},r,'tikh','aic');

%Get the time-domain fits of all signals
D1fit = K1*Pfit;
D2fit = K2*Pfit;
D3fit = K3*Pfit;

% Finally we can plot all our results and compare them to the ground truth.
subplot(221)
plot(t1,D1,'k.',t1,D1fit,'b',t2,D2+1/2,'k.',t2,D2fit+1/2,'b',t3,D3+1,'k.',t3,D3fit+1,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')
set(gca,'FontSize',14)

subplot(222)
plot(r,P,'k',r,Pfit,'b','LineWidth',1.5)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('P(r) [nm^{-1}]')
set(gca,'FontSize',14)

% Studying the influence of the global weights
%-----------------------------------------------------------------------------

% By default, the global weights, i.e. how much influence each signal has on 
% the global fit, are chosen automatically based on the individual signal lengths 
% and approximate noise levels. However, this weights can be specified by passing 
% the option |'GlobalWeights|'. In this example, the signal |D1| has the largest 
% amount of points and the highest noise. Hence by default this signals will be 
% less weighted than the others. In this last section we want to study the effect 
% of setting all weights equal for all signals has on the gobal fit.

%Set all global weigths equal (w=1/3)
weights = ones(3,1);

%Run global fit of time-domain parameteric model
Pfit = fitregmodel({D1,D2,D3},{K1,K2,K3},r,'tikh','aic','GlobalWeights',weights);

%Get the time-domain fits of all signals
D1fit = K1*Pfit;
D2fit = K2*Pfit;
D3fit = K3*Pfit;

%Plotting
subplot(223)
plot(t1,D1,'k.',t1,D1fit,'r',t2,D2+1/2,'k.',t2,D2fit+1/2,'r',t3,D3+1,'k.',t3,D3fit+1,'r','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')
set(gca,'FontSize',14)

subplot(224)
plot(r,P,'k',r,Pfit,'r','LineWidth',1.5)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('P(r) [nm^{-1}]')
set(gca,'FontSize',14)