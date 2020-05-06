%============================================================================
% DeerLab Example:
% Confidence intervals for regularization results
%============================================================================

% In this example we will show a simpe example of uncertainty estimation for 
% Tikhonov regularization results. The example will cover the use of confidence intervals
% obtained from curvature matrices and boostrap analysis.

clear,clc,clf 

%===================
% Simulate the data
%===================

% Let's start by generating some data.

%Prepare signal components
t = linspace(-0.4,3.5,200);

%Use a distance-axis with less points to make analysis faster
r = linspace(2,5,200);

P = dd_gauss2(r,[3 0.3 3.5 0.4 0.6]);
B = bg_strexp(t,[0.04,1]);
lam = 0.32;

%Simulate signal
V = dipolarsignal(t,r,P,lam,B,'noiselevel',0.01);

% For the sake of simplicity, in this examples we will assume that we know the 
% background exactly. Our first step is to generate the proper dipolar kernel.

%Generate dipolar kernel with background
KB = dipolarkernel(t,r,lam,B);

%=========================================
% Curvature matrix confidence intervals
%=========================================

% We now have all the elements required to fit our distance distribution via 
% regularization. We will use the AIC to select the regularization parameter in 
% the Tikhonov regularization.

% By default, fitregmodel returns confidence intervals for the fitted
% distribution. (see |help fitregmodel|).

%Fit data via regularization
[Pfit,Pci] = fitregmodel(V,KB,r,'tikh','aic');
%Obtain time-domain fit
Vfit = KB*Pfit;

%Plot the fit results
subplot(311)
plot(t,V,'k.',t,Vfit,'r','LineWidth',1)
grid on, axis tight,box on
xlabel('time [\mus]')
ylabel('V(t)')
legend('Truth','Fit')

subplot(312)
cla,hold on
plot(r,P,'k',r,Pfit,'r','LineWidth',1)
fill([r fliplr(r)],[Pci(:,1); flipud(Pci(:,2))],'r','FaceAlpha',0.2,'LineStyle','none')
grid on, axis tight,box on
xlabel('distance [nm]')
ylabel('P(r) [nm^{-1}]')
title('Curvature Matrix CI')
legend('Truth','Fit','95%-CI')


%=========================================
% Bootstrapped confidence intervals
%=========================================

% Now we are interested in the bootstrap confidence intervals. For this, we
% need to define a boot function e.g. |mybootfcn()| which takes a signal as
% output and returns the outputs of interest (Pfit in our example).

% Launch bootstrapping
Nsamples = 100;
booci = bootan(@(V)mybootfcn(V,KB,r),V,Vfit,Nsamples);
Pci = booci{1}.ci95;

% By plotting the results, one can see that the bootstrapped confidence intervals 
% are narrower in comparison to the ones obtained via the curvature
% matrices. This is due to the inherent accurate nature of bootstrapping. 

subplot(313)
cla,hold on
plot(r,P,'k',r,Pfit,'b','LineWidth',1)
fill([r fliplr(r)],[Pci(:,1); flipud(Pci(:,2))],'b','FaceAlpha',0.2,'LineStyle','none')
grid on, axis tight,box on
xlabel('distance [nm]')
ylabel('P(r) [nm^{-1}]')
title('Bootstrapped CI')
legend('Truth','Fit','95%-CI')

function Pfit = mybootfcn(V,KB,r)
 Pfit = fitregmodel(V,KB,r,'tikh','aic');
end