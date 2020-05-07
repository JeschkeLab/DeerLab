%============================================================================
% DeerLab Example:
% Multi-Gauss fitting of a 4-pulse DEER signal
%============================================================================

% This example showcases how to fit a simple 4-pulse DEER signal with
% background using a multi-Gauss model, i.e automatically optimizing the
% number of Gaussians in the model.

clear, clc, clf

%==================
% Generate data
%==================

rng(1)
t = linspace(-0.25,4,300); % time axis, us
r = linspace(2.5,4.5,300); % distance axes
param0 = [3 0.3 0.2 3.5 0.3 0.65 4 0.2]; % parameters for three-Gaussian model
P = dd_gauss3(r,param0); % model distance distribution
lam = 0.25;
B = bg_exp(t,0.05);

% Generate dipolar signal with noise
V = dipolarsignal(t,r,P,lam,B,'NoiseLevel',0.01);

%=========================
% Run multi-Gauss fitting
%=========================
NGauss = 5; % maximum number of Gaussians

% Launch a multi-Gauss fit, including the fitting of an exponential background 
[Pfit,param,Pci,paramci,Nopt,metrics,Peval] = fitmultimodel(V,t,r,@dd_gauss,NGauss,...
    'aic','background',@bg_exp,'confidencelevel',0.95,'multistart',1);

% Construct the fitted dipolar kernel
K = dipolarkernel(t,r,param(end-1),bg_exp(t,param(end)));
Vfit = K*Pfit;

% When comparing different parametric models is always a good idea to look
% at the Akaike weights for each model. They basically tell you the
% probability of a model being the most optimal choice.

% Compute the Akaike weights
d = metrics - min(metrics);
w = 100*exp(-d/2)/sum(exp(-d/2));

% Plot results
subplot(321), cla
hold on
plot(t,V,'k.','LineWidth',1)
plot(t,Vfit,'b','LineWidth',1.5)
plot(t,(1-param(end-1))*bg_exp(t,param(end)),'r','LineWidth',1.5)
box on,legend('model','fit')
xlabel('time (\mus)'),ylabel('S(t)')
axis tight

subplot(322), cla
hold on
plot(r,P,'k','LineWidth',1.5)
plot(r,Pfit,'b','LineWidth',1.5)
fill([r fliplr(r)], [Pci(:,1); flipud(Pci(:,2))],'b','Linestyle','none','facealpha',0.25)
box on, axis tight 
legend('model','optimal fit','95%-CI')
xlabel('distance (nm)'), ylabel('P(r)')

subplot(323), cla
hold on
ax = 1:length(metrics);
plot(ax,w,'-o','LineWidth',1.5)
box on,axis tight
ylabel('Akaike Weight [%]')
xlabel('Number of Gaussians in model')

subplot(3,2,[4 6]), cla
hold on
plot(r,Peval + 2*(1:NGauss).','b-','LineWidth',1.5)
box on,axis tight
set(gca,'ytick',2:2:2*NGauss,'yticklabel',1:NGauss)
xlabel('distance (nm)')
ylabel('Number of Gaussians in model')
