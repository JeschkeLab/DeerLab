%============================================================================
% DeerLab Example:
% Multi-Gauss fitting of a 4-pulse DEER signal
%============================================================================

clear,clc

rng(1)

% Generate data
%-------------------------
t = linspace(-0.25,4,300); % time axis, us
r = linspace(2.5,4.5,300); % distance axes
param0 = [3 0.3 3.5 0.3 4 0.2 0.1 0.65]; % parameters for two-Gaussian model
P = dd_gauss3(r,param0); % model distance distribution
lam = 0.35;
B = bg_exp(t,0.08);

[r,P] = groundtruth('figure5');


% Generate dipolar signal with noise
V = dipolarsignal(t,r,P,lam,B,'NoiseLevel',0.01);

% Run multi-Gauss fitting
%-------------------------
NGauss = 5; % maximum number of Gaussians
[Pfit,param,Pci,paramci,Nopt,metrics,Peval] = fitmultigauss(V,t,r,NGauss,'aic','background',@bg_exp,'confidencelevel',0.95,'multistart',1);
K = dipolarkernel(t,r,param(end-1),bg_exp(t,param(end)));

Vfit = K*Pfit;

fprintf('The optimal number of Gaussians is: %i \n',Nopt)

% Plot results
%-------------------------
figure(8),clf

figure(1),clf
hold on
plot(t,V,'k.','LineWidth',1)
plot(t,Vfit,'b','LineWidth',1.5)
plot(t,(1-param(end-1))*bg_exp(t,param(end)),'r','LineWidth',1.5)
box on,legend('model','fit')
xlabel('time (\mus)'),ylabel('S(t)')
axis tight
set(gca,'fontsize',14)

figure(2),clf
hold on
plot(r,P,'k','LineWidth',1.5)
plot(r,Pfit,'b','LineWidth',1.5)
fill([r fliplr(r)], [Pci(:,1); flipud(Pci(:,2))],'b','Linestyle','none','facealpha',0.25)
box on, axis tight 
legend('model','optimal fit')
xlabel('distance (nm)'), ylabel('P(r)')
set(gca,'fontsize',14)
% xlim([2.5 4.5])

figure(3),clf
hold on
ax = 1:length(metrics);
w = exp(-(metrics - min(metrics))/2)/sum(exp(-(metrics - min(metrics))/2))
plot(ax,metrics - min(metrics),'-o','LineWidth',1.5)
box on,axis tight
ylabel('AIC')
xlabel('number of Gaussians in model')
set(gca,'fontsize',14)

figure(4),clf
hold on
plot(r,Peval + 2*(1:NGauss).','b-','LineWidth',1.5)
box on,axis tight
set(gca,'ytick',2:2:2*NGauss,'yticklabel',1:NGauss)
xlabel('distance (nm)')
ylabel('#Gaussians in model')
set(gca,'fontsize',14)
% xlim([2.5 4.5])
ylim(2*[0.8 7.5])
