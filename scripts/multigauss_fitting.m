%========================================================
% DeerAnalyis2
% Example: Multi-Gauss fitting
% Fit and an unknown number of Gaussian distributions to 
% a noisy signal
%========================================================

clear

% Parameters
%-------------------------
nt = 300;
dt = 0.008; % us
param0 = [3 0.3 4 0.3 0.5]; % parameters for two-Gaussian model

% Preparation
%-------------------------
t = linspace(0,dt*nt,nt);
r = time2dist(t);
P = rd_twogaussian(r,param0);
K = dipolarkernel(t,r);

% Generate dipolar signal with noise
S = dipolarsignal(t,r,P,'NoiseLevel',0.05);

% Run multi-Gauss fitting
%-------------------------
NGauss = 6; % maximum number of Gaussians
[Pfit,param,Nopt,metrics] = multigauss(S,K,r,NGauss,'AICc');

fprintf('The optimal number of Gaussians is: %i \n',Nopt)

% Plot results
%-------------------------
figure(8),clf

subplot(131),hold on
plot(t,S,'k','LineWidth',1)
plot(t,K*Pfit,'b','LineWidth',1.5)
grid on,box on,legend('model','fit')
xlabel('time (\mus)'),ylabel('S(t)')

subplot(132),hold on
plot(r,P,'k','LineWidth',1.5)
plot(r,Pfit,'b','LineWidth',1.5)
grid on, box on, axis tight 
legend('model','fit')
xlabel('distance (nm)'), ylabel('P(r)')

subplot(133),hold on
plot(metrics,'b-o','LineWidth',1.5)
grid on,box on,axis tight
legend('AICc')
ylabel('model selection metric')
xlabel('number of Gaussians in model')
