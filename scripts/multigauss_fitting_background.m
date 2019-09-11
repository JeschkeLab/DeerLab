%========================================================
% DeerAnalyis2
% Example: Multi-Gauss fitting with background
% Fit and an unknown number of Gaussian distributions and
% a background to a noisy signal
%========================================================

clear

% Generate data
%-------------------------
t = -0.050:0.010:3; % time axis, us
r = time2dist(t); % distance axes
param0 = [3 0.3 4 0.3 0.3]; % parameters for two-Gaussian model
P = rd_twogaussian(r,param0); % model distance distribution
lam = 0.25;
B = td_exp(t,0.3);
K = dipolarkernel(t,r,lam,B);

% Generate dipolar signal with noise
V = dipolarsignal(t,r,P,'NoiseLevel',0.05,'ModDepth',lam,'Background',B);

% Run multi-Gauss fitting
%-------------------------
NGauss = 5; % maximum number of Gaussians
[Pfit,param,Nopt,metrics,Peval] = fitmultigauss(V,K,r,NGauss,'AICc');
Vfit = K*Pfit;

fprintf('The optimal number of Gaussians is: %i \n',Nopt)

% Plot results
%-------------------------
figure(8),clf

subplot(5,2,[1 3 5]),hold on
plot(t,V,'k.','LineWidth',1)
plot(t,Vfit,'b','LineWidth',1.5)
grid on,box on,legend('model','fit')
xlabel('time (\mus)'),ylabel('S(t)')
axis tight

subplot(9,2,[8 10]);
plot(t,V-Vfit,'k.','LineWidth',1.5);
grid on, box on, axis tight
xlabel('time (\mus)'), ylabel('residual');

subplot(9,2,[2 4]),hold on
plot(r,P,'k','LineWidth',1.5)
plot(r,Pfit,'b','LineWidth',1.5)
grid on, box on, axis tight 
legend('model','optimal fit')
xlabel('distance (nm)'), ylabel('P(r)')

subplot(5,2,[7 9]),hold on
plot(metrics,'b-o','LineWidth',1.5)
grid on,box on,axis tight
legend('AICc')
ylabel('model selection metric')
xlabel('number of Gaussians in model')

subplot(5,2,[8 10]),hold on
plot(r,Peval + (1:NGauss).','b-','LineWidth',1.5)
grid on,box on,axis tight
set(gca,'ytick',1:NGauss,'yticklabel',1:NGauss)
xlabel('distance (nm)')
ylabel('#Gaussians in model')

