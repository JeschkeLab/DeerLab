%========================================================
% DeerAnalyis2
% Example: Time-domain parametric model fitting
% Construct the a time-domain model of the DEER signal 
% including background and fit distribution+background at 
% the same time.
%========================================================

clear, clc

% Model parameters
%----------------------------------------------
r1 = 6; w1 = 0.3; % center and width of first Gaussian, nm
r2 = 4; w2 = 0.3; % center and width of second Gaussian, nm
amp1 = 0.3; % amplitude of first Gaussian
lam = 0.3; % modulation amplitude
k = 0.3; % beckground decay constant
sigma = 0.01; % noise level

% Generate signal
%----------------------------------------------
t = linspace(0,5,251);
r = time2dist(t);
K = dipolarkernel(t,r);
P = rd_twogaussian(r,[r1 w1 r2 w2 amp1]);
B = td_exp(t,k);
V = dipolarsignal(t,r,P,'ModDepth',lam,'Background',B,'Noiselevel',sigma);

% Define model
%----------------------------------------------
% Construct time-domain model function including background
mymodel = @(t,p) td_exp(t,p(2)).*((1- p(1)) + p(1)*K*rd_twogaussian(r,p(3:end)));

% Define lower/upper bounds and initial guess of parameters
upper = [1 200 20 5 20 5 1];
lower = [0 0 1.0 0.05 1.0 0.05 0];
param0 = [0.5 0.35 6 0.2 3.5 0.4 0.4];

% Convert model function to valid parametric model
model = paramodel(mymodel,param0,lower,upper);

% Fit the model to time-domain signal
%----------------------------------------------
[param,Vfit] = fitparamodel(V,mymodel,t,param0,'Upper',upper,'Lower',lower);
Pfit = rd_twogaussian(r,param(3:end));

% Plotting
%----------------------------------------------
figure(1),clf

subplot(2,1,1)
plot(t,V,'.',t,Vfit,'LineWidth',1.5)
xlabel('time (\mus)')
ylabel('V(t)')
grid on,axis tight, box on
legend('Data','Fit')

subplot(2,2,3)
plot(r,P,r,Pfit,'LineWidth',1.5);
xlabel('distance (nm)')
ylabel('P(r)')
grid on,axis tight, box on
legend('Model','Fit')

subplot(2,2,4)
trueparam = [lam k r1 w1 r2 w2 amp1];
barh(100*(1 - param./trueparam))
tags = {'\lambda','k','<r_1>','\sigma_1','<r_2>','\sigma_2','A_1',};
set(gca,'yticklabel',tags)
ylabel('Model Parameters')
xlabel('relative fit error (%)')


