%================================================================
% DeerAnalyis2
% Example: 5-pulse DEER fitting
% Fit data obtained from 5-pulse DEER experiments without any correction
%================================================================

clear,clc,clf

%Preparation
%----------------------------------------------
N = 100;
t = linspace(-0.3,2,N);
r = time2dist(t);
lambda = 1;
%Distribution with two Gaussians
P = rd_twogaussian(r,[2.5 0.3 4 0.2 0.4]);

%Generate dipolar signal without artefact
S = dipolarsignal(t,r,P);
V0 = (1-lambda) + lambda*S;

%Generate dipolar signal with artefact
tshift = max(t/2); %Time shift of the artefact
Amp = 0.4; %Relative amplitude of the artefact
[~,pos] = min(abs(t - tshift));
t5pulse = t - tshift;
S5pulse = dipolarsignal(t5pulse,r,P);
V = S + Amp*S5pulse;
V = V/max(V);
V = (1-lambda) + lambda*V;


%Add noise to both signals
V = V + whitegaussnoise(N,0.02);
V0 = V0 + whitegaussnoise(N,0.02);

fitartefact = false;
[Vfit0,Pfit0] = my5pDEER(t,[],r,V0,fitartefact);

%Fitting (with 5-pulse DEER artefact)
%----------------------------------------------
fitartefact = true;
%Create function handle depending on r and param from the custom model
fcnhandle = @(t,param)my5pDEER(t,param,r,V,fitartefact);
%Launch the fitting of the B-parametric model + Tikhonov regularization
[~,parafit] = fitparamodel(V,fcnhandle,t,[0.3],'Lower',[0],'Upper',[1]);
%Obtain the fitted signal and distance distribution
[Vfit,Pfit] = my5pDEER(t,parafit,r,V,fitartefact);


%Plot results
%----------------------------------------------
subplot(221)
plot(r,P,'k',r,Pfit0,'b','LineWidth',1.5)
box on, grid on, axis tight
xlabel('Distance [nm]')
ylabel('P(r)')
title('4-pulse DEER')
subplot(222)
plot(t,V0,'k',t,Vfit0,'b','LineWidth',1.5)
xlabel('Time [\mus]')
ylabel('V(t)')
box on, grid on, axis tight
title('4-pulse DEER')

subplot(223)
plot(r,P,'k',r,Pfit,'b','LineWidth',1.5)
box on, grid on, axis tight
xlabel('Distance [nm]')
ylabel('P(r)')
title('5-pulse DEER')

subplot(224)
plot(t,V,'k',t,Vfit,'b','LineWidth',1.5)
xlabel('Time [\mus]')
ylabel('V(t)')
box on, grid on, axis tight
title('5-pulse DEER')

%Definition of the custom model
%----------------------------------------------

function [Vfit,Pfit] = my5pDEER(t,param,r,V,fitartefact)

%Fit the modulation depth as first parameter...

%Construct a kernel with the fitted background
if fitartefact
    K = dipolarkernel(t,r,'FivePulseCoeff',param(1));
else
    K = dipolarkernel(t,r);
end
%Prepare regularization
L = regoperator(length(V),2);
alphas = regparamrange(K,L);
alpha = selregparam(alphas,V,K,L,'tikh','aic');
%Regularize the data using the fitted backgorund
Pfit = fitregmodel(V,K,r,L,'tikhonov',alpha);
%Get the signal for comparison in time-domain
Vfit = K*Pfit;

end


