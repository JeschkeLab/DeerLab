% %============================================================================
% DeerLab Example:
% Comparison of background treatment approaches
%============================================================================

% This script compares different approaches for treating the background
% (see the following publication for detailed information
%  https://doi.org/10.1039/C9CP06111H ).

clear, clc, clf

%===================
% Signal simulation
%===================

% Let's consider a dipolar signal arising from some Gaussian-shaped distance 
% distribution with a relatively fast decaying stretched exponential background. 
% We will set the modulation depth to 40% and add some noise to it.

% Model input parameters
decay = 0.4;     % decay constant for background (us^-1)
strfact = 1;     % stretching factor
lam = 0.4;       % modulation depth
r01 = 3.5;       % center of 1st Gaussian (nm)
w1 = 0.25;       % width of 1st Gaussian (nm)
A1 = 0.45;       % amplitude of 1st Gaussian
r02 = 4.0;       % center of 2nd Gaussian (nm)
w2 = 0.45;       % width of 2nd Gaussian (nm)

% Construct time and distance axes
t = linspace(0,5,251);
r = linspace(2,5,100);

% Generate distance distribution, background, and dipolar signal
P = dd_gauss2(r,[r01 w1 A1 r02 w2]);
B = bg_strexp(t,[decay strfact]);
V = dipolarsignal(t,r,P,lam,B,'noiselevel',0.01);

%======================
% Background treatment
%======================

% Now in order to treat the background we first need to know its parameters. 
% Usually we would have to fit the background either as a two-way approach (using 
% |fitbackground|) or by means of bilevel optimization. However, in this example 
% we will assume that we know the background parameters exactly to remove the 
% effects of background parameters misfits. 

%==========
% Division
%==========

% Knowing the background function we can first consider the approach of background 
% division. Once the signal has been processed into a dipolar-evolution function, 
% we us Tikhonov regularization to infer the distance distribution.

% "Correct" the signal by diving by background and rescaling
Vdiv = (V./B- (1-lam))/lam;

% Generate the dipolar kernel
K0 = dipolarkernel(t,r);

% Run Tikhonov regularization and get the time-domain fit
Pdiv = fitregmodel(Vdiv,K0,r,'tikh','aic');
Vfit = K0*Pdiv;

% Plotting
subplot(521)
plot(t,Vdiv,'k.',t,Vfit,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')

subplot(522)
plot(r,P,'k',r,Pdiv,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('V_i(t)')

%======================
% Subtraction
%======================
% Next we will try background subtraction. 

% "Correct" the signal for its background via subtraction
Vsub = (V - (1-lam)*B)/lam;

% Generate the dipolar kernel
K = dipolarkernel(t,r);

% Run Tikhonov regularization and get the time-domain fit
Psub = fitregmodel(Vsub,K,r,'tikh','aic');
Vfit = K*Psub;

% Plotting
subplot(523)
plot(t,Vsub,'k.',t,Vfit,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')

subplot(524)
plot(r,P,'k',r,Psub,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('V_i(t)')

%=======================
% Kernel with background
%=======================

% Next we includie the background information into the dipolar kernel to fit the
% primary data directly without the need for background correction.

% Generate the dipolar kernel with the background
K = dipolarkernel(t,r,lam,B);

% Run Tikhonov regularization and get the time-domain fit
Pker = fitregmodel(V,K,r,'tikh','aic');
Vfit = K*Pker;

%Plotting
subplot(525)
plot(t,V,'k.',t,Vfit,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')

subplot(526)
plot(r,P,'k',r,Pker,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('V_i(t)')

%============================================
% Kernel with square-root of background
%============================================

% As mentioned here, another approach involves the partial correction 
% of the dipolar signal by square-root of the backgrund and introduction of the 
% other squared-root of the background into the kernel.  

% Correct the signal, and generate the dipolar kernel with the background
Vsqrt = V./sqrt(B);
K = dipolarkernel(t,r,lam,sqrt(B));

% Run Tikhonov regularization and get the time-domain fit
Pker = fitregmodel(Vsqrt,K,r,'tikh','aic');
Vfit = K*Pker;

% Plotting
subplot(527)
plot(t,Vsqrt,'k.',t,Vfit,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')

subplot(528)
plot(r,P,'k',r,Pker,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('V_i(t)')

%============================================
% One-step fit
%============================================

% Finally, we analyze the time-domain signal in one step, without separating
% the background first. This is the preferred procedure.

[V1,P1,B1] = fitsignal(V,t,r,'P',@bg_strexp);

subplot(529)
plot(t,V,'k.',t,V1,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('time [\mus]')
ylabel('V_i(t)')

subplot(5,2,10)
plot(r,P,'k',r,P1,'b','LineWidth',1)
grid on,axis tight,box on
xlabel('distance [nm]')
ylabel('V_i(t)')
