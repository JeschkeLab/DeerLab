%% Fitting a custom time-domain model of a 4-pulse DEER signal
% Author: Luis Fabregas
% In this example we will construct a time-domain parametric model describing 
% our full 4-pulse DEER signal (that means distribution + background + other time-domain 
% parameters) using the built-in parametric models. We will then use this model 
% to fit our full parametric model to our signal in one step.
%% Generating a dipolar signal
% Let's start by simulating a dipolar signal, whose parameters we know. For 
% this example we will use a dipolar signal with a modulation depth of 27% originating 
% from a Gaussian distance distribution (mean distance of 3.5 $nm$ and width of 
% 0.3 $nm$) with a stretched exponential background model with a decay rate of 
% 0.15$\mu s^{-1}$ and a dimensionality of 2.8.

clc,clf,clear

%Define the parameters
rmean = 3.5; %nm
w = 0.3; %nm
lam = 0.27;
k = 0.15; %us-1
d = 2.8;

par = [lam rmean w k d];
%% 
% We will simulate the dipolar signal on a time-axis defined between -0.2$\mu 
% s$ and 4.0$\mu s$

t = linspace(-0.2,4,300);
r = time2dist(t);
%% 
% Now we can generate the signal with some added noise

%Generate the distance distribution
P = rd_onegaussian(r,[rmean w]);
%Generate the background
B = td_strexp(t,[k d]);
%Generate the signal
V = dipolarsignal(t,r,P,'moddepth',lam,'background',B,'noiselevel',0.01);

plot(t,V)
xlabel('time [\mus]'),ylabel('V(t)')
%% Defining the time-domain model
% This is probably the most crucial step: defining the model. How well a signal 
% can be fitted depends mainly on how well the model can describe our signal. 
% In order to fit our time-domain dipolar signal in one step, we require a model 
% which contains the information on the distance distribution, the background 
% and modlation depth. Even tough such a mode is not implemented as a function 
% in the program, we can easily define our own model using the built-in parametric 
% models.
% 
% In this case we know the groun truth, so we will construct our time-domain 
% model according to 
% 
% $$\mathit{\mathbf{V}}=\left(\left(1-\lambda \;\right)+\;\lambda \textrm{KP}\right)\odot 
% \mathit{\mathbf{B}}\;$$
% 
% where the background $\mathbf{B}$ will be modeled by a stretched exponential 
% time-domain model (|td_strexp|) and the distribution $\mathbf{P}$ will be modeled 
% by a single Gaussian distribution (|rd_onegaussian|)
% 
% We can build the model by defining an anonymous function handle |model|, which 
% must first accept the time axis |t| and the parameters |par|. The parameters 
% of our model must be defined on a vector |par|
%% 
% * |par(1) -> lam|
% * |par(2) -> rmean|
% * |par(3) -> w|
% * |par(4) -> k|
% * |par(5) -> d|

%Define our time-domain model
model = @(t,par) (1 - par(1) + par(1)*dipolarkernel(t,r)*rd_onegaussian(r,par(2:3))).*td_strexp(t,par(4:5));
%% Fitting the model
% At this point our model is ready to be fitted. The only remaining step remaining 
% is to define the initial values of the model parameters. Since it is a custom 
% model, all parameter values will be unbounded so it is recommended to set some 
% upper/lower boundaries for them. All three values can be defined as vectors 
% with the same order as in the |par| vector. All these values can be well-estimated 
% by visual inspection of the data.

% Parameter lam rmean  w   k   d
% Index      1    2    3   4   5
par0 =     [0.35 4.0 0.4  0.1 3.0];
lower =    [0.10 2.0 0.1  0.0 2.0];
upper =    [0.50 7.0 0.5  0.3 4.0];
%% 
% Finally we can launch the fitting function to obtain our fitted model parameters.

%Launch the fit
fitpar = fitparamodel(V,model,t,par0,'lower',lower,'upper',upper);
%% 
% We can now obtain the fitted models, and compare the fit parameters to our 
% input parameters:

%Get the fits
Vfit = model(t,fitpar);
Pfit = rd_onegaussian(r,fitpar(2:3));

%Plot results
subplot(121)
plot(t,V,'k.',t,Vfit)
subplot(122)
plot(r,P,'k',r,Pfit)

round(par,2)
round(fitpar,2)