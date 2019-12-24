%% Computing the Pake pattern of a dipolar signal
% Author: Luis Fabregas
% Reviewed by: (Pending)
% A very basic example for displaying the frequency-domain Pake pattern of a 
% given dipolar signal.
%% Loading & Pre-processing
% Let's start by extracting the primary data from the spectrometer file using 
% |deerload. |

%Load experimental data
[traw,Vraw] = deerload('data/experiment_example.DTA');
traw = traw/1000; %ns -> us
%Pre-processing
V = correctphase(Vraw);
t = correctzerotime(V,traw);
V = correctscale(V,t);
%% Prepare the signal
% Since experimental dipolar signals contain the background, this must be fitted 
% removed prior to Fourier transform.
% 
% First we proceed to fit the background function using some time-domain parametric 
% model. In this example we will use a stretched exponential function (|td_strexp|). 
% Using the |fitbackground| function we obtain the fitted background as well as 
% the fitted modulation depth.

[B,lambda] = fitbackground(V,t,@td_strexp); 
%% 
% Now we can use these fitted variables to isolate the dipolar evolution function 
% from the primary data. Removal of the background via division leads to a noise 
% increase at later times and thus to an approximation |Dcorr| of the real dipolar 
% evolution function.

Dcorr = (V./B - (1 - lambda))/lambda;
%% Generating the Pake pattern
% Now that the signal has the appropiate structure for Fourier transform it, 
% we can call the |fftspec| function to obtained the Pake pattern.

[nu,pake] = fftspec(t,Dcorr);
%% 
% In order to avoid truncation ripples in the Fourier spectrum and at the same 
% time to compensate for the increase of noise, we recommend the use of apodization 
% using the appropiate option in |fftspec|.

[nu,pake] = fftspec(t,Dcorr,'Apodization',true);
%% 
% We do not need to worry about the zero-filling since |fftspec| takes care 
% of setting it to twice the amount of points in the signal, to preserve all information. 
% Adding more points will artificially increase the resolution of the Pake pattern. 
% The improvement will only be visual as no further information can be gained 
% from additional zero-filling.