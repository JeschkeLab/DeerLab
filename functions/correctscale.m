% 
% CORRECTSCALE Amplitude scale correction
%
%       [Vc,V0] = CORRECTSCALE(V,t) 
%       [Vc,V0] = CORRECTSCALE(V,t,tmax) 
%       [Vc,V0] = CORRECTSCALE(V,t,tmax,model) 
%
%  Takes the experimental dipolar signal (V) on a given time axis (t) and
%  rescales it so that it is 1 at time zero, taking into account the noise.
%  It fits a model specified in (model) over the interval |t|=<tmax around time zero.
%  It returns the rescaled signal (Vc) and the scaling factor (V0).
%
%  Inputs:
%    V     experimenal signal (vector)
%    t     time axis, in microseconds (vector)
%    tmax  time cutoff for fit range, in microseconds (scalar), default 0.2
%    model model to fit over |t|<tmax
%            'deer' DEER model with Gaussian distribution (default)
%            'gauss'    raised Gaussian
%
%  Outputs:
%    Vc    rescaled signal
%    V0    rescaling factor, such that Vc = V/V0
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [Vc,V0] = correctscale(V,t,tmax,model)

if ~isreal(V)
   error('Input signal cannot be complex.') 
end

if nargin<3 || isempty(tmax)
    tmax = 0.2;
end
if nargin<4 || isempty(model)
    model = 'deer';
end

% Validate input
validateattributes(V,{'numeric'},{'nonempty','vector'},mfilename,'V')
validateattributes(t,{'numeric'},{'nonempty','increasing','vector'},mfilename,'t')
validateattributes(tmax,{'numeric'},{'scalar','positive'},mfilename,'tmax')
validateattributes(model,{'char'},{'nonempty'},mfilename,'model')

% Use column vectors
V = V(:);
t = t(:);

% Generate distance axis
r = time2dist(t);

% Approximate scaling
Amp0 = max(V);
V_ = V/Amp0;

% Limit fit to region around time zero
idx = abs(t)<=tmax;
V_ = V_(idx);
t_ = t(idx);

switch lower(model)
    case 'deer'
        % DEER model with Gaussian distibution
        fitmodel = @(tt,p) p(1)*bg_exp(tt,p(5)).*((1-p(2)) + p(2)*dipolarkernel(tt,r)*dd_gauss(r,p(3:4)));
        par0 = [1 0.5 3 0.3 0.2];
        lb = [1e-3 1e-5 0 0 0];
        ub = [10 1 20 5 100];
    case 'gauss'
        % Gaussian function
        fitmodel = @(t,p) (p(1)-p(2))+p(2)*exp(-t.^2/p(3)^2);
        par0 = [1 1 1];
        lb = [1e-3 1e-3 1e-3];
        ub = [10 1000 10000];
    otherwise
        error('Unknown model ''%s''.',modeltype);
end

% Run the parametric model fitting
parfit = fitparamodel(V_,fitmodel,t_,par0,'Upper',ub,'Lower',lb);

% Get the fitted signal amplitude and scale the signal
V0 = Amp0*parfit(1);
Vc = V/V0;

% Plotting
doPlotting = nargout==0;
if doPlotting
    t_ = linspace(min(t_),max(t_),numel(t_)*4-3);
    Vfit_ = fitmodel(t_,parfit);
    plot(t,V,'.',t_,Amp0*Vfit_);
    yline(V0,'r');
    d = max(V)-min(V);
    axis tight
    ylim([min(V) max(V)]+[-1 1]*d/20);
    xlabel('time (us)');
    ylabel('V');
    grid on
end

end
