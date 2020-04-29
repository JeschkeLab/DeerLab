%
% BG_STREXP Stretched exponential background model
%
%   info = BG_STREXP
%   Returns an (info) structure containing the specifics of the model.
%
%   B = BG_STREXP(t,param)
%   B = BG_STREXP(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% ------------------------------------------------------------------
% PARAM(1) kappa   0.25        0           200         decay rate
% PARAM(2)   d       1         0            6          stretch factor
% ------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.



function output = bg_strexp(t,param,lambda)

nParam = 2;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'stretched exponential';
    info.nparam  = nParam;
    info.parameters(1).name = 'decay rate kappa';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 0.25;
    info.parameters(1).units = 'us^-1';
    
    info.parameters(2).name = 'fractal dimension d';
    info.parameters(2).range = [0 6];
    info.parameters(2).default = 1;
    info.parameters(2).units = ' ';
    
    output = info;
    return
end

if nargin<3
    lambda = 1;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% If necessary inputs given, compute the model background
kappa = param(1);
d = param(2);
Background = exp(-lambda*kappa*abs(t).^d);
if ~iscolumn(Background)
    Background = Background';
end
output = Background;

return