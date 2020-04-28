%
% BG_EXP Exponential background model
%
%   info = BG_EXP
%   Returns an (info) structure containing the specifics of the model.
%
%   B = BG_EXP(t,param)
%   B = BG_EXP(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% ------------------------------------------------------------------
% PARAM(1) kappa    0.35      0            200        decay rate
% ------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.



function output = bg_exp(t,param,lambda)

nParam = 1;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = 'Exponential';
    info.nparam  = nParam;
    info.parameters(1).name = 'Decay rate kappa';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 0.35;
    info.parameters(1).units = 'us^-1';
    
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

% If necessary inputs given, compute the model distance distribution
kappa = param(1);
B = exp(-lambda*kappa*abs(t));
B = B(:);
output = B;


return