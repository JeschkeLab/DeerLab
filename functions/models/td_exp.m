%
% TD_EXP Exponential background model
%
%   info = TD_EXP
%   Returns an (info) structure containing the specifics of the model.
%
%   B = TD_EXP(t,param)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% ------------------------------------------------------------------
% PARAM(1)  k    0.35      0            200        decay rate
% ------------------------------------------------------------------
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.



function output = td_exp(t,param)

nParam = 1;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Exponential';
    info.nparam  = nParam;
    info.parameters(1).name = 'Decay rate k';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 0.35;
    info.parameters(1).units = 'us^-1';
    
    output = info;
    return
end

%If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

%If necessary inputs given, compute the model distance distribution
t = abs(t);
Background = exp(-param(1)*t);
if ~iscolumn(Background)
    Background = Background';
end
output = Background;


return