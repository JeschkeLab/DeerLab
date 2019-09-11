%
% TD_STREXP Stretched exponential background model
%
%   info = TD_STREXP
%   Returns an (info) structure containing the specifics of the model.
%
%   B = TD_STREXP(t,param)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% ------------------------------------------------------------------
% PARAM(1)  k     0.25      0            200        decay rate
% PARAM(2)  d      3       0            6          fractal dimension
% ------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function output = td_strexp(t,param)

nParam = 2;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Stretched exponential';
    info.nparam  = nParam;
    info.parameters(1).name = 'Decay rate k';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 0.25;
    info.parameters(1).units = 'us^-1';
    
    info.parameters(2).name = 'Fractal dimension d';
    info.parameters(2).range = [0 6];
    info.parameters(2).default = 3;
    info.parameters(2).units = ' ';
    
    output = info;
    return
end

%If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

%If necessary inputs given, compute the model distance distribution
t = abs(t);
Background = exp(-(param(1)*t).^(param(2)/3));
if ~iscolumn(Background)
    Background = Background';
end
output = Background;

return