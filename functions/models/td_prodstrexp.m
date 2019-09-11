%
% TD_PRODSTREXP Product of two stretched exponentials background model
%
%   info = TD_SUMSTREXP
%   Returns an (info) structure containing the specifics of the model.
%
%   B = TD_SUMSTREXP(t,param)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% PARAM(1)  k1     3.5      0            200        1st strexp decay rate
% PARAM(2)  d1      3       0            6          1st strexp fractal dimension
% PARAM(3)  k2     3.5      0            200        2nd strexp decay rate
% PARAM(4)  d2      3       0            6          2nd strexp fractal dimension
% --------------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function output = td_prodstrexp(t,param)

nParam = 4;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Product of two stretched exponentials';
    info.Equation  = 'exp(-(k1*t)^(d1/3))*exp(-(k2*t)^(d2/3))';
    info.nparam  = nParam;
    info.parameters(1).name = 'Decay rate k1 of 1st stretched exponential';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'us^-1';
    
    info.parameters(2).name = 'Fractal dimension d1 of 1st stretched exponential';
    info.parameters(2).range = [0 6];
    info.parameters(2).default = 3;
    info.parameters(2).units = ' ';
    
    info.parameters(3).name = 'Decay rate k2 of 2nd stretched exponential';
    info.parameters(3).range = [0 200];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'us^-1';
    
    info.parameters(4).name = 'Fractal dimension d2 of 2nd stretched exponential';
    info.parameters(4).range = [0 6];
    info.parameters(4).default = 3;
    info.parameters(4).units = ' ';
    
    output = info;
    
end


%If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

%If necessary inputs given, compute the model distance distribution
t = abs(t);
StretchedExp1 = exp(-(param(1)*t).^(param(2)/3));
StretchedExp2 = exp(-(param(3)*t).^(param(4)/3));
Background = StretchedExp1.*StretchedExp2;
if ~iscolumn(Background)
    Background = Background';
end
output = Background;


return