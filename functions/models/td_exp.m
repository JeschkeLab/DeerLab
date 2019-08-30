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
% PARAM(1)  k     3.5      0            200        decay rate
% ------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function output = td_exp(t,param)

nParam = 1;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Exponential';
    info.Equation  = ['exp(-kt)'];
    info.nParam  = nParam;
    info.parameters(1).name = 'Decay rate k';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'us^-1';
    
    output = info;
    
elseif nargin == 2
    
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
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end

return