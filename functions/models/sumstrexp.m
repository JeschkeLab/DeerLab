%
% SUMSTREXP Sum of two stretched exponentials background model
%
%   INFO = SUMSTREXP
%   Returns an INFO structure containing the specifics of the model.
%
%   X = SUMSTREXP(T,PARAM)
%   Comptues the N-point model from the N-point time axis T according to
%   the paramteres array PARAM. The required parameters can also be found
%   in the INFO structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% PARAM(1)  k1     3.5      0            200        1st strexp decay rate
% PARAM(2)  d1      3       0            6          1st strexp fractal dimension
% PARAM(3)  k2     3.5      0            200        2nd strexp decay rate
% PARAM(4)  d3      3       0            6          2nd strexp fractal dimension
% PARAM(5)  A1      0.5     0            1          Relative amplitude
% --------------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function output = sumstrexp(t,param)

nParam = 5;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Sum of two stretched exponentials';
    info.Equation  = 'A1*exp(-(k1*t)^(d1/3)) + (1-A1)*exp(-(k2*t)^(d2/3))';
    info.nParam  = nParam;
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
    
    info.parameters(5).name = 'Relative amplitude of 1st stretched exponential';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    info.parameters(5).units = ' ';
    
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end
    
    %If necessary inputs given, compute the model distance distribution
    StretchedExp1 = exp(-(param(1)*t).^(param(2)/3));
    StretchedExp2 = exp(-(param(3)*t).^(param(4)/3));
    Background = param(5)*StretchedExp1 + (1-param(5))*StretchedExp2;
    if ~iscolumn(Background)
        Background = Background';
    end
    output = Background;
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end

return