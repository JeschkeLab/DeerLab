function output = rd_onegaussian(r,param)
%
% ONEGAUSSIAN Gaussian distribution parametric model
%
%   info = ONEGAUSSIAN
%   Returns an (info) structure containing the specifics of the model.
%
%   P = ONEGAUSSIAN(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <r>    3.5     1.0         20         mean distance
% param(2)   s     0.5     0.02        5          standard deviation
% --------------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

nParam = 2;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Single Gaussian distribution';
    info.Equation  = ['exp(-(r-<r>)²/(',char(963),'*sqrt(2))²)'];
    info.nParam  = nParam;
    info.parameters(1).name = 'Mean distance <r>';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963)];
    info.parameters(2).range = [0.02 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end    
    
    %If necessary inputs given, compute the model distance distribution
    Distribution = exp(-((r-param(1))/(param(2))).^2);
    if ~iscolumn(Distribution)
        Distribution = Distribution';
    end
    Distribution = Distribution/sum(Distribution)/mean(diff(r));
    output = Distribution;
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end

return