function output = rd_twogaussian(r,param)
%
% TWOGAUSSIAN Sum of two Gaussian distributions parametric model
%
%   info = TWOGAUSSIAN
%   Returns an (info) structure containing the specifics of the model.
%
%   P = TWOGAUSSIAN(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <r1>   2.5     1.5         20         1st mean distance
% param(2)  s(r1)  0.5     0.2         5          std. dev. of 1st distance
% param(3)  <r2>   3.5     1.5         20         2nd mean distance
% param(4)  s(r2)  0.5     0.2         5          std. dev. of 2nd distance
% param(5)  A      0.5     0           1          fraction of pairs at 1st distance
% --------------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

nParam = 5;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Two-Gaussian distribution';
    info.Equation  = ['A*exp(-(r-<r1>)²/(',char(963),'1*sqrt(2))²) + (1-A)*exp(-(r-<r2>)²/(',char(963),'2*sqrt(2))²)'];
    info.nParam  = nParam;
    info.parameters(1).name = 'Mean distance <r1> 1st Gaussian';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963),'1 1st Gaussian'];
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'Mean distance <r2> 2nd Gaussian';
    info.parameters(3).range = [1 20];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = ['Standard deviation ',char(963),'2 2nd Gaussian'];
    info.parameters(4).range = [0.2 5];
    info.parameters(4).default = 0.5;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Relative amplitude A 1st Gaussian';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end    
    
    %If necessary inputs given, compute the model distance distribution
    Gaussian1 = exp(-((r-param(1))/(param(2))).^2);
    Gaussian2 = exp(-((r-param(3))/(param(4))).^2);
    Distribution = param(5)*Gaussian1 + max(1 - param(5),0)*Gaussian2;
    if ~iscolumn(Distribution)
        Distribution = Distribution';
    end
    %Normalize
    Distribution = Distribution/sum(Distribution)/mean(diff(r));
    output = Distribution;
else
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end

return