%
% RD_FOURGAUSSIAN Sum of four Gaussian distributions parametric model
%
%   info = RD_FOURGAUSSIAN
%   Returns an (info) structure containing the specifics of the model.
%
%   P = RD_FOURGAUSSIAN(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <r1>   2.5     1.5         20         1st mean distance
% param(2)  w1     0.5     0.2         5          FWHM of 1st distance
% param(3)  <r2>   3.5     1.5         20         2nd mean distance
% param(4)  w2     0.5     0.2         5          FWHM of 2nd distance
% param(5)  <r3>   5.0     1.5         20         3rd mean distance
% param(6)  w3     0.5     0.2         5          FWHM of 3rd distance
% param(7)  <r4>   5.0     1.5         20         4th mean distance
% param(8)  w4     0.5     0.2         5          FWHM of 4th distance
% param(9)  p1     0.25    0           1          Fraction of pairs at 1st distance
% param(10) p2     0.25    0           1          Fraction of pairs at 2nd distance
% param(11) p3     0.25    0           1          Fraction of pairs at 3rd distance
%--------------------------------------------------------------------------
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function output = rd_fourgaussian(r,param)

nParam = 8;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Three-Gaussian distribution';
    info.Equation  = ['A1*exp(-(r-<r1>)²/(',char(963),'1*sqrt(2))²) + A2*exp(-(r-<r2>)²/(',char(963),'2*sqrt(2))²) + (1-A1-A2)*exp(-(r-<r3>)²/(',char(963),'3*sqrt(2))²)'];
    info.nParam  = nParam;
    info.parameters(1).name = 'Mean distance <r1> of 1st Gaussian';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w1 of 1st Gaussian';
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'Mean distance <r2>  of 2nd Gaussian';
    info.parameters(3).range = [1 20];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = 'FWHM w2 of 2nd Gaussian';
    info.parameters(4).range = [0.2 5];
    info.parameters(4).default = 0.5;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Mean distance <r3> of 3rd Gaussian';
    info.parameters(5).range = [1 20];
    info.parameters(5).default = 3.5;
    info.parameters(5).units = 'nm';
    
    info.parameters(6).name = 'FWHM w3 of 3rd Gaussian';
    info.parameters(6).range = [0.2 5];
    info.parameters(6).default = 0.5;
    info.parameters(6).units = 'nm';
    
    info.parameters(7).name = 'Mean distance <r4> of 4th Gaussian';
    info.parameters(7).range = [1 20];
    info.parameters(7).default = 3.5;
    info.parameters(7).units = 'nm';
    
    info.parameters(8).name = 'FWHM w4 of 4th Gaussian';
    info.parameters(8).range = [0.2 5];
    info.parameters(8).default = 0.5;
    info.parameters(8).units = 'nm';
    
    info.parameters(9).name = 'Relative amplitude A1 of 1st Gaussian';
    info.parameters(9).range = [0 1];
    info.parameters(9).default = 0.25;
    
    info.parameters(10).name = 'Relative amplitude A2 of 2nd Gaussian';
    info.parameters(10).range = [0 1];
    info.parameters(10).default = 0.25;
    
    info.parameters(11).name = 'Relative amplitude A3 of 3rd Gaussian';
    info.parameters(11).range = [0 1];
    info.parameters(11).default = 0.25;
       
    output = info;
    return
end

if nargin~=2
    error('Model requires two input arguments.')
end
    
% Check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Compute the model distance distribution
Gaussian1 = exp(-((r(:)-param(1))/(param(2))).^2);
Gaussian2 = exp(-((r(:)-param(3))/(param(4))).^2);
Gaussian3 = exp(-((r(:)-param(5))/(param(6))).^2);
Gaussian4 = exp(-((r(:)-param(7))/(param(8))).^2);
P = param(9)*Gaussian1 + param(10)*Gaussian2 + param(11)*Gaussian3 + max(1 - param(9) - param(10) - param(11),0)*Gaussian4;

% Normalize
dr = r(2)-r(1);
P = P/sum(P)/dr;

output = P;

return