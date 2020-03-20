%
% DD_THREEGAUSS Sum of three Gaussian distributions parametric model
%
%   info = DD_THREEGAUSS
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_THREEGAUSS(r,param)
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
% param(7)  p1     0.3     0           1          fraction of pairs at 1st distance
% param(8)  p2     0.3     0           1          fraction of pairs at 2nd distance
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function output = dd_threegauss(r,param)

nParam = 8;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Three-Gaussian distribution';
    info.nparam  = nParam;
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
    
    info.parameters(4).name = ['FWHM w2 of 2nd Gaussian'];
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
    
    info.parameters(7).name = 'Relative amplitude A1 of 1st Gaussian';
    info.parameters(7).range = [0 1];
    info.parameters(7).default = 0.3;
    
    info.parameters(8).name = 'Relative amplitude A2 of 2nd Gaussian';
    info.parameters(8).range = [0 1];
    info.parameters(8).default = 0.3;
       
    output = info;
    return
end
    
% Check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

%Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
Lam1 = (param(2)/sqrt(2*log(2)));
Lam2 = (param(4)/sqrt(2*log(2)));
Lam3 = (param(6)/sqrt(2*log(2)));
Gaussian1 = sqrt(2/pi)*1/Lam1*exp(-2*((r(:) - param(1))/Lam1).^2);
Gaussian2 = sqrt(2/pi)*1/Lam2*exp(-2*((r(:) - param(3))/Lam2).^2);
Gaussian3 = sqrt(2/pi)*1/Lam3*exp(-2*((r(:) - param(5))/Lam3).^2);
P = param(7)*Gaussian1 + param(8)*Gaussian2 + max(1 - param(7) - param(8),0)*Gaussian3;

% Normalize
dr = r(2)-r(1);
if ~all(P==0)
P = P/sum(P)/dr;    
end
output = P;

return