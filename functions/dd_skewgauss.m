%
% DD_SKEWGAUSS Skew Gaussian distribution parametric model
%
%   info = DD_SKEWGAUSS
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_SKEWGAUSS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)   r0    3.5     1.0         20         location (nm)
% param(2)   w     0.5     0.2         5          FWHM (nm)
% param(3)  alpha  5.0     -15         15         skewness
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_skewgauss(r,param)

nParam = 3;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = 'Skew Gaussian distribution';
    info.nparam  = nParam;
    
    info.parameters(1).name = 'Location r0';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w';
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
 
    info.parameters(3).name = 'skewness alpha';
    info.parameters(3).range = [-25 25];
    info.parameters(3).default = 5;
    info.parameters(3).units = '';
    
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
r0 = param(1);
fwhm = param(2);
alpha = param(3);

sigma = fwhm/(2*sqrt(2*log(2)));
x = (r(:) - r0)/sigma/sqrt(2);

P = 1/sqrt(2*pi)*exp(-x.^2).*(1+erf(alpha*x));

if ~all(P==0)
    P = P/trapz(r,P);    
end
output = P;

return
