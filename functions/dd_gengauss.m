%
% DD_GENGAUSS Generalized Gaussian distribution parametric model
%
%   info = DD_GENGAUSS
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_GENGAUSS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)   r0    3.5     1.0         20         location (nm)
% param(2)   w     0.5     0.2         5          FWHM (nm)
% param(3)  beta   5.0     0.25        15         kurtosis
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_gengauss(r,param)

nParam = 3;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Skew Gaussian distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'Mean distance <r>';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w';
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'kurtosis beta';
    info.parameters(3).range = [0.25 15];
    info.parameters(3).default = 5;
    info.parameters(3).units = '';
    
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
  error('The number of input parameters does not match the number of model parameters.')
end

%Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
sigma = (param(2)/(2*sqrt(2*log(2))));
beta = param(3);
x = abs(r(:) - param(1))/sigma;
P = beta/(2*sigma*gamma(1/beta))*exp(-x.^beta);
dr = r(2)-r(1);
if ~all(P==0)
P = P/sum(P)/dr;    
end
output = P;

return