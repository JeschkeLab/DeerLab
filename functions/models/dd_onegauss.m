
%
% DD_ONEGAUSS Gaussian distribution parametric model
%
%   info = DD_ONEGAUSS
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_ONEGAUSS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <r>    3.5     1.0         20         mean distance
% param(2)   w     0.5     0.2         5          FWHM
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function output = dd_onegauss(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Single Gaussian distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'Mean distance <r>';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w';
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
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
P = sqrt(2/pi)*1/sigma*exp(-((r(:) - param(1))/(sqrt(2)*sigma)).^2);
dr = r(2)-r(1);
if ~all(P==0)
P = P/sum(P)/dr;    
end
output = P;

return