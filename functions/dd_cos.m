%
% DD_COS Raised-cosine parametric model
%
%   info = DD_COS
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_COS(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  r0     3.0     0.1        20         center
% param(2)  fwhm   0.5     0.1         5         FWHM
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_cos(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Raised-cosine distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'center r0';
    info.parameters(1).range = [0.1 20];
    info.parameters(1).default = 3.0;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w';
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
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

phi = (r-r0)/fwhm*pi;
P = (1+cos(phi))/2/fwhm;
P(r<r0-fwhm | r>r0+fwhm) = 0;

if any(P~=0)
    P = P/sum(P)/mean(diff(r));
end
P = P(:);
output = P;

return
