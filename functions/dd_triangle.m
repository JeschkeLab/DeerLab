%
% DD_TRIANGLE Triangle distribution parametric model
%
%   info = DD_TRIANGLE
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_TRIANGLE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  r0     3.5     1.0         20         mode
% param(2)  wL     0.3     0.1         5          left width
% param(3)  wR     0.3     0.1         5          right width
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_triangle(r,param)

nParam = 3;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Triangle distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'Center distance r0';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'width left wL';
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.3;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'width right wR';
    info.parameters(3).range = [0.1 5];
    info.parameters(3).default = 0.3;
    info.parameters(3).units = 'nm';
    
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
wL = abs(param(2));
wR = abs(param(3));
rL = r0 - wL;
rR = r0 + wR;
idxL = r>=r0-wL & r<=r0;
idxR = r<=r0+wR & r>=r0;
P = zeros(numel(r),1);
if wL>0
    P(idxL) = (r(idxL)-rL)/wL/(wL+wR);
end
if wR>0
    P(idxR) = -(r(idxR)-rR)/wR/(wL+wR);
end

if any(P~=0)
    P = P/trapz(r,P);
end

output = P;

return