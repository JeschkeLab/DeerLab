%
% DD_UNIFORM Uniform distribution parametric model
%
%   info = DD_UNIFORM
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_UNIFORM(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  rL     2.5     0.1         6           left edge
% param(2)  rR     3.0     0.2         20          right edge
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = dd_uniform(r,param)

nParam = 2;

if nargin~=0 && nargin~=2 
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Uniform distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'left edge rL';
    info.parameters(1).range = [0.1 6];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'right edge rR';
    info.parameters(2).range = [0.2 20];
    info.parameters(2).default = 3.0;
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
rL = min(abs(param));
rR = max(abs(param));
P = zeros(numel(r),1);
P(r>=rL & r<=rR) = 1;
P = P/sum(P)/mean(diff(r));

output = P;

return