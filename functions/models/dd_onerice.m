function output = dd_onerice(r,param)
%
% DD_ONERICE 3D-Rice distribution parametric model
%
%   info = DD_ONERICE
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_ONERICE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name     symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  nu      3.5     1.0         10          center
% param(2)  sigma   0.7     0.1          5          width
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


nParam = 2;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = 'Single 3D-Rice distribution';
    info.nparam  = nParam;
    info.parameters(1).name = ['Center ',char(957)];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Width ',char(963)];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    output = info;
    return
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute non-central chi distribution with 3 degrees of freedom (a 3D Rician)
nu = param(1);
sig = param(2);
a = 1;
P = multirice3d(r,nu,sig,a);

if ~all(P==0)
    P = P/sum(P)/mean(diff(r));
end

output = P;

return
