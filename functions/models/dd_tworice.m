function output = dd_tworice(r,param)
%
% TWORICE Sum of two rician distributions parametric model
%
%   info = TWORICE
%   Returns an (info) structure containing the specifics of the model.
%
%   P = TWORICE(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  nu1      2.5     1.0        10         mean distance, 1st component
% param(2)  sigma1   0.4     0.1        5          standard deviation, 1st component
% param(3)  nu2      4.0     1.0        10         mean distance, 2nd component
% param(4)  sigma2   0.4     0.1        5          standard deviation, 2nd component
% param(5)  p1       0.5     0          1          fraction of 1st component
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


nParam = 5;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = 'Two Rice/Rician distributions';
    info.nparam  = nParam;
    info.parameters(1).name = ['Mean distance ',char(957),'1 1st Rician'];
    info.parameters(1).range = [1 10];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = ['Standard deviation ',char(963),'1 1st Rician'];
    info.parameters(2).range = [0.1 5];
    info.parameters(2).default = 0.7;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = ['Mean distance ',char(957),'2 2nd Rician'];
    info.parameters(3).range = [1 10];
    info.parameters(3).default = 4.0;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = ['Standard deviation ',char(963),'2 2nd Rician'];
    info.parameters(4).range = [0.1 5];
    info.parameters(4).default = 0.7;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Relative amplitude A 1st Rician';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    
    output = info;
    return;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute non-central chi distribution with 3 degrees of freedom (a 3D Rician)
nu1 = param(1);
sig1 = param(2);
nu2 = param(3);
sig2 = param(4);
p1 = param(5);
p2 = 1-p1;
P = p1*rice3d(r,nu1,sig1) + p2*rice3d(r,nu2,sig2);
P = P(:);

if ~all(P==0)
    P = P/sum(P)/mean(diff(r));
end

output = P;

return