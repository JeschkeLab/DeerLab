%
% DD_GAUSS2 Sum of two Gaussian distributions parametric model
%
%   info = DD_GAUSS2
%   Returns an (info) structure containing the specifics of the model.
%
%   P = DD_GAUSS2(r,param)
%   Computes the N-point model (P) from the N-point distance axis (r) according to 
%   the paramteres array (param). The required parameters can also be found 
%   in the (info) structure.
%
% PARAMETERS
% name      symbol default lower bound upper bound
% --------------------------------------------------------------------------
% param(1)  <r1>   2.5     1.0         20         center of 1st Gaussian
% param(2)  w1     0.5     0.2         5          FWHM of 1st Gaussian
% param(3)  <r2>   3.5     1.0         20         center of 2nd Gaussian
% param(4)  w2     0.5     0.2         5          FWHM of 2nd Gaussian
% param(5)  A      0.5     0           1          amplitude of 1st Gaussian
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function output = dd_gauss2(r,param)

nParam = 5;

if nargin~=0 && nargin~=2
    error('Model requires two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'Two-Gaussian distribution';
    info.nparam  = nParam;
    info.parameters(1).name = 'Mean distance <r1> of 1st Gaussian';
    info.parameters(1).range = [1 20];
    info.parameters(1).default = 2.5;
    info.parameters(1).units = 'nm';
    
    info.parameters(2).name = 'FWHM w1 of 1st Gaussian';
    info.parameters(2).range = [0.2 5];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = 'nm';
    
    info.parameters(3).name = 'Mean distance <r2> 2nd Gaussian';
    info.parameters(3).range = [1 20];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'nm';
    
    info.parameters(4).name = 'FWHM w2 of 2nd Gaussian';
    info.parameters(4).range = [0.2 5];
    info.parameters(4).default = 0.5;
    info.parameters(4).units = 'nm';
    
    info.parameters(5).name = 'Relative amplitude A 1st Gaussian';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    output = info;
    return
end

% Check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Compute the model distance distribution
fwhm = param([2 4]);
r0 = param([1 3]);
a = param(5);
a(2) = max(1-a,0);
P = multigaussfun(r,r0,fwhm,a);

output = P;

return
