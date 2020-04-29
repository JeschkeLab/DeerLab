%
% BG_SUMSTREXP Sum of two stretched exponentials background model
%
%   info = BG_SUMSTREXP
%   Returns an (info) structure containing the specifics of the model.
%
%   B = BG_SUMSTREXP(t,param)
%   B = BG_SUMSTREXP(t,param,lambda)
%   Computes the N-point model (B) from the N-point time axis (t) according to
%   the paramteres array (param). The required parameters can also be found
%   in the (info) structure. The pathway amplitude (lambda) can be
%   included, if not given the default lambda=1 will be used.
%
% PARAMETERS
% name    symbol default lower bound upper bound
% --------------------------------------------------------------------------
% PARAM(1) kappa1   3.5      0            200        1st strexp decay rate
% PARAM(2)  d1      3        0            6          1st strexp fractal dimension
% PARAM(3) kappa2   3.5      0            200        2nd strexp decay rate
% PARAM(4)  d2      3        0            6          2nd strexp fractal dimension
% PARAM(5)  A1      0.5      0            1          Relative amplitude
% --------------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function output = bg_sumstrexp(t,param,lambda)

nParam = 5;

if all(nargin~=[0 2 3])
    error('Model requires at least two input arguments.')
end

if nargin==0
    %If no inputs given, return info about the parametric model
    info.model  = 'sum of two stretched exponentials';
    info.nparam  = nParam;
    info.parameters(1).name = 'decay rate kappa1 of 1st stretched exponential';
    info.parameters(1).range = [0 200];
    info.parameters(1).default = 3.5;
    info.parameters(1).units = 'us^-1';
    
    info.parameters(2).name = 'fractal dimension d1 of 1st stretched exponential';
    info.parameters(2).range = [0 6];
    info.parameters(2).default = 3;
    info.parameters(2).units = ' ';
    
    info.parameters(3).name = 'decay rate kappa2 of 2nd stretched exponential';
    info.parameters(3).range = [0 200];
    info.parameters(3).default = 3.5;
    info.parameters(3).units = 'us^-1';
    
    info.parameters(4).name = 'fractal dimension d2 of 2nd stretched exponential';
    info.parameters(4).range = [0 6];
    info.parameters(4).default = 3;
    info.parameters(4).units = ' ';
    
    info.parameters(5).name = 'relative amplitude of 1st stretched exponential';
    info.parameters(5).range = [0 1];
    info.parameters(5).default = 0.5;
    info.parameters(5).units = ' ';
    
    output = info;
    return

end

if nargin<3
    lambda = 1;
end

% If user passes them, check that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% If necessary inputs given, compute the model background
kappa1 = param(1);
d1 = param(2);
kappa2 = param(3);
d2 = param(4);
w1 = param(5);
strexp1 = exp(-lambda*kappa1*abs(t).^(d1/3));
strexp2 = exp(-lambda*kappa2*abs(t).^(d2/3));
B = w1*strexp1 + (1-w1)*strexp2;
B = B(:);
output = B;


return