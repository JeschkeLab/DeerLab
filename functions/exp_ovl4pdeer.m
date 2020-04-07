%
% EXP_OVL4PDEER 4-pulse DEER with band overlap experiment model 
%
%   info = EXP_OVL4PDEER(t)
%   Returns an (info) structure containing the specifics of the model.
%
%   [K,B] = EXP_OVL4PDEER(t,r,param,Bmodel)
%   Computes the NxM dipolar kernel (K) and N-point multi-pathway background (B) 
%   from the N-point time axis (t) and M-point distance axis (r) according to
%   the paramteres array (param) and background model (Bmodel). The required
%   parameters can also be found in the (info) structure.
%
%
% PARAMETERS
% name     symbol  default lower bound upper bound
% -----------------------------------------------------------------------
% PARAM(1)  lam0    0.1       0            1        Unmodulated pathway amplitude
% PARAM(2)  lam1    0.8       0            1        1st Modulated pathway amplitude
% PARAM(3)  lam2    0.1       0            1        2nd Modulated pathway amplitude
% PARAM(4)  T02    max(t)  max(t)-2      max(t)+2   2nd Modulated pathway refocusing time
% -----------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function varargout = exp_ovl4pdeer(t,r,param,Bmodel)

nParam = 4;

if nargin~=1 && nargin~=4
    error('Model requires one or four input arguments.')
end

if nargin==1
    % If no inputs given, return info about the parametric model
    info.model  = '4-pulse DEER experiment (single pathway)';
    info.nparam  = nParam;
    info.parameters(1).name = 'Unmodulated pathway amplitude';
    info.parameters(1).range = [0 1];
    info.parameters(1).default = 0.4;
    info.parameters(1).units = '';
    
    info.parameters(2).name = '1st Modulated pathway amplitude';
    info.parameters(2).range = [0 1];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = '';
    
    info.parameters(3).name = '2nd Modulated pathway amplitude';
    info.parameters(3).range = [0 1];
    info.parameters(3).default = 0.1;
    info.parameters(3).units = '';
    
    info.parameters(4).name = '2nd Modulated pathway refocusing time';
    info.parameters(4).range = [max(t) - 2 max(t) + 2];
    info.parameters(4).default = max(t);
    info.parameters(4).units = '';
    varargout{1} = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Extract parameter
lambda = param(1:3);
lambda = lambda/sum(lambda);
T0 = [0 param(4)];

% Dipolar pathways
pathway(1,:) = [lambda(1) NaN];
pathway(2,:) = [lambda(2) T0(1)];
pathway(3,:) = [lambda(3) T0(2)];

% Construct multi-pathway kernel and background
K = dipolarkernel(t,r,pathway,Bmodel);
B = Bmodel(lambda(2)*(t-T0(1))).*Bmodel(lambda(3)*(t-T0(2)));

varargout{1} = K;
varargout{2} = B;

end