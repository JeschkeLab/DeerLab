%
% EXP_4PDEER Single-pathway 4-pulse DEER experiment model 
%
%   info = EXP_4PDEER(t)
%   Returns an (info) structure containing the specifics of the model.
%
%   [K,B] = EXP_4PDEER(t,r,param)
%   [K,B] = EXP_4PDEER(t,r,param,Bmodel)
%   Computes the NxM dipolar kernel (K) and N-point multi-pathway background (B) 
%   from the N-point time axis (t) and M-point distance axis (r) according to
%   the paramteres array (param) and background model (Bmodel). The required
%   parameters can also be found in the (info) structure.
%
%
% PARAMETERS
% name     symbol default lower bound upper bound
% -----------------------------------------------------------------------
% PARAM(1)  lam     0.3       0            1    Modulated pathway amplitude (modulation depth)
% -----------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function varargout = exp_4pdeer(t,r,param,Bmodel)

nParam = 1;

if nargin==2
    error('Model doesn''t work with 2 input arguments.')
end
if nargin>4
    error('Model takes at most 4 input arguments, not %d.',nargin)
end

if nargin<4
  Bmodel = [];
end

if nargin==0
    % If no inputs given, return info about the parametric model
    info.model  = '4-pulse DEER experiment (single pathway)';
    info.nparam  = nParam;
    info.parameters(1).name = 'Modulation depth';
    info.parameters(1).range = [0 1];
    info.parameters(1).default = 0.3;
    info.parameters(1).units = '';
    
    varargout{1} = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Parse input
validateattributes(r,{'numeric'},{'nonnegative','increasing','nonempty'},mfilename,'r')

% Dipolar pathways
lambda = param(1);
pathway(1,:) = [1-lambda NaN];
pathway(2,:) = [lambda 0];

% Construct multi-pathway kernel and background
if ~isempty(Bmodel)
    K = dipolarkernel(t,r,pathway,Bmodel);
    B = Bmodel(lambda*t);
else
    K = dipolarkernel(t,r,pathway);
    B = ones(size(t));
end
varargout{1} = K;
varargout{2} = B;

end
