%
% EX_4PDEER Single-pathway 4-pulse DEER experiment model 
%
%   info = EX_4PDEER(t)
%   Returns an (info) structure containing the specifics of the model, including
%   a list of parameters.
%
%   pathways = EX_4PDEER(t,param)
%   Computes the dipolar pathway information array according to the paramater
%   array (param).
%
%
% PARAMETERS
% name     symbol default lower bound upper bound
% -----------------------------------------------------------------------
% PARAM(1)  lam     0.3       0            1    modulated pathway amplitude (modulation depth)
% -----------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = ex_4pdeer(t,param)

nParam = 1;

if nargin>2
    error('Model requires one or two input arguments.')
end

if nargin==1
    % If no inputs given, return info about the parametric model
    info.model  = '4-pulse DEER experiment (single pathway)';
    info.nparam  = nParam;
    info.parameters(1).name = 'modulation depth';
    info.parameters(1).range = [0 1];
    info.parameters(1).default = 0.3;
    info.parameters(1).units = '';
    
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Dipolar pathways
lambda = param(1);
pathway(1,:) = [1-lambda NaN];
pathway(2,:) = [lambda 0];
output = pathway;

end
