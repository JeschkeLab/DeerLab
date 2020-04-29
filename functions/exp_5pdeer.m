%
% EXP_5PDEER 5-pulse DEER experiment model 
%
%   info = EXP_5PDEER(t)
%   Returns an (info) structure containing the specifics of the model.
%
%   pathways = EXP_5PDEER(t,param)
%   Computes the dipolar pathway informatio matrix according to the paramters
%   array (param) and the specified experiment. The required parameters can
%   also be found in the (info) structure.
%
%
% PARAMETERS
% name     symbol  default lower bound upper bound
% -----------------------------------------------------------------------
% PARAM(1)  lam0    0.4       0            1        Unmodulated pathway amplitude
% PARAM(2)  lam1    0.4       0            1        1st Modulated pathway amplitude
% PARAM(3)  lam2    0.2       0            1        2nd Modulated pathway amplitude
% PARAM(4)  T02  max(t)/2  max(t)/2-2  max(t)/2+2   2nd Modulated pathway refocusing time
% -----------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = exp_5pdeer(t,param)

nParam = 4;

if nargin==1 && nargin>2
    error('Model requires at two input arguments.')
end

if nargin==1
    % If no inputs given, return info about the parametric model
    info.model  = '5-pulse DEER experiment (two modulated pathways)';
    info.nparam  = nParam;
    info.parameters(1).name = 'Unmodulated pathway amplitude';
    info.parameters(1).range = [0 1];
    info.parameters(1).default = 0.4;
    info.parameters(1).units = '';
    
    info.parameters(2).name = '1st Modulated pathway amplitude';
    info.parameters(2).range = [0 1];
    info.parameters(2).default = 0.4;
    info.parameters(2).units = '';
    
    info.parameters(3).name = '2nd Modulated pathway amplitude';
    info.parameters(3).range = [0 1];
    info.parameters(3).default = 0.2;
    info.parameters(3).units = '';
    
    info.parameters(4).name = '2nd Modulated pathway refocusing time';
    info.parameters(4).range = [max(t)/2 - 2 max(t)/2 + 2];
    info.parameters(4).default = max(t)/2;
    info.parameters(4).units = '';
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Extract parameter
lambda = param(1:3);
lambda = lambda/sum(lambda);
T0 = [0 param(4)];

% Dipolar pathways
pathways(1,:) = [lambda(1) NaN];
pathways(2,:) = [lambda(2) T0(1)];
pathways(3,:) = [lambda(3) T0(2)];
output = pathways;

end
