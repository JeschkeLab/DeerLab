%
% EXP_5PDEER 7-pulse DEER experiment model 
%
%   info = EXP_7PDEER(t)
%   Returns an (info) structure containing the specifics of the model.
%
%   pathways = EXP_7PDEER(t,param)
%   Computes the dipolar pathway informatio matrix according to the paramters
%   array (param) and the specified experiment. The required parameters can
%   also be found in the (info) structure.
%
%
% PARAMETERS
% name     symbol  default lower bound upper bound
% -----------------------------------------------------------------------
% PARAM(1)  lam0    0.4       0            1            Unmodulated pathway amplitude
% PARAM(2)  lam1    0.4       0            1            1st Modulated pathway amplitude
% PARAM(3)  lam2    0.2       0            1            2nd Modulated pathway amplitude
% PARAM(4)  lam3    0.2       0            1            2nd Modulated pathway amplitude
% PARAM(5)  T02  max(t)/5   max(t)/5-2     max(t)/5+2   2nd Modulated pathway refocusing time
% PARAM(6)  T03  max(t)*2/5  max(t)*2/5-2  max(t)2/5+2  3rd Modulated pathway refocusing time
% -----------------------------------------------------------------------
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function output = exp_7pdeer(t,param)

nParam = 6;

if nargin==1 && nargin>2
    error('Model requires at two input arguments.')
end

if nargin==1
    % If no inputs given, return info about the parametric model
    info.model  = '7-pulse DEER experiment (three modulated pathways)';
    info.nparam  = nParam;
    info.parameters(1).name = 'Unmodulated pathway amplitude';
    info.parameters(1).range = [0 1];
    info.parameters(1).default = 0.3;
    info.parameters(1).units = '';
    
    info.parameters(2).name = '1st Modulated pathway amplitude';
    info.parameters(2).range = [0 1];
    info.parameters(2).default = 0.5;
    info.parameters(2).units = '';
    
    info.parameters(3).name = '2nd Modulated pathway amplitude';
    info.parameters(3).range = [0 1];
    info.parameters(3).default = 0.3;
    info.parameters(3).units = '';
    
    info.parameters(4).name = '3rd Modulated pathway amplitude';
    info.parameters(4).range = [0 1];
    info.parameters(4).default = 0.2;
    info.parameters(4).units = '';
    
    info.parameters(5).name = '2nd Modulated pathway refocusing time';
    info.parameters(5).range = [max(t)/5 - 2, max(t)/5 + 2];
    info.parameters(5).default = max(t)/5;
    info.parameters(5).units = 'us';
    
    info.parameters(6).name = '3rd Modulated pathway refocusing time';
    info.parameters(6).range = [max(t)*2/5 - 2, max(t)*2/5 + 2];
    info.parameters(6).default = max(t)*2/5;
    info.parameters(6).units = 'us';
    
    output = info;
    return
end

% Assert that the number of parameters matches the model
if length(param)~=nParam
    error('The number of input parameters does not match the number of model parameters.')
end

% Extract parameters
lambda = param(1:4);
lambda = lambda/sum(lambda);
T0 = [0 param(5:6)];

% Dipolar pathways
pathways(1,:) = [lambda(1) NaN];
pathways(2,:) = [lambda(2) T0(1)];
pathways(3,:) = [lambda(3) T0(2)];
pathways(4,:) = [lambda(4) T0(3)];
output = pathways;

end
