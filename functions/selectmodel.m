%
% SELECTMODEL Optimal parametric model selection
%
%   opt = SELECTMODEL({@model1,...,@modelN},S,r,K,{'aic',...})
%   Evaluates the fits of the parametric models (model1,...,modelN) to a
%   signal (S) according to the dipolar kernel (K) and distance axis (r).
%   The models must be passed as a cell array of function handles. Each fit
%   is then evaluated according to the model selection criterions 
%   ('aic','aicc','bic') specified in the last input argument M-point cell array. 
%   Function returns a M-point array containing the optimal models
%   according to each selection method.
%
%   [opt,f,params] = SELECTMODEL(...)
%   Returns a cell array the method selector functionals (f) for the
%   different methods and a cell array (params) with the fitted parameters
%   for each of the evaluated models.
%
%   opt = SELECTMODEL(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%   See "help fitparamodel" for a detailed list of the property-value pairs
%   accepted by the function.
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [optima,functionals,fitparams] = selectmodel(Models,S,r,K,Methods,varargin)

% Input validation
%-------------------------------------------------------------------------------
if ~iscell(Models)
  error('First input must be a cell array of model functions.');
end
if ~iscell(Methods)
   Methods = {Methods}; 
end
if length(varargin)==1
   varargin = varargin{1}; 
end

allowedMethodInputs = {'aic','aicc','bic'};
for i = 1:length(Methods)
    if strcmp(Methods{i},'all')
        Methods = allowedMethodInputs;
        break;
    end
    validateattributes(Methods{i},{'char'},{'nonempty'})
    Methods{i} = validatestring(Methods{i},allowedMethodInputs);
end

% Convert distance axis to nanometers if given in Angstrom
if ~isnanometer(r)
   r = r/10; % A -> nm
end

% Run all parametric model fits and evaluate selection metrics
%-------------------------------------------------------------------------------
nMethods = length(Methods);
N = numel(S);
AICc = zeros(nMethods,1);
BIC = zeros(nMethods,1);
AIC = zeros(nMethods,1);
for i = 1:length(Models)
    [paramfit,Pfit] = fitparamodel(S,Models{i},r,K,varargin);
    fitparams{i} = paramfit;
    
    nParams = numel(paramfit);
    Q = nParams + 1;
    SSR = sum((S(:)-K*Pfit(:)).^2);
    AIC(i) =  N*log(SSR/N) + 2*Q;
    AICc(i) = N*log(SSR/N) + 2*Q + 2*Q*(Q+1)/(N-Q-1);
    BIC(i) =  N*log(SSR/N) + Q*log(N);
end

% Identify optimal models based on selection criteria
%-------------------------------------------------------------------------------
optima = zeros(nMethods,1);
functionals = cell(nMethods,1);
for i = 1:nMethods
    switch Methods{i}
        case 'aic'
            functional = AIC;
        case 'aicc'
            functional = AICc;
        case 'bic'
            functional = BIC;
    end
    [~,optimum] = min(functional);
    functionals{i} = functional;
    optima(i) = optimum;
end

if length(functionals)==1
   functionals = functionals{1}; 
end
