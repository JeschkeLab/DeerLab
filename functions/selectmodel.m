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
%   [opt,f] = SELECTMODEL(...)
%   Returns the method selector functionals for the different methods.
%
%   alphas = SELECTMODEL(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%   See "help fitparamodel" for a detailed list of the property-value pairs
%   accepted by the function.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function [optima,functionals] = selectmodel(Models,Signal,DistanceAxis,Kernel,Methods,varargin)

%Input validation
if ~iscell(Methods)
   Methods = {Methods}; 
end
if ~iscolumn(DistanceAxis)
    DistanceAxis = DistanceAxis.';
end
if length(varargin)==1
   varargin = varargin{1}; 
end
allowedMethodInputs = {'aic','aicc','bic'};
if iscell(Methods)
    for i=1:length(Methods)
        if strcmp(Methods{i},'all')
            Methods = allowedMethodInputs;
            break;
        end
        validateattributes(Methods{i},{'char'},{'nonempty'})
        Methods{i} = validatestring(Methods{i},allowedMethodInputs);
    end
else
    validateattributes(Methods,{'char'},{'nonempty'})
    if strcmp(Methods,'all')
        Methods = allowedMethodInputs;
    else
        Methods = validatestring(Methods,allowedMethodInputs);
        Methods = {Methods};
    end
end
%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
Modelsstring = [];
for i=1:length(Models)
   Modelsstring = [Modelsstring func2str(Models{i})]; 
end 
    hashKey = datahash({Modelsstring,Signal,DistanceAxis,Kernel,Methods,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    [optima,functionals] = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Pre-allocate vectors
N = length(Signal);
aicc = zeros(length(Models),1);
bic = zeros(length(Models),1);
aic = zeros(length(Models),1);
%Run the different parametric model fits
for i=1:length(Models)
    currentModel = Models{i};
    Info = currentModel();
    [FitDistribution] = fitparamodel(Signal,Kernel,DistanceAxis,currentModel,[],varargin);
    K = Info.nParam + 1;
    aicc(i) = N*log(sum(Kernel*FitDistribution - Signal).^2/N) + 2*K + (2*K*(K+1))/(N - K - 1);
    aic(i) = N*log(sum(Kernel*FitDistribution - Signal).^2/N) + 2*K;
    bic(i) = N*log(sum(Kernel*FitDistribution - Signal).^2/N) + K*log(N);
end

%Apply the requested selection methods
optima = zeros(length(Methods),1);
functionals = cell(length(Methods),1);
for i=1:length(Methods)
    currentMethod = Methods{i};
    switch currentMethod
        case 'aic'
            functional = aic;
            [~,optimum] = min(functional);
        case 'aicc'
            functional = aicc;
            [~,optimum] = min(functional);
        case 'bic'
            functional = bic;
            [~,optimum] = min(functional);
    end
    functionals{i} = functional;
    optima(i) = optimum;
end

if length(functionals)==1
   functionals = functionals{1}; 
end

%Store output result in the cache
Output = {optima,functionals};
cachedData = addcache(cachedData,hashKey,Output);
