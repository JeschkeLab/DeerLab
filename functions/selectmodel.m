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

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({Models,Signal,DistanceAxis,Kernel,Methods,varargin});
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
