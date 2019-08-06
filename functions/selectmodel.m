function [optima,aicc,bic] = selectmodel(Models,Signal,DistanceAxis,Kernel,Methods,varargin)

if ~iscell(Methods)
   Methods = {Methods}; 
end
if ~iscolumn(DistanceAxis)
    DistanceAxis = DistanceAxis.';
end
if length(varargin)==1
   varargin = varargin{1}; 
end
N = length(Signal);
aicc = zeros(length(Models),1);
for i=1:length(Models)
    currentModel = Models{i};
    Info = currentModel();
    [FitDistribution] = fitparamodel(Signal,Kernel,DistanceAxis,currentModel,[],varargin);
    K = Info.nParam + 1;
    aicc(i) = N*log(sum(Kernel*FitDistribution - Signal).^2/N) + 2*K + (2*K*(K+1))/(N - K - 1);
    bic(i) = N*log(sum(Kernel*FitDistribution - Signal).^2/N) + K*log(N);
    
end

optima = zeros(length(Methods),1);
for i=1:length(Methods)
    currentMethod = Methods{i};
    switch currentMethod
        case 'aicc'
            [~,optimum] = min(aicc);
        case 'bic'
            [~,optimum] = min(bic);
    end
    optima(i) = optimum;
end
