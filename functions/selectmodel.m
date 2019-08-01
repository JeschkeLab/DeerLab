function [optimum,aicc] = selectmodel(Models,Signal,DistanceAxis,Kernel)

N = length(Signal);

aicc = zeros(length(Models),1);
for i=1:length(Models)
    
    currentModel = Models{i};
    Info = currentModel();
    [FitDistribution] = fitparamodel(Signal,Kernel,DistanceAxis,currentModel,[]);
    K = Info.nParam + 1;
    aicc(i) = N*log(sum(Kernel*FitDistribution - Signal).^2/N) + 2*K + (2*K*(K+1))/(N - K - 1);
    
end

[~,optimum] = min(aicc);

