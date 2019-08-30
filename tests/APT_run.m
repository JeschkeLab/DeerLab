function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

DistDomainSmoothing = 0.1;
%Test apt using a 1GHz excitation bandwidth
aptKernel = aptkernel(TimeAxis,'ExcitationBandwidth',1000);
[aptDistribution,aptDistanceAxis] = apt(DipEvoFcn,aptKernel,DistDomainSmoothing);

error = abs(aptDistribution - Distribution);
err = any(error>9e-1);
data = [];
maxerr = max(error);

if opt.Display
    figure(8),clf
    hold on
    plot(DistanceAxis,Distribution)
    plot(aptDistanceAxis,aptDistribution)
    
end

end