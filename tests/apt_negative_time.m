function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(-0.6,2,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

DistDomainSmoothing = 0.05;
%Test apt using a 1GHz excitation bandwidth
aptKernel = aptkernel(TimeAxis,'ExcitationBandwidth',1000);
[aptDistribution,aptDistanceAxis] = apt(DipEvoFcn,aptKernel,DistDomainSmoothing);

error = abs(aptDistribution - Distribution);
err = any(error>9e-1);
data = [];
maxerr = max(error);

if opt.Display
    figure(8),clf
    subplot(121)
    hold on
    plot(TimeAxis,DipEvoFcn)
    plot(TimeAxis,dipolarkernel(TimeAxis,aptDistanceAxis)*aptDistribution)
    subplot(122)
    hold on
    plot(DistanceAxis,Distribution)
    plot(aptDistanceAxis,aptDistribution)
    
end

end