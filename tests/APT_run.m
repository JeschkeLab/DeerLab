function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 150;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussfcn(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

DistDomainSmoothing = 0.2;
%Test apt using a 1GHz excitation bandwidth
aptKernel = aptkernel(TimeAxis,'ExcitationBandwidth',1000);
[aptDistribution,DistanceAxis] = apt(DipEvoFcn,aptKernel,DistDomainSmoothing);

error = abs(aptDistribution - Distribution);
err = any(error>1e-1);
data = [];
maxerr = max(error);

if opt.Display
plot(DistanceAxis,aptDistribution)
end

end