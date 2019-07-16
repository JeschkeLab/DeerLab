function [err,data] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 150;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2distAxis(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

DistDomainSmoothing = 0.2;
%Test APT using a 1GHz excitation bandwidth
APTKernel = getAPTkernel(TimeAxis,'ExcitationBandwidth',1000);
outDist = APT(DipEvoFcn,APTKernel,DistDomainSmoothing);

err = any(abs(outDist.Distribution - Distribution)>1e-1);
data = [];

if opt.Display
outDist.plot(DistanceAxis,Distribution)
end

end