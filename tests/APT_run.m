function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 150;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

DistDomainSmoothing = 0.2;
%Test APT using a 1GHz excitation bandwidth
APTKernel = getAPTkernel(TimeAxis,'ExcitationBandwidth',1000);
[APTDistribution,DistanceAxis] = APT(DipEvoFcn,APTKernel,DistDomainSmoothing);

error = abs(APTDistribution - Distribution);
err = any(error>1e-1);
data = [];
maxerr = max(error);

if opt.Display
plot(DistanceAxis,APTDistribution)
end

end