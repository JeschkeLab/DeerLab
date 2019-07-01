function [err,data] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 150;
TimeStep = 0.008;
rmin = (4*TimeStep*52.04/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
Length = length(TimeAxis);

DistanceAxis = linspace(rmin,rmax,Dimension);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution'/sum(Distribution);
Background = exp(-0.05*TimeAxis);
Kernel = getAPTkernel(Dimension,TimeStep);

DipEvoFcn = Kernel*Distribution;
ExpData = (DipEvoFcn + 10).*Background;

Signal = pdsdata('TimeAxis',TimeAxis,'ExpData',ExpData);
Signal = prepare(Signal);

Opts = daopts('DistDomainSmoothing',0.05);

outDist = APT(Signal,Opts);
outDist.plotDistr;
err = any(abs(outDist.Distribution - Distribution)>1e-2);

data = [];

end