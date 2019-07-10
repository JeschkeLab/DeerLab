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
Distribution = Distribution/sum(Distribution);
Background = exp(-0.05*TimeAxis);
Kernel = getKernel(Dimension,TimeStep*1000,rmin,rmax);

DipEvoFcn = Kernel*Distribution;
DipEvoFcn = DipEvoFcn';
ExpData = (DipEvoFcn + 1).*Background;

Signal = DAsignal('TimeAxis',TimeAxis,'ExpData',ExpData);
Signal = prepare(Signal);

Opts = DAoptions('DistDomainSmoothing',0.2);

outDist = APT(Signal,Opts);
Distribution = Distribution/sum(Distribution);

err(1) = any(abs(outDist.Distribution - Distribution)>1e-1);
err(2) = any(abs(outDist.FitSignal - Signal.DipEvoFcn')>1e-1);
err(3) = ~strcmp(outDist.Method,'APT');

err = any(err);
data = [];
if opt.Display
outDist.plot(DistanceAxis,Distribution)
end

end