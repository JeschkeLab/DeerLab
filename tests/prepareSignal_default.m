function [err,data] = test(opt,olddata)

%======================================================
% Test data preparation function
%======================================================

%Parameters
Dimension = 200;
TimeStep = 0.008;
rmin = (4*TimeStep*52.04/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
Length = length(TimeAxis);

DistanceAxis = linspace(rmin,rmax,Dimension);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution'/sum(Distribution);
Background = exp(-0.05*TimeAxis);
Kernel = getKernel(Dimension,TimeStep/1000);

DipEvoFcn = Kernel*Distribution;
DipEvoFcn = DipEvoFcn';
ExpData = (DipEvoFcn + 10).*Background;

Signal = pdsdata('TimeAxis',TimeAxis,'ExpData',ExpData);
Signal1 = Signal.prepare;
Signal2 = Signal.prepare();
Signal3 = prepare(Signal);

err(1) = Signal.TimeStep ~= TimeStep*1000;
err(2) = any(abs(Signal1.DipEvoFcn - DipEvoFcn)>1e-2);
err(3) = any(abs(Signal2.DipEvoFcn - DipEvoFcn)>1e-2);
err(4) = any(abs(Signal3.DipEvoFcn - DipEvoFcn)>1e-2);

err = any(err);
data = [];

end