function [err,data] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Dimension = 200;
TimeStep = 0.008;
rmin = (4*TimeStep*52.04/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);
TimeAxis = linspace(0,TimeStep,Dimension);
DistanceAxis = linspace(rmin,rmax,Dimension);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution'/sum(Distribution);
Background = exp(-150*TimeAxis);
Kernel = getKernel(Dimension,TimeStep);

Trace  = Kernel*Distribution;
Trace = (Trace + 10).*Background';

Background = Background/Trace(1);
Trace = Trace/Trace(1);

KernelB = getKernel(Dimension,TimeStep,[],[],Background);

TraceB  = KernelB*Distribution;

err = any(abs(TraceB - Trace)>1e-10);
data = [];

end