function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,Dimension*TimeStep,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);
Background = exp(-0.05*TimeAxis);
ExcitationBandwidth = 10; %MHz

Kernel = getKernel(TimeAxis,DistanceAxis,[],'ExcitationBandwidth',ExcitationBandwidth);

Trace  = Kernel*Distribution;
Trace = (Trace + 2).*Background';

Background = Background*(1-1/Trace(1));
Trace = Trace/Trace(1);

KernelB = getKernel(TimeAxis,DistanceAxis,Background,'KernelBType','full','ExcitationBandwidth',ExcitationBandwidth);

TraceB  = KernelB*Distribution;

err = any(abs(TraceB - Trace)>1e-10);
maxerr = max(abs(TraceB - Trace));
data = [];

end