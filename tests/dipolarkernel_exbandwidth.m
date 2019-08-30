function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,Dimension*TimeStep,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = rd_onegaussian(DistanceAxis,[3,0.5]);
Distribution = Distribution/sum(Distribution);
ExcitationBandwidth = 0.05; %MHz

Kernel = dipolarkernel(TimeAxis,DistanceAxis,[],'KernelBType','full','ExcitationBandwidth',ExcitationBandwidth);

Trace  = Kernel*Distribution;

err = any(abs(Trace - 0*Trace)>1e-5);
maxerr = max(abs(Trace -  0*Trace));
data = [];
end