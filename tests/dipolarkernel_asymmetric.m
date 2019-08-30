function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Ntime = 100;
Ndist = 200;

TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Ntime,Ntime);
[~,rmin,rmax] = time2dist(TimeAxis);
DistAxis = linspace(rmin,rmax,Ndist);
kernelF = dipolarkernel(TimeAxis,DistAxis,'KernelCalcMethod','fresnel');
kernelE = dipolarkernel(TimeAxis,DistAxis,'KernelCalcMethod','explicit');
kernelF = kernelF/mean(diff(DistAxis));
kernelE = kernelE/mean(diff(DistAxis));
err(1) = any(any(abs(kernelF - kernelE)>1e-2));
err(2) = any(size(kernelF) ~= [Ntime Ndist]);
err = any(err);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end