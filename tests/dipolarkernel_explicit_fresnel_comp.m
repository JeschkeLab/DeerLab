function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

N = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*N,N);
DistAxis = time2dist(TimeAxis);
kernelF = dipolarkernel(TimeAxis,DistAxis,'KernelCalcMethod','fresnel');
kernelE = dipolarkernel(TimeAxis,DistAxis,'KernelCalcMethod','explicit');

err = any(abs(kernelF - kernelE)>1e-2);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end