function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

N = 200;
dt = 0.008;
t = linspace(0,dt*N,N);
DistAxis = time2dist(t);
kernelF = dipolarkernel(t,DistAxis,'KCalcMethod','fresnel');
kernelE = dipolarkernel(t,DistAxis,'KCalcMethod','explicit');

err = any(abs(kernelF - kernelE)>1e-2);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end