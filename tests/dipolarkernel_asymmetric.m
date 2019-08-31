function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Ntime = 100;
Ndist = 200;

dt = 0.008;
t = linspace(0,dt*Ntime,Ntime);
[~,rmin,rmax] = time2dist(t);
DistAxis = linspace(rmin,rmax,Ndist);
kernelF = dipolarkernel(t,DistAxis,'KCalcMethod','fresnel');
kernelE = dipolarkernel(t,DistAxis,'KCalcMethod','explicit');
kernelF = kernelF/mean(diff(DistAxis));
kernelE = kernelE/mean(diff(DistAxis));
err(1) = any(any(abs(kernelF - kernelE)>1e-2));
err(2) = any(size(kernelF) ~= [Ntime Ndist]);
err = any(err);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end