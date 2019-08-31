function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

load(fullfile('comparison','oldDAkernel500'));

N = 500;
dt = 0.008;
t = linspace(0,dt*N,N);
DistAxis = time2dist(t);
kernelOut = dipolarkernel(t,DistAxis);
kernelOut = kernelOut/mean(diff(DistAxis));

err = any(abs(kernelOut - kernel)>7e-3);
maxerr = max(max(abs(kernelOut - kernel)));
data = [];

end