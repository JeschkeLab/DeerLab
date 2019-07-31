function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

load(fullfile('comparison','oldDAkernel500'));

N = 500;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*N,N);
DistAxis = time2dist(TimeAxis);
kernelOut = dipolarkernel(TimeAxis,DistAxis);

err = any(abs(kernelOut - kernel)>1e-2);
maxerr = max(max(abs(kernelOut - kernel)));
data = [];

end