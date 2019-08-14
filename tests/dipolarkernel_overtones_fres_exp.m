
function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

N = 100;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*N,N);
DistAxis = time2dist(TimeAxis);
Tmix = 50; %us
T1 = 88; %us

coefficients = overtones(3,Tmix,T1);

kernelF = dipolarkernel(TimeAxis,DistAxis,[],'KernelCalcMethod','fresnel','OvertoneCoeffs',coefficients);
kernelE = dipolarkernel(TimeAxis,DistAxis,[],'KernelCalcMethod','explicit','OvertoneCoeffs',coefficients);
kernelF = kernelF/mean(diff(DistAxis));
kernelE = kernelE/mean(diff(DistAxis));

err = any(abs(kernelF - kernelE)>1e-3);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end