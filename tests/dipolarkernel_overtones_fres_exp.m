
function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

N = 100;
dt = 0.008;
t = linspace(0,dt*N,N);
DistAxis = time2dist(t);
Tmix = 50; %us
T1 = 88; %us

coefficients = overtones(3,Tmix,T1);

kernelF = dipolarkernel(t,DistAxis,'Method','fresnel','OvertoneCoeffs',coefficients);
kernelE = dipolarkernel(t,DistAxis,'Method','explicit','OvertoneCoeffs',coefficients);
kernelF = kernelF/mean(diff(DistAxis));
kernelE = kernelE/mean(diff(DistAxis));

err = any(abs(kernelF - kernelE)>1e-3);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end