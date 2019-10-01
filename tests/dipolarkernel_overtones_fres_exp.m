
function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

N = 100;
t = linspace(0,3,N);
r = time2dist(t);
Tmix = 50; %us
T1 = 88; %us

coefficients = overtones(3,Tmix,T1);

kernelF = dipolarkernel(t,r,'Method','fresnel','OvertoneCoeffs',coefficients);
kernelE = dipolarkernel(t,r,'Method','explicit','OvertoneCoeffs',coefficients);
kernelF = kernelF/mean(diff(r));
kernelE = kernelE/mean(diff(r));

err = any(abs(kernelF - kernelE)>1e-3);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end