function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

N = 200;
dt = 0.008;
t = linspace(0,dt*N,N);
DistAxis = time2dist(t);

KF = dipolarkernel(t,DistAxis,'Method','fresnel');
KG = dipolarkernel(t,DistAxis,'Method','grid');

delta = abs(KF - KG);
err = any(delta(:)>1e-2);
maxerr = max(max(delta));

data = [];

end