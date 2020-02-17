function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check whether kernel matrix has correct size
%======================================================

Ntime = 50;
Ndist = 70;

t = linspace(-0.5,2,Ntime);
r = linspace(1,5,Ndist);

KF = dipolarkernel(t,r,'Method','fresnel');
KG = dipolarkernel(t,r,'Method','grid');

err(1) = any(size(KF) ~= [Ntime Ndist]);
err(2) = any(size(KG) ~= [Ntime Ndist]);
err = any(err);

maxerr = 0;
data = [];

end