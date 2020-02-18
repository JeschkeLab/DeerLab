function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Compare kernel calculation methods
%======================================================

N = 50;
t = linspace(0,3,N);
r = linspace(1,5,N);

KF = dipolarkernel(t,r,'Method','fresnel');
KG = dipolarkernel(t,r,'Method','grid');

delta = abs(KF-KG);
err = any(delta(:)>1e-4);

maxerr = max(delta(:));
data = [];

end