function [pass,maxerr] = test(opt)

% Check non-square kernel matrix construction

Ntime = 50;
Ndist = 70;
t = linspace(-0.5,2,Ntime);
r = linspace(1,5,Ndist);
KF = dipolarkernel(t,r,'Method','fresnel');
KG = dipolarkernel(t,r,'Method','grid');
KI = dipolarkernel(t,r,'Method','integral');

% Pass 1: fresnel-constructed kernel has right dimensions
pass(1) = all(size(KF) == [Ntime Ndist]);
% Pass 2: fresnel-constructed kernel has right dimensions
pass(2) = all(size(KG) == [Ntime Ndist]);
% Pass 2: integral-constructed kernel has right dimensions
pass(3) = all(size(KI) == [Ntime Ndist]);

pass = all(pass);

maxerr = NaN;
 

end