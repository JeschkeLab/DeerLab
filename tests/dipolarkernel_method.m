function [pass,maxerr] = test(opt)

% Check that the different kernel calculation methods lead to same result

t = linspace(0,3,50);
r = linspace(1,5,50);

KF = dipolarkernel(t,r,'Method','fresnel');
KG = dipolarkernel(t,r,'Method','grid');
KI = dipolarkernel(t,r,'Method','integral');

delta1 = abs(KF-KG);
delta2 = abs(KF-KI);
delta3 = abs(KI-KG);

% Pass 1: fresnel and integral methods give same result
pass(1) = all(delta1(:) < 1e-4);
% Pass 1: fresnel and grid methods give same result
pass(2) = all(delta2(:) < 1e-4);
% Pass 1: grid and integral methods give same result
pass(3) = all(delta3(:) < 1e-4);

pass = all(pass);

maxerr = max(delta1(:));
 

end