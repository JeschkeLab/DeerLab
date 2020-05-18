function [pass,maxerr] = test(opt)

% Check that dipolarkernel works with non-uniform distance axes

rng(1)

t = linspace(0,4,400);
r = sqrt(linspace(1,7^2,200));

P = dd_gauss(r,[4.5 0.5]);

K = dipolarkernel(t,r);

V = K*P;

% Pass 1: normalization works allright
pass(1) = abs(V(t==0) - 1) < 1e-6;

pass = all(pass);

maxerr = max(V(t==0) - 1);


end