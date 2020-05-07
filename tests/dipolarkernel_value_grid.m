function [pass,maxerr] = test(opt)

% Test whether kernel matrix element (calculated using grid average) is correct.

% Generate kernel numerically
t = 1; % us
r = 1; % nm
K = dipolarkernel(t,r,'Method','grid','nKnots',2e6);

% Kernel value for 1us and 1nm computed using Mathematica (FresnelC and
% FresnelS) and CODATA 2018 values for ge, muB, mu0, and h.
Kref = 0.024697819895260188;

% Pass: numerical kernel is accurately computed
pass = abs(K-Kref) < 1e-6;

maxerr = max(abs(K-Kref));
 

end
