function [pass,maxerr] = test(opt)

% Check that selregparam works with non-square kernel matrices

t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

alpha1 = selregparam(S,K,r,'tikhonov','aic');

t = linspace(0,3,400);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

alpha2 = selregparam(S,K,r,'tikhonov','aic');

% Pass: alpha2 compensates for larger condition number
pass = alpha2 > alpha1;

maxerr = NaN;
 

end