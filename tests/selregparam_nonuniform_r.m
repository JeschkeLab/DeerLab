function [pass,maxerr] = test(opt)

% Check that selregparam works with nonuniform distance axes

rng(1)

t = linspace(0,3,200);

r = sqrt(linspace(1,7^2,200));
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
V = K*P + whitegaussnoise(t,0.02);

alpha = selregparam(V,K,r,'tikhonov','aic');
alpharef = 0.409459068846533;

% Pass: alpha2 compensates for larger condition number
pass = ~isequal(alpha,alpharef);

maxerr = NaN;
 

end