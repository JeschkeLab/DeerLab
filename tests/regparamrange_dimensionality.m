function [pass,maxerr] = test(opt)

% Check indifference of winlowpass() towards input dimensionality

t = linspace(0,5,80);
r = time2dist(t);
K = dipolarkernel(t,r);
L = regoperator(r,2);

alphas = regparamrange(K,sparse(L));

% Pass: the alphas vector is a column vector
pass = iscolumn(alphas);

maxerr = NaN;
 

end