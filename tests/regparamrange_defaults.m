function [pass,maxerr] = test(opt)

% Check that regparamrange works with the defaults

t = linspace(0,3,200);
r = time2dist(t);
K = dipolarkernel(t,r);
L = regoperator(r,2);

alpha = regparamrange(K,sparse(L));

% Pass: the right number of elements are returned
pass = length(alpha) == 84;

maxerr = NaN;
 

end