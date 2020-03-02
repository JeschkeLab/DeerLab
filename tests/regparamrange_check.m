function [pass,maxerr] = test(opt)

% Check that regparamrange resolution can be set up manually

t = linspace(0,3,200);
r = time2dist(t);
K = dipolarkernel(t,r);
L = regoperator(r,2);
lgRes = 0.1;

alpha = regparamrange(K,L,'Resolution',lgRes);

% Pass: the right number of elements are returned
pass = length(alpha) == 85;

maxerr = NaN;
 

end