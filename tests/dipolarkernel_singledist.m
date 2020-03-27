function [pass,maxerr] = test(opt)

% Check that one can generate kernels with one single distance-domain point

t = linspace(0,3,100);
K = dipolarkernel(t,3);

% Pass: distance-domain dimension is a singlet
pass = size(K,2)==1;

maxerr = NaN;
 

end