function [pass,maxerr] = test(opt)

% Check that one can generate kernels with one single distance-domain point

r = linspace(1,5,100);
K = dipolarkernel(0.5,r);

% Pass: time-domain dimension is a singlet
pass = size(K,1)==1;

maxerr = NaN;
 

end