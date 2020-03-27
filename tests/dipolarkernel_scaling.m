function [pass,maxerr] = test(opt)

% Check whether K matrix elements scale properly with t and r

t = 2; % nm
r = 3; % us
c = 1.2; % distance scaling factor

K1 = dipolarkernel(t,r);
K2 = dipolarkernel(t*c^3,r*c);

maxerr = max(abs(K1-K2));

% Pass: kernel is properly scaled
pass = maxerr < 1e-15;
 
end