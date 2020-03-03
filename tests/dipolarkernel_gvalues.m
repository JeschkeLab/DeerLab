function [pass,maxerr] = test(opt)

% Check whether K matrix elements scale properly with t and r

t = 2; % nm
r = 3; % us

% isotropic g values of two spins
g1 = 2.1;
g2 = 2.4;

K1 = dipolarkernel(t,r,'g',[g1 g2]);
K2 = dipolarkernel(t,r/(g1*g2/gfree^2)^(1/3));

maxerr = max(abs(K1-K2));

% Pass: kernel is properly scaled
pass = maxerr < 1e-15;
 
end